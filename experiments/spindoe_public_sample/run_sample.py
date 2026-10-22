"""Run the unchanged upstream estimator on its public sample, on CPU."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ["MPLBACKEND"] = "Agg"
os.environ["MPLCONFIGDIR"] = str(ROOT / ".cache" / "matplotlib")
os.environ["WANDB_MODE"] = "disabled"
os.environ["WANDB_SILENT"] = "true"
manifest = json.loads((ROOT / "source_manifest.json").read_text())
for entry in manifest["files"]:
    data = (ROOT / "upstream" / entry["path"]).read_bytes()
    digest = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
    if digest != entry["git_blob_sha1"]:
        raise ValueError(f"Changed upstream file: {entry['path']}")

sys.path.insert(0, str(ROOT / "upstream" / "python"))
# Avoid unused training/data-loader dependencies. Original files remain intact;
# only these two imports are removed from transparent runtime copies.
runtime = ROOT / ".cache" / "inference_modules"
runtime.mkdir(parents=True, exist_ok=True)
for filename, unused_import in {"doe.py": "from dotdataset import read_ball_images\n", "dotdetector.py": "import wandb\n"}.items():
    source = (ROOT / "upstream/python" / filename).read_text()
    assert source.count(unused_import) == 1
    (runtime / filename).write_text(source.replace(unused_import, ""))
sys.path.insert(0, str(runtime))
os.chdir(ROOT / "upstream" / "python")
import cv2
import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial.transform import Rotation as R
import torch
from doe import DOE
from dotdetector import DotDetector
from spin_regressor import SpinRegressor
from spindoe import SpinDOE

torch.set_num_threads(4)
np.random.seed(0)
torch.manual_seed(0)

class CPUDOE(DOE):
    def init_dot_detector(self, model_path):
        # Same network and weights; bypass Lightning's GPU/training-state loader.
        checkpoint = torch.load(model_path, map_location="cpu", weights_only=True)
        self.dot_detector = DotDetector(**checkpoint["hyper_parameters"])
        self.dot_detector.load_state_dict(checkpoint["state_dict"], strict=True)
        self.dot_detector.eval()

# Bypass only the upstream constructor's hard-coded CUDA selection and GUI CLI.
estimator = SpinDOE.__new__(SpinDOE)
estimator.doe = CPUDOE(ROOT / "upstream/data/model/spindoe.ckpt", use_gpu=False)
estimator.spin_regressor = SpinRegressor()
paths = sorted((ROOT / "upstream/data/test").glob("*.png"), key=lambda p: int(p.stem))
timestamps_ns = np.array([int(p.stem) for p in paths], dtype=np.int64)
times = timestamps_ns * 1e-9  # Preserve the upstream timestamp convention.
relative_times = (timestamps_ns - timestamps_ns[0]) * 1e-9
imgs = [np.array(Image.open(p).convert("RGB")) for p in paths]
start = time.perf_counter()
with torch.inference_mode():
    spin, rots, heatmaps, valid_idx = estimator.estimate(times, imgs)
elapsed = time.perf_counter() - start
baseline = {"spin_available": spin is not None, "doe_pose_count": sum(r is not None for r in rots), "pose_frames": [i for i, r in enumerate(rots) if r is not None], "estimated_rps": None if spin is None else float(np.linalg.norm(spin)/(2*np.pi))}
np.savez_compressed(OUT / "baseline_estimates.npz", timestamps_ns=timestamps_ns, heatmaps=np.stack(heatmaps), quaternions_xyzw=np.array([r.as_quat() if r is not None else [np.nan]*4 for r in rots]))
diagnostic_ordering = False
if spin is None:
    diagnostics = {"pose_frames": [i for i, r in enumerate(rots) if r is not None], "frames": [{"min": float(h.min()), "max": float(h.max()), "dots": len(estimator.doe.extract_dots(h)[1])} for h in heatmaps]}
    (OUT / "failure_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    print(json.dumps(diagnostics, indent=2), flush=True)
    # Explicit diagnostic only, not a successful unmodified reproduction.
    # Geometry uses the first two detections as a basis: test larger blobs first.
    original_extract = estimator.doe.extract_dots
    def size_sorted_extract(heatmap):
        mask, keypoints = original_extract(heatmap)
        return mask, sorted(keypoints, key=lambda point: -point.size)
    estimator.doe.extract_dots = size_sorted_extract
    np.random.seed(0)
    with torch.inference_mode():
        spin, rots, heatmaps, valid_idx = estimator.estimate(times, imgs)
    diagnostic_ordering = True
    if spin is None:
        raise RuntimeError("Both the upstream baseline and size-sorted diagnostic produced no spin estimate")
valid_idx = np.asarray(valid_idx, dtype=int)
pose_idx = np.array([i for i, r in enumerate(rots) if r is not None])
anchor = int(valid_idx[0])
pred_raw = estimator.predicted_rots(rots[anchor], times[anchor], times, spin)
pred_reverse = estimator.predicted_rots(rots[anchor], times[anchor], times, -spin)

def errors(predictions, indexes):
    return np.array([np.degrees((predictions[i] * rots[i].inv()).magnitude()) for i in indexes])

raw_err = errors(pred_raw, valid_idx)
reverse_err = errors(pred_reverse, valid_idx)
seed_results = []
for seed in range(5):
    np.random.seed(seed)
    candidate, subset = estimator.spin_regressor.RANSAC_regress(times[pose_idx], [rots[i] for i in pose_idx])
    seed_results.append({"seed": seed, "rps": None if candidate is None else float(np.linalg.norm(candidate)/(2*np.pi)), "spin_raw_rad_s": None if candidate is None else candidate.tolist(), "inliers": None if subset is None else pose_idx[subset].tolist()})

# Independent known-motion input tests the sign against the documented forward model.
synthetic_t = np.arange(25) * 0.001
known_spin = np.array([20.0, 30.0, -40.0])
initial_rot = R.from_euler("xyz", [0.3, -0.2, 0.5])
known_rots = R.from_rotvec(synthetic_t[:, None] * known_spin) * initial_rot
np.random.seed(0)
synthetic_estimate, synthetic_inliers = SpinRegressor().RANSAC_regress(synthetic_t, known_rots)
synthetic_check = {"input_rad_s": known_spin.tolist(), "upstream_output_rad_s": synthetic_estimate.tolist(), "raw_error_rad_s": float(np.linalg.norm(synthetic_estimate-known_spin)), "sign_reversed_error_rad_s": float(np.linalg.norm(-synthetic_estimate-known_spin))}

frames = []
font = ImageFont.load_default(size=18)
small = ImageFont.load_default(size=15)
for i, img in enumerate(imgs):
    canvas = Image.new("RGB", (960, 420), (18, 22, 29))
    draw = ImageDraw.Draw(canvas)
    draw.text((20, 12), f"Size-sorted diagnostic | {np.linalg.norm(spin)/(2*np.pi):.2f} rps candidate | baseline: failed", font=font, fill="white")
    panels = [("Input", img), ("Diagnostic: original sign", estimator.plot_dots(img.copy(), pred_raw[i], (255, 70, 150))), ("Diagnostic: reversed sign", estimator.plot_dots(img.copy(), pred_reverse[i], (50, 230, 240)))]
    for col, (label, panel) in enumerate(panels):
        x = col * 320 + 22
        draw.text((x, 52), label, font=small, fill=(220, 225, 235))
        canvas.paste(Image.fromarray(panel).resize((272, 272), Image.Resampling.NEAREST), (x, 80))
    state = "RANSAC inlier" if i in valid_idx else ("pose rejected by RANSAC" if rots[i] is not None else "no DOE pose")
    draw.text((22, 370), f"Frame {i+1}/{len(imgs)} | t = {relative_times[i]*1000:.3f} ms | {state}", font=font, fill="white")
    draw.text((22, 396), "60x slow motion; actual timestamps; no independent spin ground truth.", font=small, fill=(180, 190, 205))
    frames.append(canvas)

duration = relative_times[-1] + np.median(np.diff(relative_times))
playback_times = np.arange(0, duration, 1/(30*60))
indexes = np.minimum(np.searchsorted(relative_times, playback_times, side="right") - 1, len(imgs)-1)
with imageio.get_writer(str(OUT / "diagnostic_60x.mp4"), fps=30, codec="libx264", quality=8, macro_block_size=4) as writer:
    for i in indexes:
        writer.append_data(np.array(frames[i]))
frames[8].save(OUT / "comparison_frame.png")
contact = Image.new("RGB", (960, 420*4))
for row, index in enumerate([0, 7, 14, 21]):
    contact.paste(frames[index], (0, row*420))
contact.save(OUT / "contact_sheet.png")
np.savez_compressed(OUT / "estimates.npz", timestamps_ns=timestamps_ns, heatmaps=np.stack(heatmaps), quaternions_xyzw=np.array([r.as_quat() if r is not None else [np.nan]*4 for r in rots]), ransac_indices=valid_idx, spin_raw_rad_s=spin)
summary = {
    "unmodified_algorithm_baseline": baseline,
    "diagnostic_keypoint_size_sort_applied": diagnostic_ordering,
    "source_commit": manifest["commit"], "source_files_hash_verified": True,
    "frame_count": len(paths), "image_sizes_wh": sorted(set((im.shape[1], im.shape[0]) for im in imgs)),
    "time_span_ms": float(relative_times[-1]*1000), "median_interval_ms": float(np.median(np.diff(relative_times))*1000),
    "max_interval_ms": float(np.max(np.diff(relative_times))*1000), "nominal_sampling_hz": float(1/np.median(np.diff(relative_times))),
    "doe_pose_count": int(len(pose_idx)), "ransac_inlier_count": int(len(valid_idx)), "ransac_indices": valid_idx.tolist(),
    "estimated_rps": float(np.linalg.norm(spin)/(2*np.pi)), "upstream_spin_rad_s": spin.tolist(),
    "sign_reversed_spin_rad_s": (-spin).tolist(),
    "frame_convention": "DOE observation coordinates: x image-right, y image-up, visible hemisphere z>0; not table coordinates",
    "forward_pose_residual_deg_on_ransac_inliers": {"upstream_median": float(np.median(raw_err)), "upstream_mean": float(np.mean(raw_err)), "sign_reversed_median": float(np.median(reverse_err)), "sign_reversed_mean": float(np.mean(reverse_err))},
    "known_motion_sign_check": synthetic_check, "random_seed_results": seed_results,
    "cpu_estimation_seconds": elapsed, "independent_ground_truth_available": False,
    "limitations": ["One short official cropped-frame sequence; not event-camera data", "Pose residuals are internal consistency, not true spin accuracy", "Reversed sign shown as a separate diagnostic; upstream sources unchanged", "No table/world extrinsic calibration"],
}
(OUT / "results.json").write_text(json.dumps(summary, indent=2)+"\n")
print(json.dumps(summary, indent=2), flush=True)
