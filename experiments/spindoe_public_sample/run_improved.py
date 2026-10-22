"""Compare multiple-basis matching against the saved official-sample baseline."""
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs/improved"
OUT.mkdir(exist_ok=True)
os.environ["MPLBACKEND"] = "Agg"
os.environ["MPLCONFIGDIR"] = str(ROOT / ".cache/matplotlib")
sys.path.insert(0, str(ROOT / "upstream/python"))

import cv2
import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial.transform import Rotation
from baygeohasher import BayGeoHasher
from corrected_spin import BoundedSpinRegressor
from point_matcher import match_points

# Reuse the SAME network heatmaps to isolate the geometry change.
data = np.load(ROOT / "outputs/baseline_estimates.npz")
timestamps = data["timestamps_ns"]
times = (timestamps - timestamps[0]) * 1e-9
heatmaps = data["heatmaps"]
reference = np.loadtxt(ROOT / "upstream/cad/ref_pattern.csv", delimiter=",")
hasher = BayGeoHasher()
hasher.gen_hash_table(reference)

# Exact blob extraction parameters from upstream DOE.init_blob_detector().
params = cv2.SimpleBlobDetector_Params()
params.filterByCircularity = params.filterByConvexity = params.filterByInertia = False
params.filterByArea = params.filterByColor = True
params.minArea, params.maxArea = 1, 30
params.blobColor, params.minDistBetweenBlobs = 255, 1
detector = cv2.SimpleBlobDetector_create(params)

poses, frame_results, points_by_frame = [], [], []
start = time.perf_counter()
for i, heatmap in enumerate(heatmaps):
    keypoints = detector.detect(cv2.inRange(heatmap, 0.1, 255))
    xy = np.array([[2*(k.pt[0]/60-0.5), -2*(k.pt[1]/60-0.5)] for k in keypoints])
    points = np.c_[xy, np.sqrt(np.maximum(0, 1-np.sum(xy*xy, axis=1)))] if len(xy) else np.empty((0,3))
    result = match_points(points, hasher)
    accepted = result is not None and result["rmse"] <= 0.12
    poses.append(result["rotation"] if accepted else None)
    points_by_frame.append(points)
    frame_results.append({"frame": i, "detected_points": len(points), "accepted": accepted,
                          "rmse": result["rmse"] if result else None,
                          "competing_rmse": result["competing_rmse"] if result else None,
                          "bases_tested": result["bases_tested"] if result else 0})
pose_seconds = time.perf_counter()-start
valid = np.array([i for i, pose in enumerate(poses) if pose is not None], dtype=int)
np.random.seed(0)
# Explicit hypothesis range for this recorded sample; it does not prove that
# higher, aliased speeds are absent. Use original capture intervals, not inliers.
regressor = BoundedSpinRegressor(max_spin_rps=0.5/np.median(np.diff(times)))
spin, inlier_positions = regressor.RANSAC_regress(times[valid], [poses[i] for i in valid])
if spin is None:
    inliers = np.array([], dtype=int)
    predictions = [None]*len(times)
else:
    inliers = valid[inlier_positions]
    anchor = inliers[0]
    predictions = Rotation.from_rotvec((times-times[anchor])[:,None]*spin) * poses[anchor]

def pose_residual(indices):
    if spin is None or len(indices) == 0:
        return None
    return float(np.median([np.degrees((predictions[i]*poses[i].inv()).magnitude()) for i in indices]))

# Verify order robustness on the actual observations, not only synthetic points.
rng = np.random.default_rng(2026)
order_check = []
for i in (0, len(times)//2, len(times)-1):
    result = match_points(rng.permutation(points_by_frame[i]), hasher)
    original = frame_results[i]
    if result is None:
        consistent = original["rmse"] is None
    else:
        consistent = original["rmse"] is not None and abs(result["rmse"]-original["rmse"]) < 1e-12
        if poses[i] is not None:
            consistent = consistent and (result["rotation"]*poses[i].inv()).magnitude() < 1e-9
    order_check.append({"frame": i, "unchanged": bool(consistent)})

baseline = json.loads((ROOT / "outputs/results.json").read_text())
summary = {"frame_count": len(times), "same_network_heatmaps": True,
           "baseline_doe_pose_count": baseline["unmodified_algorithm_baseline"]["doe_pose_count"],
           "baseline_spin_available": baseline["unmodified_algorithm_baseline"]["spin_available"],
           "size_sorted_diagnostic_pose_count": baseline["doe_pose_count"],
           "improved_pose_count": len(valid), "improved_spin_inlier_count": len(inliers),
           "estimated_rps": None if spin is None else float(np.linalg.norm(spin)/(2*np.pi)),
           "spin_observation_frame_rad_s": None if spin is None else spin.tolist(),
           "forward_pose_residual_median_deg_inliers": pose_residual(inliers),
           "forward_pose_residual_median_deg_all_accepted": pose_residual(valid),
           "order_checks": order_check, "pose_matching_seconds": pose_seconds,
           "phase_fit": regressor.diagnostics,
           "independent_ground_truth_available": False, "logo_detector_implemented": False,
           "frames": frame_results}
(OUT / "results.json").write_text(json.dumps(summary, indent=2)+"\n")
np.savez_compressed(OUT / "estimates.npz", timestamps_ns=timestamps,
                    quaternions_xyzw=np.array([pose.as_quat() if pose else [np.nan]*4 for pose in poses]),
                    inlier_indices=inliers, spin_rad_s=np.array([]) if spin is None else spin)

def overlay(image, rotation):
    image = image.copy()
    if rotation is None:
        return image
    for point in rotation.apply(reference):
        if point[2] > 0:
            x = int((point[0]+1)*0.5*image.shape[1])
            y = int((1-point[1])*0.5*image.shape[0])
            cv2.circle(image, (x,y), 1, (50,230,240), -1)
    return image

font = ImageFont.load_default(size=17)
small = ImageFont.load_default(size=14)
frames=[]
for i, timestamp in enumerate(timestamps):
    original = np.array(Image.open(ROOT / f"upstream/data/test/{timestamp}.png").convert("RGB"))
    canvas=Image.new("RGB",(960,420),(18,22,29))
    draw=ImageDraw.Draw(canvas)
    candidate = "no spin estimate" if spin is None else f"{summary['estimated_rps']:.2f} rps candidate"
    draw.text((20,12),f"Multiple-basis matching | {candidate} | no independent ground truth",font=font,fill="white")
    for column,(label, img) in enumerate([("Input",original),("Per-frame point matching",overlay(original,poses[i])),("Temporal spin prediction",overlay(original,predictions[i]))]):
        x=column*320+22
        draw.text((x,52),label,font=small,fill="white")
        canvas.paste(Image.fromarray(img).resize((272,272),Image.Resampling.NEAREST),(x,80))
    draw.text((22,370),f"Frame {i+1}/22 | {times[i]*1000:.3f} ms | pose valid: {poses[i] is not None} | spin inlier: {i in inliers}",font=font,fill="white")
    draw.text((22,398),"60x slow motion; same detector; multiple bases + bounded circular-phase fitting.",font=small,fill=(180,190,205))
    frames.append(canvas)
playback=np.arange(0,times[-1]+np.median(np.diff(times)),1/(30*60))
with imageio.get_writer(str(OUT/"improved_60x.mp4"),fps=30,codec="libx264",quality=8,macro_block_size=4) as writer:
    for index in np.searchsorted(times,playback,side="right")-1:
        writer.append_data(np.array(frames[index]))
frames[8].save(OUT/"comparison_frame.png")
print(json.dumps({k:v for k,v in summary.items() if k!="frames"},indent=2))
