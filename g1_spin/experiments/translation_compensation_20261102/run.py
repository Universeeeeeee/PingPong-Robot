"""Controlled geometry experiment, not a physical DVS or a real G1 benchmark."""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"g1_spin/src"))
from g1_spin.synthetic import generate
from g1_spin.spherical_cmax import estimate, lift, warp
from translation import compensate, centers_at, load_centers

OUT = Path(__file__).resolve().parent
META = dict(image_plane="G1 image plane", center_semantics="ball_center",
            timebase="synthetic integer microseconds from first generated sample",
            source="synthetic.generate ideal point geometry plus imposed translation",
            synchronization_verified=True,
            synchronization_basis="constructed using exactly the same integer event times",
            center_error="exact imposed synthetic center; real center error unknown",
            data_kind="controlled_synthetic_not_G1_recording")


def save_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fit(points, times, truth):
    start = perf_counter()
    result = estimate(points, times, max_rpm=1500)
    runtime = perf_counter()-start
    omega = result["omega_rad_s"] if result["valid"] else result["quality"].get("rejected_candidate_omega_rad_s")
    omega = np.asarray(omega, float) if omega is not None else np.full(3, np.nan)
    speed = np.linalg.norm(omega)
    axis_error = np.rad2deg(np.arccos(np.clip(np.dot(omega, truth)/(speed*np.linalg.norm(truth)), -1, 1)))
    return dict(valid=result["valid"], failure_reason=result["failure_reason"],
                rpm=float(speed*60/(2*np.pi)), axis_error_deg=float(axis_error),
                omega_error_rad_s=float(np.linalg.norm(omega-truth)), omega_rad_s=omega.tolist(),
                relative_improvement=result["quality"].get("relative_improvement"),
                runtime_s=runtime)


def main():
    radius, center, rpm = 100., np.array([128., 128.]), 1000.
    axis = np.array([.3, .7, .6]); axis /= np.linalg.norm(axis)
    truth = axis*rpm*2*np.pi/60
    batch = generate(rpm=rpm, axis=axis, duration_s=.04, sample_hz=3000, seed=0,
                     center=center, radius=radius)
    xy, times = np.column_stack((batch.x, batch.y)), batch.t_us
    reference = int((int(times[0])+int(times[-1]))//2)
    dt = np.array([int(t)-reference for t in times])*1e-6
    direction = np.array([1., .4]); direction /= np.linalg.norm(direction)
    speeds = [0., 100., 500.]
    moved = {v: xy+dt[:, None]*v*direction for v in speeds}
    raw_valid = {v: np.sum(((coords-center)/radius)**2, axis=1) < 1-1e-12
                 for v, coords in moved.items()}
    common = np.logical_and.reduce(list(raw_valid.values()))
    event_ids = np.flatnonzero(common)
    common_xy, common_times, common_dt = xy[common], times[common], dt[common]
    baseline_points = lift(common_xy[:, 0], common_xy[:, 1], *center, radius)
    base = fit(baseline_points, common_dt, truth)
    # Sparse linear samples at 0/10/20/30/40 ms; the last event is earlier than 40 ms.
    sample_times = np.arange(0, 40001, 10000, dtype=np.int64)
    summaries, csv_rows, restoration = [], [], {}
    for speed in speeds:
        sample_centers = center+(sample_times-reference)[:, None]*1e-6*speed*direction
        csv_path = OUT/f"synthetic_centers_{int(speed)}px_s.csv"
        save_csv(csv_path, [dict(time_us=int(t), cx_px=float(c[0]), cy_px=float(c[1]))
                            for t, c in zip(sample_times, sample_centers)])
        (OUT/"synthetic_centers_metadata.json").write_text(json.dumps(META, indent=2)+"\n")
        loaded_t, loaded_xy, meta = load_centers(csv_path, OUT/"synthetic_centers_metadata.json")
        start = perf_counter()
        corrected, reference_center = compensate(moved[speed][common], common_times,
                        loaded_t, loaded_xy, reference_us=reference, metadata=meta,
                        event_timebase=META["timebase"], max_gap_us=10000)
        correction_s = perf_counter()-start
        point_error = np.linalg.norm(corrected-common_xy, axis=1)
        restoration[speed] = corrected
        variants = {"uncompensated": moved[speed][common], "compensated": corrected}
        methods = {}
        for name, coords in variants.items():
            points = lift(coords[:, 0], coords[:, 1], *reference_center, radius)
            fitted = fit(points, common_dt, truth)
            methods[name] = fitted
            csv_rows.append(dict(speed_px_s=speed, method=name, event_count=len(event_ids),
                         **{k: v for k, v in fitted.items() if k != "omega_rad_s"},
                         max_restore_error_px=float(point_error.max()) if name == "compensated" else None))
        summaries.append(dict(speed_px_s=speed, full_input_outside_fixed_sphere=int((~raw_valid[speed]).sum()),
                              correction_runtime_s=correction_s, max_restore_error_px=float(point_error.max()),
                              rms_restore_error_px=float(np.sqrt(np.mean(point_error**2))), methods=methods))
    # A 10% velocity underestimate leaves exactly 10% of imposed displacement.
    speed = 500.
    wrong_centers = center+(sample_times-reference)[:, None]*1e-6*speed*.9*direction
    wrong, _ = compensate(moved[speed][common], common_times, sample_times, wrong_centers,
                       reference_us=reference, metadata=META, event_timebase=META["timebase"], max_gap_us=10000)
    wrong_error = np.linalg.norm(wrong-common_xy, axis=1)
    # Constant clock offset cancels for linear translation; acceleration makes it observable.
    fine_t = np.arange(0, 42001, 100, dtype=np.int64)
    acceleration = np.array([2000., -1000.])
    def quadratic(t):
        rel = (t-reference)*1e-6
        return center+rel[:, None]*speed*direction+.5*rel[:, None]**2*acceleration
    true_curve = quadratic(common_times)
    accelerated_events = common_xy+true_curve-center
    shifted_query = np.concatenate([common_times+1000, np.array([reference+1000])])
    shifted_centers = centers_at(fine_t, quadratic(fine_t), shifted_query, max_gap_us=100)
    shifted_corrected = accelerated_events-(shifted_centers[:-1]-shifted_centers[-1])
    sync_error = np.linalg.norm(shifted_corrected-common_xy, axis=1)
    expected_sync = np.abs(common_dt)*.001*np.linalg.norm(acceleration)
    # Actual lift also subtracts the (wrong) absolute reference center.
    # Its pixel-space input error is c(t)-c(t+delta), including v*delta.
    prelift_error = (shifted_corrected-shifted_centers[-1])-(common_xy-center)
    expected_prelift = quadratic(common_times)-quadratic(common_times+1000)
    linear_samples = center+(fine_t-reference)[:, None]*1e-6*speed*direction
    linear_shift = centers_at(fine_t, linear_samples, shifted_query, max_gap_us=100)
    linear_aligned = moved[speed][common]-(linear_shift[:-1]-linear_shift[-1])
    linear_prelift_error = (linear_aligned-linear_shift[-1])-(common_xy-center)
    # Independently verify SO(3) and IWE/estimator did not remove rotation.
    # Pixel restoration implies identical geometry; estimate comparison is separately measured.
    corrected_points = lift(restoration[500.][:, 0], restoration[500.][:, 1], *center, radius)
    warp_difference = np.max(np.linalg.norm(warp(corrected_points, common_dt, truth)-
                                            warp(baseline_points, common_dt, truth), axis=1))
    validation = dict(input_event_count=len(times), common_event_count=len(event_ids),
                common_excluded_count=int((~common).sum()), reference_us=reference,
                common_event_id_sha256=hashlib.sha256(event_ids.astype("<i8").tobytes()).hexdigest(),
                truth_rpm=rpm, truth_axis_camera=axis.tolist(), truth_omega_rad_s=truth.tolist(),
                baseline=base, cases=summaries, sensitivity=dict(
                    velocity_underestimate_fraction=.1, velocity_error_max_px=float(wrong_error.max()),
                    velocity_error_rms_px=float(np.sqrt(np.mean(wrong_error**2))),
                    clock_offset_us=1000, acceleration_px_s2=acceleration.tolist(),
                    accelerated_clock_error_max_px=float(sync_error.max()),
                    accelerated_clock_error_rms_px=float(np.sqrt(np.mean(sync_error**2))),
                    clock_closed_form_max_difference_px=float(np.max(np.abs(sync_error-expected_sync))),
                    accelerated_clock_prelift_error_max_px=float(np.linalg.norm(prelift_error, axis=1).max()),
                    accelerated_clock_prelift_closed_form_max_difference_px=float(np.max(np.abs(prelift_error-expected_prelift))),
                    linear_clock_relative_alignment_max_error_px=float(np.linalg.norm(linear_aligned-common_xy, axis=1).max()),
                    linear_clock_prelift_error_max_px=float(np.linalg.norm(linear_prelift_error, axis=1).max())),
                ground_truth_warp_max_baseline_difference= float(warp_difference))
    assert max(c["max_restore_error_px"] for c in summaries) < 1e-12
    assert warp_difference < 1e-12
    assert base["valid"], base
    for case in summaries:
        comp = case["methods"]["compensated"]
        assert comp["valid"], comp
        assert np.linalg.norm(np.array(comp["omega_rad_s"])-base["omega_rad_s"]) < .05
    assert np.max(np.abs(sync_error-expected_sync)) < 3e-6
    assert np.max(np.abs(prelift_error-expected_prelift)) < 3e-6
    assert np.max(np.abs(np.linalg.norm(linear_prelift_error, axis=1)-.5)) < 1e-12
    (OUT/"results.json").write_text(json.dumps(validation, indent=2, allow_nan=False)+"\n")
    save_csv(OUT/"results.csv", csv_rows)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)
    for ax, coords, title in zip(axes, [common_xy, moved[500.][common], restoration[500.]],
                                ["Fixed-center baseline", "+500 px/s translation", "Known-center correction"]):
        scatter = ax.scatter(*coords.T, c=common_dt*1000, s=2, cmap="viridis", rasterized=True)
        ax.add_patch(plt.Circle(center, radius, fill=False, color="gray", linewidth=.8))
        ax.set_title(title); ax.set_aspect("equal"); ax.set_xlabel("Synthetic image x (px)")
        ax.invert_yaxis()
    axes[0].set_ylabel("Synthetic image y (px)")
    fig.colorbar(scatter, ax=axes, label="Time from reference (ms)", shrink=.7)
    fig.suptitle(f"Ideal rotating points, not physical DVS | same {len(event_ids)} event IDs in every panel")
    fig.savefig(OUT/"compensation_comparison.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    sources = [Path(__file__), OUT/"translation.py", OUT/"test_translation.py",
               ROOT/"g1_spin/src/g1_spin/synthetic.py", ROOT/"g1_spin/src/g1_spin/spherical_cmax.py",
               ROOT/"g1_spin/src/g1_spin/events.py"]
    env = dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
               numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__,
               git_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
               source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
               synthetic_input_sha256=hashlib.sha256(b"".join(a.tobytes() for a in
                              [batch.x, batch.y, batch.t_us, batch.polarity])).hexdigest(),
               synthetic_parameters=dict(rpm=rpm, axis=axis.tolist(), duration_s=.04,
                                         sample_hz=3000, seed=0, radius_px=radius, center_px=center.tolist()))
    (OUT/"environment.json").write_text(json.dumps(env, indent=2)+"\n")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
