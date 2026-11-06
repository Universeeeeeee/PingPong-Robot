"""Real-event sensitivity to unverified center assumptions; no physical truth."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
GEOM = OUT.parent/'geometry_observation_20261101'
DATA = ROOT/'g1_spin/data/decoded/local_analysis_222403'
sys.path.insert(0, str(ROOT/'g1_spin/src'))
sys.path.insert(0, str(OUT.parent/'translation_compensation_20261102'))
from g1_spin.spherical_cmax import estimate, lift, iwe, warp
from translation import centers_at

SEGMENTS = {'translation': (.30, .60), 'early_turn': (1.60, 1.90),
            'middle': (6.20, 6.50), 'late': (11.20, 11.50)}
VARIANTS = ['A_fixed', 'B_50ms', 'C_100ms', 'B_cx_plus2', 'B_cx_minus2',
            'B_cy_plus2', 'B_cy_minus2', 'B_velocity_075', 'B_velocity_125']
RADIUS = 48.


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_write(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def csv_write(name, rows):
    with (OUT/name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def file_hashes():
    paths = [ROOT/'g1_spin/src/g1_spin/spherical_cmax.py',
             OUT.parent/'translation_compensation_20261102/translation.py', DATA/'manifest.json']
    for folder in [GEOM, OUT.parent/'boundary_support_20261101',
                   OUT.parent/'translation_compensation_20261102']:
        paths += sorted(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    manifest = json.loads((DATA/'manifest.json').read_text())
    paths += [DATA/p for p in manifest['parts']]
    return {str(p.relative_to(ROOT)): digest(p) for p in paths}


def prepare():
    with (GEOM/'observations.csv').open() as f:
        observations = list(csv.DictReader(f))
    models, diagnostics = {}, []
    fig, axes = plt.subplots(4, 2, figsize=(12, 12))
    for row_idx, (segment, (start, end)) in enumerate(SEGMENTS.items()):
        anchor = (start+end)/2
        models[segment] = {}
        fitted_rows = {}
        for duration in [50, 100]:
            rows = [r for r in observations if r['segment'] == segment and
                    r['duration_ms'] == str(duration) and r['mode'] == 'fixed48' and
                    r['threshold_policy'] == 'fixed_0.6']
            rows.sort(key=lambda r: float(r['midpoint_s']))
            assert len(rows) == 31 and all(r['optimizer_valid'] == 'True' for r in rows)
            t = np.array([float(r['midpoint_s']) for r in rows])
            xy = np.array([[float(r['cx']), float(r['cy'])] for r in rows])
            design = np.column_stack((np.ones(len(t)), t-anchor))
            coeff = np.linalg.lstsq(design, xy, rcond=None)[0]
            fitted = design@coeff
            model = dict(anchor_s=anchor, center_anchor_px=coeff[0].tolist(),
                         velocity_px_s=coeff[1].tolist(), coverage_s=[float(t[0]), float(t[-1])],
                         duration_ms=duration, threshold_policy='fixed_0.6', radius_px=RADIUS,
                         hypothesis_only=True, physics_validated=False)
            models[segment][str(duration)] = model
            fitted_rows[duration] = xy
            residual = np.linalg.norm(xy-fitted, axis=1)
            diagnostics.append(dict(segment=segment, observation_ms=duration,
                cx_span_px=float(np.ptp(xy[:, 0])), cy_span_px=float(np.ptp(xy[:, 1])),
                vx_px_s=float(coeff[1, 0]), vy_px_s=float(coeff[1, 1]),
                predicted_50ms_motion_px=float(np.linalg.norm(coeff[1])*.05),
                line_residual_rms_px=float(np.sqrt(np.mean(residual**2))),
                line_residual_max_px=float(residual.max()), hypothesis_only=True, physics_validated=False))
            for dim in range(2):
                axes[row_idx, dim].plot(t, xy[:, dim], '.-', alpha=.55, label=f'{duration}ms candidate')
                axes[row_idx, dim].plot(t, fitted[:, dim], '--', label=f'{duration}ms linear hypothesis')
                axes[row_idx, dim].set(title=f'{segment}: c{"xy"[dim]} (unverified)', xlabel='Relative event time (s)', ylabel='px')
                axes[row_idx, dim].legend(fontsize=7)
        distance = np.linalg.norm(fitted_rows[50]-fitted_rows[100], axis=1)
        for d in diagnostics[-2:]:
            d.update(paired_50_100_median_px=float(np.median(distance)), paired_50_100_max_px=float(distance.max()))
    fig.suptitle('Unverified center hypotheses; nominal midpoints are not instantaneous measurements')
    fig.tight_layout(); fig.savefig(OUT/'center_hypotheses.png', dpi=150); plt.close(fig)
    json_write('models.json', models); csv_write('geometry_diagnostics.csv', diagnostics)
    environment = dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
        numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__,
        protected_input_sha256=file_hashes(), protocol_sha256=digest(OUT/'protocol.md'),
        run_source_sha256=digest(Path(__file__)), hypothesis_only=True, physics_validated=False)
    json_write('environment.json', environment)
    return models


def curve(models, segment, variant, query_us, reference_us):
    model = models[segment]['100' if variant == 'C_100ms' else '50']
    anchor = model['anchor_s']
    center = np.array(model['center_anchor_px'])
    velocity = np.array(model['velocity_px_s'])
    if variant == 'A_fixed':
        center = center+velocity*(reference_us/1e6-anchor); velocity = np.zeros(2)
    elif variant.startswith('B_cx'):
        center[0] += 2 if 'plus' in variant else -2
    elif variant.startswith('B_cy'):
        center[1] += 2 if 'plus' in variant else -2
    elif variant.startswith('B_velocity_'):
        center = center+velocity*(reference_us/1e6-anchor)
        anchor = reference_us/1e6
        velocity *= .75 if variant == 'B_velocity_075' else 1.25
    sample_us = np.rint(np.array(model['coverage_s'])*1e6).astype(np.int64)
    sample_xy = center+(sample_us/1e6-anchor)[:, None]*velocity
    result = centers_at(sample_us, sample_xy, np.asarray(query_us, dtype=np.int64), max_gap_us=300000)
    independent = center+(np.asarray(query_us)/1e6-anchor)[:, None]*velocity
    assert np.max(np.abs(result-independent)) < 1e-10
    return result


def candidate_omega(result):
    omega = result['omega_rad_s'] if result['valid'] else result['quality'].get('rejected_candidate_omega_rad_s')
    return None if omega is None else np.asarray(omega, float)


def angle(a, b):
    if a is None or b is None or min(np.linalg.norm(a), np.linalg.norm(b)) < 1e-12:
        return None
    return float(np.rad2deg(np.arccos(np.clip(np.dot(a, b)/np.linalg.norm(a)/np.linalg.norm(b), -1, 1))))


def optimize(models):
    environment = json.loads((OUT/'environment.json').read_text())
    assert environment['protocol_sha256'] == digest(OUT/'protocol.md')
    assert environment['run_source_sha256'] == digest(Path(__file__))
    assert environment['protected_input_sha256'] == file_hashes()
    manifest = json.loads((DATA/'manifest.json').read_text())
    events = np.concatenate([np.load(DATA/p) for p in manifest['parts']])
    origin_us = int(events['t'][0]); rel_us = events['t'].astype(np.int64)-origin_us
    results, audits, detailed, checks = [], [], [], []
    plots = []
    calls = 0
    for segment, (start, end) in SEGMENTS.items():
        if segment == 'translation':
            continue
        for offset in [.05, .20]:
            lo, hi = round((start+offset)*1e6), round((start+offset+.05)*1e6)
            reference = (lo+hi)//2
            a, b = np.searchsorted(rel_us, [lo, hi])
            raw, raw_t = events[a:b], rel_us[a:b]
            raw_xy = np.column_stack((raw['x'], raw['y'])).astype(float)
            roi = ((raw_xy[:, 0]>=475)&(raw_xy[:, 0]<638)&(raw_xy[:, 1]>=143)&(raw_xy[:, 1]<285))
            window = f'{segment}_{lo/1e6:.2f}_{hi/1e6:.2f}'
            centers = {v: curve(models, segment, v, raw_t, reference) for v in VARIANTS}
            distance = {v: np.linalg.norm(raw_xy-c, axis=1) for v, c in centers.items()}
            common = roi & (raw_xy[:, 1]<240)
            for v in VARIANTS:
                common &= distance[v] <= .8*RADIUS
            full_ids = np.arange(a, b, dtype=np.int64)
            common_ids = full_ids[common]
            selected = np.flatnonzero(common)
            if len(selected)>5000:
                selected = selected[np.linspace(0, len(selected)-1, 5000).astype(int)]
            ids, xy, times = full_ids[selected], raw_xy[selected], raw_t[selected]
            assert len(ids) >= 100
            np.savez_compressed(OUT/f'{window}_events.npz', raw_global_event_ids=full_ids,
                common_global_event_ids=common_ids, optimizer_global_event_ids=ids,
                x=xy[:, 0], y=xy[:, 1], relative_time_us=times, polarity=raw['p'][selected])
            event_hash = hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest()
            for v in VARIANTS:
                audits.append(dict(window=window, variant=v, raw_window_events=len(raw), roi_events=int(roi.sum()),
                    raw_y_lt240_events=int((roi&(raw_xy[:, 1]<240)).sum()),
                    outside_sphere_whole_window=int((distance[v]>RADIUS).sum()),
                    outside_sphere_roi=int((roi&(distance[v]>RADIUS)).sum()),
                    outside_inner_roi=int((roi&(distance[v]>.8*RADIUS)).sum()),
                    common_events=len(common_ids), common_fraction_whole_window=len(common_ids)/len(raw),
                    common_fraction_roi=len(common_ids)/int(roi.sum()), optimizer_events=len(ids),
                    optimizer_id_sha256=event_hash, hypothesis_only=True, physics_validated=False))
            fitted = {}
            for v in VARIANTS:
                centers_now = centers[v][selected]
                ref_center = curve(models, segment, v, [reference], reference)[0]
                aligned = xy-centers_now+ref_center
                points = lift(aligned[:, 0], aligned[:, 1], *ref_center, RADIUS)
                assert np.max(np.abs(points[:, :2]-(xy-centers_now)/RADIUS)) < 1e-13
                dt = (times-reference)*1e-6
                tick = perf_counter(); result = estimate(points, dt, max_rpm=600, initial=None)
                elapsed = perf_counter()-tick; calls += 1
                fitted[v] = result
                omega = candidate_omega(result)
                q = result['quality']
                row = dict(window=window, segment=segment, start_s=lo/1e6, end_s=hi/1e6,
                    reference_s=reference/1e6, variant=v, valid=result['valid'], failure_reason=result['failure_reason'],
                    candidate_kind='accepted' if result['valid'] else 'rejected_optimizer_candidate',
                    wx_rad_s=None if omega is None else float(omega[0]),
                    wy_rad_s=None if omega is None else float(omega[1]),
                    wz_rad_s=None if omega is None else float(omega[2]),
                    candidate_rpm=None if omega is None else float(np.linalg.norm(omega)*60/(2*np.pi)),
                    contrast_before=q.get('contrast_before'), contrast_after=q.get('contrast_after'),
                    contrast_relative_gain=q.get('relative_improvement'), runtime_s=elapsed,
                    event_count=len(ids), event_id_sha256=event_hash,
                    hypothesis_only=True, physics_validated=False)
                for base in ['A_fixed', 'B_50ms']:
                    base_result = fitted[base] if base in fitted else None
                    base_omega = None if base_result is None else candidate_omega(base_result)
                    row[f'directed_angle_to_{base}_deg'] = angle(omega, base_omega)
                    row[f'omega_delta_to_{base}_rad_s'] = None if omega is None or base_omega is None else float(np.linalg.norm(omega-base_omega))
                    row[f'rpm_delta_to_{base}'] = None if omega is None or base_omega is None else float((np.linalg.norm(omega)-np.linalg.norm(base_omega))*60/(2*np.pi))
                    row[f'contrast_after_ratio_to_{base}'] = None if base_result is None else q.get('contrast_after')/base_result['quality']['contrast_after']
                results.append(row)
                detailed.append(dict(window=window, variant=v, hypothesis_only=True, physics_validated=False, result=result))
                print(f'{calls}/54 {window} {v} valid={result["valid"]} reason={result["failure_reason"]}', flush=True)
                if v == 'B_50ms':
                    plots.append((window, xy.copy(), aligned.copy(), dt.copy()))
            # Fill A-to-B comparison after B is available, using the same candidates.
            arow = results[-len(VARIANTS)]
            aw, bw = candidate_omega(fitted['A_fixed']), candidate_omega(fitted['B_50ms'])
            arow['directed_angle_to_B_50ms_deg'] = angle(aw, bw)
            arow['omega_delta_to_B_50ms_rad_s'] = None if aw is None or bw is None else float(np.linalg.norm(aw-bw))
            arow['rpm_delta_to_B_50ms'] = None if aw is None or bw is None else float((np.linalg.norm(aw)-np.linalg.norm(bw))*60/(2*np.pi))
            arow['contrast_after_ratio_to_B_50ms'] = fitted['A_fixed']['quality']['contrast_after']/fitted['B_50ms']['quality']['contrast_after']
            if not checks:
                velocity = np.array([120., -70.])
                sample_t = np.array([lo, hi], dtype=np.int64)
                known_translation = (sample_t-reference)[:, None]*1e-6*velocity
                imposed = (times-reference)[:, None]*1e-6*velocity
                moved = xy+imposed
                queried = centers_at(sample_t, known_translation, times, max_gap_us=50000)
                recovered = moved-queried
                error = np.linalg.norm(recovered-xy, axis=1)
                assert error.max() < 1e-12
                checks.append(dict(window=window, kind='coordinate_equivariance_only_not_physical_DVS',
                    event_count=len(ids), event_id_sha256=event_hash, imposed_velocity_px_s=velocity.tolist(),
                    max_restoration_px=float(error.max()), rms_restoration_px=float(np.sqrt(np.mean(error**2))),
                    timestamps_and_ids_unchanged=True, physics_validated=False))
                fig, ax = plt.subplots(1, 3, figsize=(12, 4))
                for axis, coords, title in zip(ax, [xy, moved, recovered], ['Original events', 'Known coordinate translation', 'Known term removed']):
                    axis.scatter(*coords.T, c=dt*1e3, s=1, cmap='viridis'); axis.set(title=title, aspect='equal'); axis.invert_yaxis()
                fig.suptitle('Same real event IDs; coordinate self-check only'); fig.tight_layout()
                fig.savefig(OUT/'coordinate_selfcheck.png', dpi=150); plt.close(fig)
            csv_write('results.csv', results); csv_write('mask_audit.csv', audits)
            json_write('detailed_results.json', detailed)
    assert calls == 54
    assert environment['protected_input_sha256'] == file_hashes(), 'Protected input changed'
    json_write('validation.json', dict(calls=calls, coordinate_selfcheck=checks,
        protected_inputs_unchanged=True, same_event_ids_all_variants=True,
        independent_linear_interpolation_max_tolerance_px=1e-10,
        aligned_vs_direct_lift_xy_tolerance=1e-13, hypothesis_only=True, physics_validated=False))
    fig, axes = plt.subplots(6, 2, figsize=(10, 20))
    for axes_row, (window, original, aligned, dt) in zip(axes, plots):
        for ax, coords, title in zip(axes_row, [original, aligned], ['Raw common events', 'B nominal center alignment']):
            ax.scatter(*coords.T, c=dt*1e3, s=.6, cmap='viridis'); ax.set(title=f'{window}: {title}', aspect='equal'); ax.invert_yaxis()
    fig.suptitle('Unverified B hypothesis; same selected events per row, no physical accuracy claim'); fig.tight_layout()
    fig.savefig(OUT/'real_event_alignment.png', dpi=140); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for window in dict.fromkeys(r['window'] for r in results):
        rr = [r for r in results if r['window']==window]
        axes[0].plot(range(9), [r['directed_angle_to_A_fixed_deg'] for r in rr], '.-', label=window)
        axes[1].plot(range(9), [r['candidate_rpm'] for r in rr], '.-', label=window)
    for ax in axes:
        ax.set_xticks(range(9), VARIANTS, rotation=60, ha='right'); ax.legend(fontsize=6)
    axes[0].set_ylabel('Directed candidate axis difference to A (deg)')
    axes[1].set_ylabel('Optimizer candidate RPM (includes rejected)')
    fig.suptitle('Method sensitivity only; rejected candidates are not valid spin estimates'); fig.tight_layout()
    fig.savefig(OUT/'axis_sensitivity.png', dpi=150); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--prepare', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        prepare()
    else:
        optimize(json.loads((OUT/'models.json').read_text()))
