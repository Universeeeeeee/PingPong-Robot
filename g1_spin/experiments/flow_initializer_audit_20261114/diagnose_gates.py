"""Replay local time-plane gates on saved original-pixel event windows."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
FRONT = ROOT / 'g1_spin/experiments/drill_phase4_20261106/real_frontend'
sys.path.insert(0, str(ROOT / 'g1_spin/src'))
from g1_spin.optical_flow import solve_constraints


def distribution(values):
    return dict(count=len(values), median=float(np.median(values)) if values else None,
                p95=float(np.percentile(values, 95)) if values else None)


def extract(data, half_span, ratio, representation):
    xy = data['xy']
    times = data['times']
    polarity = data['polarity']
    radius = float(data['radius'])
    tree = cKDTree(xy)
    points, gradients, accepted_spatial = [], [], []
    failures = Counter()
    residual_ratios, lower_ratios, neighbors, spatial_speeds = [], [], [], []
    seeds = np.linspace(0, len(times) - 1, min(300, len(times)), dtype=int)
    for seed in seeds:
        indices = np.asarray(tree.query_ball_point(xy[seed], 4), dtype=int)
        indices = indices[(np.abs(times[indices] - times[seed]) <= half_span) &
                          (polarity[indices] == polarity[seed])]
        if representation == 'latest' and len(indices):
            order = np.lexsort((times[indices], xy[indices, 1], xy[indices, 0]))
            sorted_indices = indices[order]
            pixels = xy[sorted_indices]
            last = np.r_[np.any(pixels[1:] != pixels[:-1], axis=1), True]
            indices = sorted_indices[last]
        neighbors.append(len(indices))
        if len(indices) < 10:
            failures['too_few_neighbors'] += 1
            continue
        coords = xy[indices]
        center = coords.mean(axis=0)
        centered = coords - center
        singular = np.linalg.svd(centered, compute_uv=False)
        if singular[0] < 1e-6 or singular[1] / singular[0] < .08:
            failures['spatial_degeneracy'] += 1
            continue
        dt = times[indices] - times[seed]
        std = float(np.std(dt))
        if std < 1e-6:
            failures['time_degeneracy'] += 1
            continue
        design = np.column_stack([centered, np.ones(len(indices))])
        coeff, _, rank, _ = np.linalg.lstsq(design, dt, rcond=None)
        rms = float(np.sqrt(np.mean((design @ coeff - dt)**2)))
        gate = max(20e-6, ratio * std)
        residual_ratios.append(rms / std)
        _, inverse, counts = np.unique(coords, axis=0, return_inverse=True, return_counts=True)
        means = np.bincount(inverse, weights=dt) / counts
        lower = float(np.sqrt(np.mean((dt - means[inverse])**2)))
        lower_ratios.append(lower / gate)
        if rank < 3 or rms > gate:
            failures['plane_residual'] += 1
            continue
        sx, sy = (center - data['center']) / radius
        if sx*sx + sy*sy >= .9**2:
            failures['sphere_boundary'] += 1
            continue
        norm = np.linalg.norm(coeff[:2])
        if norm < 1e-8:
            failures['vanishing_gradient'] += 1
            continue
        points.append([sx, sy, -np.sqrt(1-sx*sx-sy*sy)])
        gradients.append(coeff[:2])
        accepted_spatial.append(center.tolist())
        spatial_speeds.append(1 / norm)
    solve = solve_constraints(np.array(points).reshape(-1,3),
                              np.array(gradients).reshape(-1,2), radius)
    candidate = None
    if len(points) >= 3:
        matrix = radius * np.cross(points, np.c_[gradients, np.zeros(len(gradients))])
        omega, _, rank, singular = np.linalg.lstsq(matrix, np.ones(len(points)), rcond=None)
        speed = np.linalg.norm(omega)
        candidate = dict(omega_rad_s=omega.tolist(), rpm=float(speed*60/(2*np.pi)),
                         axis=(omega/speed).tolist() if speed > 1e-12 else None,
                         rank=int(rank), condition=float(singular[0]/singular[-1]) if singular[-1]>1e-12 else None,
                         constraint_rms=float(np.sqrt(np.mean((matrix@omega-1)**2))))
    return dict(seeds=len(seeds), accepted=len(points), failures=dict(failures),
                local_rms_over_std=distribution(residual_ratios),
                lower_bound_exceeds_gate=sum(x>1 for x in lower_ratios),
                eligible_for_lower_bound=len(lower_ratios),
                median_neighbors=float(np.median(neighbors)) if neighbors else 0,
                normal_speed_px_s=distribution(spatial_speeds),
                accepted_spatial_centers=accepted_spatial, production_solver=solve,
                unchecked_candidate=candidate)


def main():
    manifest = json.loads((FRONT/'input_manifest.json').read_text())
    configs = [(span, .15) for span in [.0005, .001, .002]] + [(.002, ratio) for ratio in [.3, .5, .8]]
    rows, hashes = [], {}
    for record in manifest:
        path = FRONT/record['input_file']
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        data = np.load(path)
        for representation in ['all', 'latest']:
            for span, ratio in configs:
                result = extract(data, span, ratio, representation)
                rows.append(dict(case=record['case'], condition=record['condition'],
                                 source=record['source'], source_mask=record['source_mask'],
                                 events=len(data['times']), half_span_ms=span*1000,
                                 residual_ratio_limit=ratio, representation=representation, **result))
        print(record['case'], len(data['times']), flush=True)
    output = dict(rows=rows, sources=hashes,
                  scope='Offline sensitivity diagnostic on old fixed center geometry; changing gates is not validation of the resulting motion',
                  latest_definition='Latest timestamp per original pixel and polarity inside the same symmetric neighborhood; offline noncausal representation')
    (OUT/'gate_diagnostics.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h for p,h in hashes.items())
    print('DONE',len(rows),flush=True)

if __name__ == '__main__':
    main()
