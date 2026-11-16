"""Fixed coarse spatial-normal initialization followed by a single CMAX search."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

OUT = Path(__file__).resolve().parent
EXPERIMENTS = OUT.parents[1]
PHASE3 = EXPERIMENTS / 'drill_phase3_20261106'
PHASE4 = EXPERIMENTS / 'drill_phase4_20261106/real_frontend'
MODELS = EXPERIMENTS / 'drill_phase5_20261106/local_motion/models.py'
sys.path.insert(0, str(PHASE3))
import methods

spec = importlib.util.spec_from_file_location('spatial_normal_models', MODELS)
models = importlib.util.module_from_spec(spec)
spec.loader.exec_module(models)
CASES = ['run1_01495000', 'run2_02095000', 'run3_01495000',
         'run4_02495000', 'run5_04095000', 'run5_01525000']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def normal_seed(xy, times, polarity, center, radius, patch_keys):
    patch_ids = np.floor(xy / 12).astype(int)
    observations = []
    for patch_number, key in enumerate(patch_keys):
        in_patch = np.all(patch_ids == key, axis=1)
        for segment in range(5):
            for pol in [0, 1]:
                selected = in_patch & (times >= segment * .002) & (times < (segment + 1) * .002) & (polarity == pol)
                pixels = xy[selected]
                record = dict(patch=patch_number, segment=segment, polarity=pol,
                              events=int(selected.sum()), eligible=False)
                if len(np.unique(pixels, axis=0)) < 12:
                    record['failure_reason'] = 'fewer_than_12_distinct_pixels'
                    observations.append(record)
                    continue
                fit = models.fit(pixels, times[selected], 'M_spatial')
                record['fit'] = fit
                if not fit['available']:
                    record['failure_reason'] = fit['failure_reason']
                    observations.append(record)
                    continue
                normalized = (np.asarray(fit['center']) - center) / radius
                if normalized @ normalized >= .9 ** 2:
                    record['failure_reason'] = 'outside_inner_sphere'
                    observations.append(record)
                    continue
                point = np.r_[normalized, -np.sqrt(1 - normalized @ normalized)]
                normal = np.r_[fit['normal'], 0.]
                record.update(eligible=True, constraint=(radius * np.cross(point, normal)).tolist(),
                              rhs_px_s=fit['speed_px_s'], failure_reason=None)
                observations.append(record)
    eligible = [r for r in observations if r['eligible']]
    result = dict(observations=observations, available=False, constraints=len(eligible),
                  patches=len(set(r['patch'] for r in eligible)))
    if not eligible:
        result['failure_reason'] = 'no_local_normal_observations'
        return result
    matrix = np.array([r['constraint'] for r in eligible])
    rhs = np.array([r['rhs_px_s'] for r in eligible])
    omega, _, rank, singular = np.linalg.lstsq(matrix, rhs, rcond=None)
    condition = float(singular[0] / singular[-1]) if rank == 3 else None
    speed = float(np.linalg.norm(omega))
    residual = matrix @ omega - rhs
    result.update(rank=int(rank), condition=condition, omega_rad_s=omega.tolist(),
                  rpm=speed * 60 / (2*np.pi), axis=(omega / speed).tolist() if speed > 1e-12 else None,
                  residual_rms_px_s=float(np.sqrt(np.mean(residual**2))),
                  residual_median_abs_px_s=float(np.median(abs(residual))),
                  rhs_rms_px_s=float(np.sqrt(np.mean(rhs**2))))
    if len(eligible) < 6 or result['patches'] < 3:
        result['failure_reason'] = 'insufficient_observation_or_patch_count'
    elif rank < 3 or condition > 100:
        result['failure_reason'] = 'rank_or_condition_failure'
    elif speed < 1e-6:
        result['failure_reason'] = 'zero_initial_direction'
    else:
        result.update(available=True, failure_reason=None)
    return result


def compact_search(search):
    return {k: v for k, v in search.items() if k not in ['objective_history']}


def main():
    if (OUT / 'results.json').exists():
        raise FileExistsError('Existing experiment results are preserved')
    manifest = {r['case']: r for r in json.loads((PHASE4 / 'input_manifest.json').read_text())}
    old = {r['window']: r for r in (json.loads(line) for line in (PHASE3 / 'results.jsonl').open()) if r['method'] == 'G'}
    sources = [Path(__file__), MODELS, PHASE3 / 'methods.py', PHASE3 / 'frozen/spherical_cmax.py']
    sources += [PHASE4 / manifest[case]['input_file'] for case in CASES]
    sources += [PHASE3 / 'windows' / f'{case}.npz' for case in CASES]
    protocol = dict(cases=CASES, representation='all events in complete old common-mask archive',
                    model='phase5 M_spatial pixel-normal regression', patch_px=12, top_patches=6,
                    patch_selection='descending event count over full10ms; deterministic lexicographic tie-break',
                    intervals_s=[[i*.002,(i+1)*.002] for i in range(5)], separate_polarity=True,
                    minimum_unique_pixels=12, initialization_gates='at least6 observations,3patches,rank3,condition<=100,nonzero',
                    residual_gate=None, sphere='old B1 fixed center, old radius; no motion correction',
                    search='phase3 axis_search initial directed axis magnitude scan plus exactly one3D Nelder-Mead, bound3500RPM',
                    search_events='exact archived phase3 selected events, at most5000; full common-mask events used for flow',
                    controls='run1 and translation: shuffle initializer times with seed20261114, refine clean events; additionally shuffle whole refinement',
                    status='exploratory coarse candidates, no held-out measurement acceptance',
                    source_hashes={str(p): digest(p) for p in sources})
    save('protocol.json', protocol)
    results = []
    for case in CASES:
        entry = manifest[case]
        data = np.load(PHASE4 / entry['input_file'])
        selected = np.load(PHASE3 / 'windows' / f'{case}.npz')
        xy, times, polarity = data['xy'], data['times'], data['polarity']
        center, radius = data['center'], float(data['radius'])
        patches, counts = np.unique(np.floor(xy / 12).astype(int), axis=0, return_counts=True)
        order = np.lexsort((patches[:,1], patches[:,0], -counts))[:6]
        patch_keys = patches[order]
        targets = selected['xy']
        points = methods.cmax.lift(targets[:,0], targets[:,1], *center, radius)
        clean_target_times = (selected['time_us'] - entry['start_us']) * 1e-6
        variants = ['clean']
        if case in ['run1_01495000', 'run5_04095000']:
            variants += ['shuffled_seed_clean_refinement', 'shuffled_seed_shuffled_refinement']
        shuffled_times = np.random.default_rng(20261114).permutation(times)
        time_by_id = dict(zip(data['ids'].tolist(), shuffled_times.tolist()))
        for variant in variants:
            seed_times = times if variant == 'clean' else shuffled_times
            started = perf_counter()
            initializer = normal_seed(xy, seed_times, polarity, center, radius, patch_keys)
            init_elapsed = perf_counter() - started
            if variant == 'shuffled_seed_shuffled_refinement':
                target_times = np.array([time_by_id[int(event_id)] for event_id in selected['ids']])
                sorted_ids = np.argsort(target_times, kind='stable')
                target_times, target_points = target_times[sorted_ids], points[sorted_ids]
            else:
                target_times, target_points = clean_target_times, points
            row = dict(case=case, condition=entry['condition'], variant=variant,
                       source_events=len(xy), search_events=len(points), patches=patch_keys.tolist(),
                       initializer=initializer, initializer_s=init_elapsed,
                       old_G_measurement_rpm=old[case].get('measurement_rpm'),
                       old_G_axis=old[case].get('measurement_axis'))
            if initializer['available']:
                search = methods.axis_search(target_points, target_times, initializer['omega_rad_s'])
                decision = methods.decision(target_points, target_times, search, 'S')
                row.update(search=compact_search(search), final=decision,
                           pipeline_s=init_elapsed + search['search_s'])
            else:
                row.update(search=None, final=None, pipeline_s=init_elapsed)
            if variant == 'clean' or variant == 'shuffled_seed_shuffled_refinement':
                if len(points) >= 100:
                    g = methods.g_search(target_points, target_times)
                    row['G_recomputed'] = dict(search=compact_search(g), final=methods.decision(target_points, target_times,g,'G'))
                else:
                    row['G_recomputed'] = dict(final=dict(measurement_valid=False, failure_reason='insufficient_events'))
            results.append(row)
            save('results.json', results)
            print(case, variant, 'init',initializer.get('rpm'), 'constraints',initializer['constraints'],
                  'final',row['final']['measurement_rpm'] if row['final'] else None, flush=True)
    save('verification.json', dict(source_hashes_unchanged=all(digest(p)==h for p,h in protocol['source_hashes'].items()),
                                   rows=len(results), expected_rows=10,
                                   no_direction_grid_in_flow_route=all(r['search'] is None or r['search']['nm_calls']==1 for r in results),
                                   reference_spin_axis_or_rpm_used=False, supplied_geometry_used=True, note='Old B1 geometry used; reference RPM and shaft not passed to fitting'))

if __name__ == '__main__':
    main()
