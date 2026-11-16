"""Apply the frozen coarse initializer configuration to all108 historical windows."""
import hashlib
import importlib.util
import json
from pathlib import Path
from time import perf_counter
import numpy as np

OUT = Path(__file__).resolve().parent
PARENT = OUT.parent
spec = importlib.util.spec_from_file_location('coarse_initializer', PARENT / 'run.py')
coarse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coarse)
methods = coarse.methods


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    if (OUT / 'results.json').exists():
        raise FileExistsError('Preserve completed extension')
    old_rows = [json.loads(line) for line in (coarse.PHASE3 / 'results.jsonl').open()]
    old_rows = [row for row in old_rows if row['method'] == 'G']
    protocol = dict(cases=[row['window'] for row in old_rows], fixed_config=coarse.digest(PARENT / 'protocol.json'),
                    input='phase3 window selected events for both local normal fitting and CMAX; existing masks and max5000 sampling',
                    relation_to_pilot='Pilot normal fitting used complete phase4 common-mask inputs; extension uses selected phase3 inputs',
                    no_changes_to_parameters=True, purpose='Historical expansion to seek failures; not held-out recording evaluation',
                    source_hashes={str(PARENT / 'run.py'): coarse.digest(PARENT / 'run.py'),
                                   str(Path(__file__)): coarse.digest(Path(__file__)),
                                   str(coarse.PHASE3 / 'results.jsonl'): coarse.digest(coarse.PHASE3 / 'results.jsonl')})
    for old in old_rows:
        path = coarse.PHASE3 / 'windows' / (old['window'] + '.npz')
        if path.exists():
            protocol['source_hashes'][str(path)] = coarse.digest(path)
    save('protocol.json', protocol)
    rows = []
    for number, old in enumerate(old_rows):
        path = coarse.PHASE3 / 'windows' / (old['window'] + '.npz')
        row = dict(case=old['window'], condition=old['condition'], run=old['run'],
                   old_G_valid=old['measurement_valid'], old_G_rpm=old['measurement_rpm'],
                   old_G_axis=old['measurement_axis'], old_G_failure=old['failure_reason'],
                   source_file=str(path), initializer=None, final=None, search=None, status=None)
        if not path.exists():
            row['status'] = 'geometry_unavailable'
            rows.append(row)
            continue
        data = np.load(path)
        xy = data['xy']
        times = (data['time_us'] - old['start_us']) * 1e-6
        polarity = data['polarity']
        radius = float(data['radius_px'])
        center = data['center']
        row.update(input_events=len(xy), common_events=len(data['common_ids']),
                   downsampled=len(xy) < len(data['common_ids']), radius=radius, center=center.tolist())
        if len(xy) < 100:
            row['status'] = 'insufficient_events'
            rows.append(row)
            continue
        started = perf_counter()
        patches, counts = np.unique(np.floor(xy / 12).astype(int), axis=0, return_counts=True)
        order = np.lexsort((patches[:,1], patches[:,0], -counts))[:6]
        selected_patches = patches[order]
        initializer = coarse.normal_seed(xy, times, polarity, center, radius, selected_patches)
        row.update(initializer=initializer, patches=selected_patches.tolist())
        if initializer['available']:
            points = methods.cmax.lift(xy[:,0], xy[:,1], *center, radius)
            search = methods.axis_search(points, times, initializer['omega_rad_s'])
            final = methods.decision(points, times, search, 'S')
            row.update(search=coarse.compact_search(search), final=final,
                       status='accepted' if final['measurement_valid'] else 'refinement_rejected')
        else:
            row['status'] = 'initializer_rejected'
        row['compute_s'] = perf_counter() - started
        rows.append(row)
        if number % 10 == 0:
            print(number + 1, old['window'], row['status'], flush=True)
            save('progress.json', rows)
    save('results.json', rows)
    save('verification.json', dict(rows=len(rows), expected_rows=108,
                                   source_hashes_unchanged=all(coarse.digest(p)==h for p,h in protocol['source_hashes'].items()),
                                   single_NM_no_grid=all(r['search'] is None or r['search']['nm_calls']==1 for r in rows)))
    print('COMPLETE',len(rows),flush=True)

if __name__ == '__main__':
    main()
