"""Offline paired geometry/window experiment; reference pairing stays explicit."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from time import perf_counter

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
spec = importlib.util.spec_from_file_location('phase2_readonly', OUT.parent/'drill_phase2_20261106/run.py')
p2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p2)
spec = importlib.util.spec_from_file_location('current_cmax', ROOT/'g1_spin/src/g1_spin/spherical_cmax.py')
cmax = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmax)
ENDS = {1: [1.435, 1.505, 3.605, 4.375], 4: [.905, 1.705, 2.505, 3.305],
        5: [4.105, 5.105, 6.105, 6.505]}
METHODS = ['fixed', 'linear', 'median3']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False, ensure_ascii=False))


def smooth(ts, xy):
    result = xy.copy()
    for i in range(1, len(ts)-1):
        if max(ts[i]-ts[i-1], ts[i+1]-ts[i]) <= 30000:
            result[i] = np.median(xy[i-1:i+2], axis=0)
    return result


def read_window(parts, origin, begin, end):
    chunks, ids, offset = [], [], 0
    for part in parts:
        lo, hi = np.searchsorted(part['t'], [origin+begin, origin+end])
        if hi > lo:
            chunks.append(part[lo:hi])
            ids.append(np.arange(offset+lo, offset+hi, dtype=np.int64))
        offset += len(part)
    if not chunks:
        return np.empty(0, dtype=parts[0].dtype), np.empty(0, dtype=np.int64)
    return np.concatenate(chunks), np.concatenate(ids)


def main():
    if (OUT/'results.jsonl').exists():
        raise FileExistsError('Preserve existing experiment')
    (OUT/'inputs').mkdir(exist_ok=True)
    protected = [OUT/'protocol.md', OUT/'run.py', ROOT/'g1_spin/src/g1_spin/spherical_cmax.py',
                 p2.OLD/'rotation_geometry.csv', p2.OLD/'control_geometry.csv', p2.OUT/'geometry.json',
                 ROOT/'experiments/wechat_spin_20261105/originals/results.json']
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in protected}
    reference = json.loads(protected[-1].read_text())
    save(OUT/'reference.json', {'status': 'accepted_approximate_video_reference',
        'event_recording_pairing': 'unconfirmed', 'frame_synchronization': 'unconfirmed',
        'videos': [{k:r[k] for k in ['video','requested_interval_s','marker_A_crossing_rpm',
                    'conservative_readout_range_hz']} for r in reference]})
    save(OUT/'freeze.json', hashes)
    rows = []
    with (OUT/'results.jsonl').open('x') as stream:
        for run, ends in ENDS.items():
            folder = p2.ARCHIVE/'decoded'/p2.NAMES[run]
            manifest = json.loads((folder/'manifest.json').read_text())
            parts = [np.load(folder/f, mmap_mode='r') for f in manifest['parts']]
            origin = int(parts[0]['t'][0])
            ts, xy = p2.observations(run)
            median_xy = smooth(ts, xy)
            jobs = [(end, 'rotation' if run != 5 else 'translation',
                     'development' if i < 2 else 'evaluation') for i,end in enumerate(ends)]
            if run == 5:
                jobs += [(.535,'static','development'), (1.535,'static','evaluation')]
            for end_s, condition, partition in jobs:
                end = round(end_s*1e6)
                for ms in [10,20,40]:
                    begin = end-ms*1000
                    events, ids = read_window(parts, origin, begin, end)
                    times = events['t'].astype(np.int64)-origin
                    pixels = np.column_stack([events['x'],events['y']]).astype(float)
                    name = f'run{run}_{end:08d}_{ms}ms'
                    base = dict(window=name, run=run, condition=condition, partition=partition,
                                begin_us=begin, end_us=end, window_ms=ms, raw_count=len(events),
                                common_count=0, optimizer_count=0, geometry_lookahead_ms=None)
                    try:
                        if condition == 'static':
                            centers = {m: np.tile([714.,355.],(len(events),1)) for m in METHODS}
                        else:
                            p2.centers(ts,xy,[begin,end])
                            center = p2.centers(ts,xy,[(begin+end)//2])[0]
                            centers = {'fixed':np.tile(center,(len(events),1)),
                                       'linear':p2.centers(ts,xy,times),
                                       'median3':p2.centers(ts,median_xy,times)}
                            right = min(int(np.searchsorted(ts,end))+1,len(ts)-1)
                            base['geometry_lookahead_ms'] = max(0.,(int(ts[right])+2500-end)/1000)
                        masks = [np.linalg.norm(pixels-centers[m],axis=1)<=48 for m in METHODS]
                        selected = np.flatnonzero(np.logical_and.reduce(masks))
                        common_ids = ids[selected]
                        base['common_count'] = len(selected)
                        if len(selected)>5000:
                            selected = selected[np.linspace(0,len(selected)-1,5000).astype(int)]
                        base['optimizer_count'] = len(selected)
                        np.savez_compressed(OUT/'inputs'/f'{name}.npz', ids=ids[selected],
                            common_ids=common_ids, events=events[selected], origin_us=origin,
                            **{m:centers[m][selected] for m in METHODS})
                        for method in METHODS:
                            t0 = perf_counter()
                            rel = pixels[selected]-centers[method][selected]
                            points = cmax.lift(rel[:,0],rel[:,1],0,0,60.)
                            result = cmax.estimate(points,(times[selected]-begin)*1e-6,max_rpm=3500)
                            elapsed = perf_counter()-t0
                            row = dict(**base, method=method, valid=result['valid'],
                                rpm=result['spin_rpm'], omega=result['omega_rad_s'],
                                reason=result['failure_reason'], elapsed_s=elapsed,
                                quality=result['quality'])
                            rows.append(row)
                            stream.write(json.dumps(row,allow_nan=False)+'\n'); stream.flush()
                    except ValueError as exc:
                        for method in METHODS:
                            row = dict(**base,method=method,valid=False,rpm=None,omega=None,
                                       reason=str(exc),elapsed_s=None,quality={})
                            rows.append(row);stream.write(json.dumps(row)+'\n');stream.flush()
                    print(name, base['optimizer_count'], 'events', flush=True)
            del parts
    assert len(rows)==126
    summary=[]
    for condition in ['rotation','translation','static']:
        for partition in ['development','evaluation']:
            for ms in [10,20,40]:
                for method in METHODS:
                    group=[r for r in rows if (r['condition'],r['partition'],r['window_ms'],r['method']) ==
                           (condition,partition,ms,method)]
                    values=[r['rpm'] for r in group if r['valid']]
                    elapsed=[r['elapsed_s'] for r in group if r['elapsed_s'] is not None]
                    summary.append(dict(condition=condition,partition=partition,window_ms=ms,method=method,
                        planned=len(group),valid=len(values),rpm_median=float(np.median(values)) if values else None,
                        rpm_min=min(values) if values else None,rpm_max=max(values) if values else None,
                        runtime_median_s=float(np.median(elapsed)) if elapsed else None))
    save(OUT/'summary.json',summary)
    with (OUT/'results.csv').open('w') as f:
        keys=[k for k in rows[0] if k not in ['quality','omega']]
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
        writer.writerows({k:r[k] for k in keys} for r in rows)
    assert hashes == {str(p.relative_to(ROOT)):digest(p) for p in protected}
    save(OUT/'verification.json',dict(result_rows=len(rows),protected_hashes_unchanged=True))


if __name__ == '__main__':
    main()
