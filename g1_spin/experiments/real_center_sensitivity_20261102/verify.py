"""Independent closed-form and absolute-time boolean reconstruction; no optimization."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
DATA = ROOT/'g1_spin/data/decoded/local_analysis_222403'
models = json.loads((OUT/'models.json').read_text())
results = list(csv.DictReader((OUT/'results.csv').open()))
audits = list(csv.DictReader((OUT/'mask_audit.csv').open()))
with (OUT.parent/'geometry_observation_20261101/observations.csv').open() as f:
    observations = list(csv.DictReader(f))
max_coefficient_difference = 0.
for segment, durations in models.items():
    for duration, model in durations.items():
        rows = [r for r in observations if r['segment']==segment and r['duration_ms']==duration
                and r['mode']=='fixed48' and r['threshold_policy']=='fixed_0.6']
        t = np.array([float(r['midpoint_s']) for r in rows])
        xy = np.array([[float(r['cx']), float(r['cy'])] for r in rows])
        mean_t, mean_xy = t.mean(), xy.mean(axis=0)
        v = ((t-mean_t)[:, None]*(xy-mean_xy)).sum(axis=0)/np.sum((t-mean_t)**2)
        anchor_xy = mean_xy+v*(model['anchor_s']-mean_t)
        difference = max(np.max(np.abs(v-model['velocity_px_s'])),
                         np.max(np.abs(anchor_xy-model['center_anchor_px'])))
        max_coefficient_difference = max(max_coefficient_difference, float(difference))
assert max_coefficient_difference < 1e-10
manifest = json.loads((DATA/'manifest.json').read_text())
events = np.concatenate([np.load(DATA/p) for p in manifest['parts']])
origin = int(events['t'][0])
global_ids = np.arange(len(events), dtype=np.int64)
windows = []
for window in dict.fromkeys(r['window'] for r in results):
    rr = [r for r in results if r['window']==window]
    first = rr[0]; segment = first['segment']
    lo, hi, ref = [round(float(first[k])*1e6) for k in ['start_s', 'end_s', 'reference_s']]
    # Absolute timestamp boolean mask differs from run.py's relative searchsorted path.
    selected_time = (events['t']>=origin+lo)&(events['t']<origin+hi)
    ev, ids = events[selected_time], global_ids[selected_time]
    t = (ev['t'].astype(np.int64)-origin)/1e6
    xy = np.column_stack((ev['x'], ev['y'])).astype(float)
    roi = (xy[:, 0]>=475)&(xy[:, 0]<638)&(xy[:, 1]>=143)&(xy[:, 1]<285)
    common = roi&(xy[:, 1]<240)
    curves = {}
    for row in rr:
        variant = row['variant']; model = models[segment]['100' if variant=='C_100ms' else '50']
        c, v = np.array(model['center_anchor_px']), np.array(model['velocity_px_s'])
        anchor = model['anchor_s']
        ref_c = c+v*(ref/1e6-anchor)
        if variant=='A_fixed':
            center = np.tile(ref_c, (len(t), 1))
        elif variant.startswith('B_velocity_'):
            factor = .75 if variant.endswith('075') else 1.25
            center = ref_c+factor*(t-ref/1e6)[:, None]*v
        else:
            center = c+(t-anchor)[:, None]*v
            if variant.startswith('B_cx_'):
                center[:, 0] += 2 if 'plus' in variant else -2
            if variant.startswith('B_cy_'):
                center[:, 1] += 2 if 'plus' in variant else -2
        curves[variant] = center
        squared_distance = np.sum((xy-center)**2, axis=1)
        common &= squared_distance <= (48*.8)**2
        audit = next(a for a in audits if a['window']==window and a['variant']==variant)
        assert int(audit['outside_sphere_roi']) == int((roi&(squared_distance>48**2)).sum())
    saved = np.load(OUT/f'{window}_events.npz')
    common_ids = ids[common]
    expected_ids = common_ids if len(common_ids)<=5000 else common_ids[np.linspace(0, len(common_ids)-1, 5000).astype(int)]
    assert np.array_equal(ids, saved['raw_global_event_ids'])
    assert np.array_equal(common_ids, saved['common_global_event_ids'])
    assert np.array_equal(expected_ids, saved['optimizer_global_event_ids'])
    assert np.array_equal(events['t'][expected_ids].astype(np.int64)-origin, saved['relative_time_us'])
    assert np.array_equal(events['p'][expected_ids], saved['polarity'])
    assert np.array_equal(events['x'][expected_ids], saved['x'])
    assert np.array_equal(events['y'][expected_ids], saved['y'])
    expected_hash = hashlib.sha256(expected_ids.astype('<i8').tobytes()).hexdigest()
    assert {r['event_id_sha256'] for r in rr} == {expected_hash}
    assert len(rr)==9 and all(r['physics_validated']=='False' for r in rr)
    for row in rr:
        omega = np.array([float(row[k]) for k in ['wx_rad_s', 'wy_rad_s', 'wz_rad_s']])
        for base in ['A_fixed', 'B_50ms']:
            br = next(x for x in rr if x['variant']==base)
            bw = np.array([float(br[k]) for k in ['wx_rad_s', 'wy_rad_s', 'wz_rad_s']])
            # atan2 cross/dot independently verifies the arccos angle used by run.py.
            degree = float(np.rad2deg(np.arctan2(np.linalg.norm(np.cross(omega, bw)), np.dot(omega, bw))))
            assert abs(degree-float(row[f'directed_angle_to_{base}_deg'])) < 2e-6
    windows.append(dict(window=window, raw_events=len(ids), common_events=len(common_ids),
                        optimizer_events=len(expected_ids), reconstructed_ids_exact=True))
environment = json.loads((OUT/'environment.json').read_text())
for path, expected in environment['protected_input_sha256'].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected
assert len(results)==54
report = dict(independent_ols_coefficient_max_difference=max_coefficient_difference,
    windows=windows, absolute_boolean_vs_searchsorted_exact=True,
    independently_reconstructed_common_mask_exact=True, timestamps_polarity_pixels_exact=True,
    directed_angle_atan2_vs_arccos_passed=True, protected_inputs_unchanged=True,
    physics_validated=False, verification_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'independent_validation.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
