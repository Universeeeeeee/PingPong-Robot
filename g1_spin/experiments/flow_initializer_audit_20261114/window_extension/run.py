"""Same-center10/50/250ms optimizer comparison on decoded old recordings."""
import hashlib
import importlib.util
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize

OUT = Path(__file__).resolve().parent
AUDIT = OUT.parent
EXPERIMENTS = AUDIT.parent
ROOT = EXPERIMENTS.parents[1]
ARCHIVE = ROOT / 'g1_spin/data/archive_20261105/decoded'
spec = importlib.util.spec_from_file_location('coarse_initializer_for_windows', AUDIT / 'coarse_seed/run.py')
coarse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coarse)
cmax = coarse.methods.cmax
CASES = ['run1_01495000', 'run2_02095000', 'run3_01495000', 'run4_02495000']
NAMES = {1:'DSP03F100025_20261105_201058', 2:'DSP03F100025_20261105_201123',
         3:'DSP03F100025_20261105_201143', 4:'DSP03F100025_20261105_203829'}
BOUND = 1200 * 2*np.pi/60


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def extract_raw(run, middle_us, center, radius):
    folder = ARCHIVE / NAMES[run]
    manifest = json.loads((folder / 'manifest.json').read_text())
    first = np.load(folder / manifest['parts'][0], mmap_mode='r')
    origin = int(first['t'][0])
    low, high = origin+middle_us-125000, origin+middle_us+125000
    pieces, sources = [], {str(folder/'manifest.json'):sha256(folder/'manifest.json')}
    for name in manifest['parts']:
        path = folder / name
        events = np.load(path, mmap_mode='r')
        if int(events['t'][-1]) < low or int(events['t'][0]) >= high:
            continue
        begin,end = np.searchsorted(events['t'],[low,high])
        selected = np.array(events[begin:end])
        inside = (selected['x'].astype(float)-center[0])**2 + (selected['y'].astype(float)-center[1])**2 < (.8*radius)**2
        pieces.append(selected[inside])
        sources[str(path)] = sha256(path)
    return np.concatenate(pieces), origin, sources


def quality(points, times, omega, score, zero, converged):
    speed = float(np.linalg.norm(omega))
    gain = float(score/max(zero,1e-12)-1)
    reason = None
    if speed*np.ptp(times) < .03 or gain < .08:
        reason = 'no_rotation_evidence'
    elif not converged:
        reason = 'optimization_not_converged'
    elif speed >= BOUND*.985:
        reason = 'search_boundary'
    elif np.linalg.eigvalsh(np.cov(points.T))[-2] < 1e-5:
        reason = 'spatially_degenerate_features'
    return dict(valid=reason is None, reason=reason, gain=gain,
                omega_rad_s=np.asarray(omega).tolist(), rpm=speed*60/(2*np.pi),
                axis=(np.asarray(omega)/speed).tolist() if speed>1e-12 else None,
                contrast=score, zero_contrast=zero)


def single_search(points, times, initial_axis):
    dt = times-(times[0]+times[-1])/2
    contrast_calls = 0
    objective_calls = 0
    def contrast(omega):
        nonlocal contrast_calls
        contrast_calls += 1
        return float(np.var(cmax.iwe(cmax.warp(points,dt,omega))/len(points)))
    started = perf_counter()
    zero = contrast(np.zeros(3))
    scan = []
    seed = np.zeros(3)
    if initial_axis is not None:
        axis = np.asarray(initial_axis)/np.linalg.norm(initial_axis)
        n_speed = int(np.clip(np.ceil(BOUND*np.ptp(times)/.4),8,30))
        scan.append(dict(rpm=0.,score=zero))
        for magnitude in np.linspace(BOUND/n_speed,BOUND,n_speed):
            scan.append(dict(rpm=float(magnitude*60/(2*np.pi)),score=contrast(axis*magnitude)))
        best = max(scan,key=lambda x:x['score'])
        seed = axis*best['rpm']*2*np.pi/60
    def objective(omega):
        nonlocal objective_calls
        objective_calls += 1
        speed = np.linalg.norm(omega)
        if speed > BOUND:
            return 1+(speed/BOUND-1)**2
        return -contrast(omega)/max(zero,1e-12)
    simplex = np.vstack([seed,seed+np.eye(3)*max(2.,BOUND/40)])
    fit = minimize(objective,seed,method='Nelder-Mead',options=dict(initial_simplex=simplex,maxiter=350,xatol=.05,fatol=1e-6))
    score = contrast(fit.x)
    outcome = quality(points,times,fit.x,score,zero,bool(fit.success))
    return dict(outcome=outcome,selected_seed=seed.tolist(),zero_scan_won=bool(scan and np.linalg.norm(seed)<1e-12),
                scan=scan,nm_calls=1,objective_calls=objective_calls,contrast_calls=contrast_calls,
                iterations=int(fit.nit),search_s=perf_counter()-started)


def grid_search(points,times):
    original_minimize, original_iwe = cmax.minimize,cmax.iwe
    objective_calls,contrast_calls,nm_calls = 0,0,0
    def counted_iwe(*args,**kwargs):
        nonlocal contrast_calls
        contrast_calls += 1
        return original_iwe(*args,**kwargs)
    def counted_minimize(fun,*args,**kwargs):
        nonlocal nm_calls
        nm_calls += 1
        def objective(omega):
            nonlocal objective_calls
            objective_calls += 1
            return fun(omega)
        return original_minimize(objective,*args,**kwargs)
    cmax.minimize,cmax.iwe = counted_minimize,counted_iwe
    started = perf_counter()
    try:
        result = cmax.estimate(points,times,max_rpm=1200)
    finally:
        cmax.minimize,cmax.iwe = original_minimize,original_iwe
    elapsed = perf_counter()-started
    q = result['quality']
    omega = result['omega_rad_s'] if result['valid'] else q.get('rejected_candidate_omega_rad_s')
    outcome = quality(points,times,omega,q['contrast_after'],q['contrast_before'],q['converged'])
    outcome['native_valid'] = result['valid']
    outcome['native_reason'] = result['failure_reason']
    return dict(outcome=outcome,selected_seed=q['initial_omega_rad_s'],nm_calls=nm_calls,
                objective_calls=objective_calls,contrast_calls=contrast_calls,search_s=elapsed,
                competing_solution_displacement_rad=q['competing_solution_displacement_rad'])


def main():
    if (OUT/'results.json').exists():
        raise FileExistsError('Preserve completed window comparison')
    cases = {r['window']:r for r in (json.loads(line) for line in (coarse.PHASE3/'results.jsonl').open()) if r['method']=='G'}
    protocol = dict(cases=CASES,window_ms=[10,50,250],fixed_geometry='old B1 center and radius held fixed over entire250ms',
                    roi_ratio=.8,cmax_max_events=4000,max_rpm=1200,reference_time='actual sampled endpoint midpoint for every method',
                    coarse_seed='fixed top6 12px patches; five2ms segments on central10ms only; full ROI events',
                    search='flow axis magnitude scan then free3D NM on each10/50/250ms window',
                    native_grid_rival_gate='reported separately from common quality gates',
                    interpretation='optimizer comparison under supplied fixed geometry; no physical precision claim',
                    source_hashes={str(Path(__file__)):sha256(Path(__file__)),str(AUDIT/'coarse_seed/run.py'):sha256(AUDIT/'coarse_seed/run.py'),str(coarse.MODELS):sha256(coarse.MODELS),str(coarse.PHASE3/'frozen/spherical_cmax.py'):sha256(coarse.PHASE3/'frozen/spherical_cmax.py')})
    save('protocol.json',protocol)
    rows = []
    for case in CASES:
        old = cases[case]
        window_path = coarse.PHASE3/'windows'/f'{case}.npz'
        geometry = np.load(window_path)
        center,radius = geometry['center'],float(geometry['radius_px'])
        middle_us = (old['start_us']+old['end_us'])//2
        events,origin,sources = extract_raw(old['run'],middle_us,center,radius)
        protocol['source_hashes'].update(sources)
        protocol['source_hashes'][str(window_path)] = sha256(window_path)
        relative_us = events['t'].astype(np.int64)-origin-middle_us
        central = (relative_us>=-5000)&(relative_us<5000)
        xy = np.c_[events['x'][central],events['y'][central]].astype(float)
        times = relative_us[central]*1e-6+.005
        patches,counts = np.unique(np.floor(xy/12).astype(int),axis=0,return_counts=True)
        order = np.lexsort((patches[:,1],patches[:,0],-counts))[:6]
        started = perf_counter()
        initializer = coarse.normal_seed(xy,times,events['p'][central],center,radius,patches[order])
        initializer_s = perf_counter()-started
        for duration_ms in [10,50,250]:
            chosen = np.flatnonzero((relative_us>=-duration_ms*500)&(relative_us<duration_ms*500))
            selected = chosen[np.linspace(0,len(chosen)-1,min(4000,len(chosen))).astype(int)]
            points = cmax.lift(events['x'][selected],events['y'][selected],*center,radius)
            sample_times = relative_us[selected]*1e-6
            row = dict(case=case,run=old['run'],middle_us=middle_us,window_ms=duration_ms,
                       full_roi_events=len(chosen),used_events=len(selected),center=center.tolist(),radius=radius,
                       central10ms_events=int(central.sum()),initializer=initializer,initializer_s=initializer_s,
                       sampled_events_sha256=hashlib.sha256(events[selected].tobytes()).hexdigest())
            row['grid'] = grid_search(points,sample_times)
            row['zero'] = single_search(points,sample_times,None)
            row['flow'] = single_search(points,sample_times,initializer['axis']) if initializer['available'] else None
            rows.append(row)
            save('results.json',rows)
            print(case,duration_ms,[(name,row[name]['outcome']['rpm'],row[name]['outcome']['valid']) for name in ['grid','zero','flow'] if row[name]],flush=True)
        save('protocol.json',protocol)
    save('verification.json',dict(rows=len(rows),expected_rows=12,source_hashes_unchanged=all(sha256(p)==h for p,h in protocol['source_hashes'].items()),
                                   single_routes_one_NM=all(r['zero']['nm_calls']==1 and (r['flow'] is None or r['flow']['nm_calls']==1) for r in rows)))

if __name__ == '__main__':
    main()
