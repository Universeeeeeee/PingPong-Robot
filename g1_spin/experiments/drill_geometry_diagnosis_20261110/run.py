"""Bounded 2x3 diagnosis without changing estimator or old experiment outputs."""
import importlib.util
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from route import covered,select

OUT=Path(__file__).resolve().parent
PREV=OUT.parent/'drill_candidate_validation_20261110'
s=importlib.util.spec_from_file_location('previous_trial',PREV/'run.py');prev=importlib.util.module_from_spec(s);s.loader.exec_module(prev)
ROOT=prev.ROOT;old=prev.old;p2=prev.p2;cmax=prev.cmax
ENDS={1:[3565000],3:[1465000,2665000,4335000,6265000],4:[865000],5:[4135000,6535000]}
METHODS=['fixed','linear','median3']


def main():
    if (OUT/'results.jsonl').exists():raise FileExistsError('Preserve output')
    protected=[OUT/'protocol.md',OUT/'route.py',OUT/'run.py',PREV/'run.py',PREV/'results.jsonl',PREV/'baseline_supplement.json',
               old.OUT/'run.py',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py',p2.OUT/'run.py',p2.OUT/'geometry.json',
               p2.OLD/'rotation_geometry.csv',p2.OLD/'control_geometry.csv']
    hashes={str(p.relative_to(ROOT)):old.digest(p) for p in protected};old.save(OUT/'freeze.json',hashes)
    (OUT/'inputs').mkdir(exist_ok=True)
    # Route is decided solely from geometric timestamp coverage before any CMAX execution.
    routes=[]
    previous_rows=[json.loads(line) for line in (PREV/'results.jsonl').read_text().splitlines()]
    supplements=json.loads((PREV/'baseline_supplement.json').read_text())
    for j in json.loads((PREV/'windows.json').read_text()):
        ts,xy=p2.observations(j['run']);choice='20ms_median3' if j['condition']=='static' else select(ts,j['end_us'])
        name=f"run{j['run']}_{j['end_us']:08d}";source=None
        if choice=='20ms_median3':source=next(r for r in previous_rows if r['window']==name and r['method']=='candidate')
        elif choice=='10ms_fixed_fallback':
            source=next((r for r in supplements if r['window']==name),None)
        routes.append(dict(**j,window=name,route=choice,source_file='baseline_supplement.json' if choice=='10ms_fixed_fallback' else 'results.jsonl' if source else None,
                           reused_result=source['result'] if source else None))
    old.save(OUT/'routing.json',routes)
    audits=[]
    with (OUT/'results.jsonl').open('x') as stream:
        for run,ends in ENDS.items():
            folder=p2.ARCHIVE/'decoded'/p2.NAMES[run];manifest=json.loads((folder/'manifest.json').read_text())
            parts=[np.load(folder/f,mmap_mode='r') for f in manifest['parts']];origin=int(parts[0]['t'][0])
            ts,xy=p2.observations(run);smoothed=old.smooth(ts,xy);radius=prev.RADIUS[run]
            for end in ends:
                raw20,ids20=old.read_window(parts,origin,end-20000,end)
                projection={}
                if run==3:
                    center=p2.centers(ts,xy,[end-10000])[0]
                    for key,lo,hi in [('first10',end-20000,end-10000),('last10',end-10000,end),('full20',end-20000,end)]:
                        q=raw20[(raw20['t']>=origin+lo)&(raw20['t']<origin+hi)]
                        projection[key]=prev.shaft(q,center,radius,[np.cos(np.deg2rad(35)),np.sin(np.deg2rad(35))])
                nearby=(ts>=end-60000)&(ts<=end+40000)
                audit=dict(run=run,end_us=end,geometry_times_us=ts[nearby].tolist(),geometry_xy=xy[nearby].tolist(),
                           route=select(ts,end),shaft=projection)
                audits.append(audit)
                for duration in [10,20]:
                    begin=end-duration*1000;name=f'run{run}_{end:08d}_{duration}ms'
                    if not covered(ts,begin,end):
                        for method in METHODS:stream.write(json.dumps(dict(run=run,end_us=end,window_ms=duration,method=method,window=name,
                            result=dict(valid=False,spin_rpm=None,axis_camera=None,omega_rad_s=None,failure_reason='geometry_gap_or_coverage'),elapsed_s=None))+'\n')
                        stream.flush();continue
                    raw,ids=old.read_window(parts,origin,begin,end);t=raw['t'].astype(np.int64)-origin
                    pix=np.column_stack([raw['x'],raw['y']]).astype(float)
                    fixed=p2.centers(ts,xy,[(begin+end)//2])[0]
                    centers=dict(fixed=np.tile(fixed,(len(raw),1)),linear=p2.centers(ts,xy,t),median3=p2.centers(ts,smoothed,t))
                    masks=[np.linalg.norm(pix-centers[m],axis=1)<=.8*radius for m in METHODS]
                    ix=np.flatnonzero(np.logical_and.reduce(masks));common=ids[ix]
                    if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
                    np.savez_compressed(OUT/'inputs'/f'{name}.npz',ids=ids[ix],common_ids=common,events=raw[ix],origin_us=origin,radius=radius,
                                        **{m:centers[m][ix] for m in METHODS})
                    for method in METHODS:
                        start=perf_counter();rel=pix[ix]-centers[method][ix]
                        points=cmax.lift(rel[:,0],rel[:,1],0,0,radius)
                        result=cmax.estimate(points,(t[ix]-begin)*1e-6,max_rpm=3500)
                        r=dict(run=run,end_us=end,window_ms=duration,method=method,window=name,result=result,elapsed_s=perf_counter()-start,
                               raw_count=len(raw),common_count=len(common),input_count=len(ix))
                        stream.write(json.dumps(r,allow_nan=False)+'\n');stream.flush()
                    print(name,len(ix),'events',flush=True)
            del parts
    old.save(OUT/'geometry_audit.json',audits)
    assert all(old.digest(ROOT/p)==h for p,h in hashes.items())
    assert len((OUT/'results.jsonl').read_text().splitlines())==48
    print('Completed 48 states; frozen source files unchanged')

if __name__=='__main__':main()
