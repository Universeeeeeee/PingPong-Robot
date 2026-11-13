"""Frozen candidate vs baseline on disjoint time intervals, same recordings."""
import importlib.util
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.optimize import least_squares

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
PREV=OUT.parent/'drill_window_geometry_20261110'
spec=importlib.util.spec_from_file_location('previous_helpers',PREV/'run.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
p2=old.p2;cmax=old.cmax
BLOCKS={1:[.7,2.1,3.5,5.6],2:[.8,2.,3.4,5.2],3:[1.4,2.6,4.2,6.2],4:[.8,1.6,2.4,3.2],5:[4.,4.4,5.6,6.4]}
RADIUS={1:60.,2:63.,3:59.,4:60.,5:60.}

def shaft(raw,center,radius,direction):
    d=np.asarray(direction,float);d/=np.linalg.norm(d);p=np.array([-d[1],d[0]])
    xy=np.column_stack([raw['x'],raw['y']])-center;a=xy@d;b=xy@p;pts=[]
    for lo in np.arange(radius+15,radius+170,12):
        m=(a>=lo)&(a<lo+12)&(abs(b)<25)
        if m.sum()>=3:pts.append([np.median(a[m]),np.median(b[m])])
    if len(pts)<6:return dict(valid=False,reason='insufficient_support')
    pts=np.array(pts)
    if np.ptp(pts[:,0])<70:return dict(valid=False,reason='short_span')
    fit=least_squares(lambda x:x[0]*pts[:,0]+x[1]-pts[:,1],[0.,0.],loss='soft_l1')
    rms=float(np.sqrt(np.mean((fit.x[0]*pts[:,0]+fit.x[1]-pts[:,1])**2)))
    direction=d+fit.x[0]*p;direction/=np.linalg.norm(direction)
    return dict(valid=rms<=3,rms_px=rms,direction=direction.tolist(),reason=None if rms<=3 else 'large_residual')

def main():
    if (OUT/'results.jsonl').exists():raise FileExistsError('Preserve prior output')
    jobs=[]
    previous=json.loads((OUT.parent/'drill_phase3_20261106/window_audit.json').read_text())
    previous += [dict(run=r['run'],start_us=r['begin_us'],end_us=r['end_us']) for r in [json.loads(s) for s in (PREV/'results.jsonl').read_text().splitlines()]]
    for run,blocks in BLOCKS.items():
        for block in blocks:
            for offset in [.065,.135]:
                end=round((block+offset)*1e6);begin=end-20000
                assert not any(r['run']==run and max(begin,r['start_us'])<min(end,r['end_us']) for r in previous)
                jobs.append(dict(run=run,begin_us=begin,end_us=end,condition='rotation' if run!=5 else 'translation'))
    for end in [.665,1.665,2.665,3.165]:
        begin=round(end*1e6)-20000;end=round(end*1e6)
        assert not any(r['run']==5 and max(begin,r['start_us'])<min(end,r['end_us']) for r in previous)
        jobs.append(dict(run=5,begin_us=begin,end_us=end,condition='static'))
    assert len(jobs)==44
    old.save(OUT/'windows.json',jobs)
    files=[OUT/'run.py',OUT/'protocol.md',OUT/'windows.json',PREV/'run.py',PREV/'reference.json',
           ROOT/'g1_spin/src/g1_spin/spherical_cmax.py',p2.OUT/'run.py',p2.OUT/'geometry.json',
           p2.OLD/'rotation_geometry.csv',p2.OLD/'control_geometry.csv']
    hashes={str(p.relative_to(ROOT)):old.digest(p) for p in files}
    old.save(OUT/'freeze.json',hashes)
    old.save(OUT/'reference.json',json.loads((PREV/'reference.json').read_text()))
    (OUT/'inputs').mkdir(exist_ok=True)
    with (OUT/'results.jsonl').open('x') as stream:
        for run in BLOCKS:
            folder=p2.ARCHIVE/'decoded'/p2.NAMES[run]
            manifest=json.loads((folder/'manifest.json').read_text())
            parts=[np.load(folder/f,mmap_mode='r') for f in manifest['parts']]
            origin=int(parts[0]['t'][0]);ts,xy=p2.observations(run);smoothed=old.smooth(ts,xy);radius=RADIUS[run]
            for job in [j for j in jobs if j['run']==run]:
                begin,end=job['begin_us'],job['end_us'];name=f'run{run}_{end:08d}'
                raw,ids=old.read_window(parts,origin,begin,end)
                times=raw['t'].astype(np.int64)-origin;pixels=np.column_stack([raw['x'],raw['y']]).astype(float)
                try:
                    if job['condition']=='static':
                        fixed=np.array([714.,355.]);candidate=np.tile(fixed,(len(raw),1));projection=dict(valid=False,reason='static')
                    else:
                        p2.centers(ts,xy,[begin,end])
                        fixed=p2.centers(ts,xy,[end-5000])[0]
                        candidate=p2.centers(ts,smoothed,times)
                        center=p2.centers(ts,xy,[end-10000])[0]
                        projection=shaft(raw,center,radius,[0,1] if run in [1,4] else [np.cos(np.deg2rad(35)),np.sin(np.deg2rad(35))]) if run in [1,3,4] else dict(valid=False,reason='no_shaft_prior')
                    common=(np.linalg.norm(pixels-fixed,axis=1)<=.8*radius)&(np.linalg.norm(pixels-candidate,axis=1)<=.8*radius)
                    for method in ['baseline','candidate']:
                        mask=common & (times>=end-(10000 if method=='baseline' else 20000))
                        indices=np.flatnonzero(mask);count=len(indices)
                        if count>5000:indices=indices[np.linspace(0,count-1,5000).astype(int)]
                        centers=np.tile(fixed,(len(indices),1)) if method=='baseline' else candidate[indices]
                        np.savez_compressed(OUT/'inputs'/f'{name}_{method}.npz',ids=ids[indices],events=raw[indices],centers=centers,radius=radius,origin_us=origin)
                        start=perf_counter();rel=pixels[indices]-centers
                        points=cmax.lift(rel[:,0],rel[:,1],0,0,radius)
                        result=cmax.estimate(points,(times[indices]-begin)*1e-6,max_rpm=3500)
                        r=dict(**job,window=name,method=method,radius=radius,common_count=count,event_count=len(indices),
                               elapsed_s=perf_counter()-start,shaft=projection,result=result)
                        stream.write(json.dumps(r,allow_nan=False)+'\n');stream.flush()
                except ValueError as exc:
                    for method in ['baseline','candidate']:
                        stream.write(json.dumps(dict(**job,window=name,method=method,radius=radius,common_count=0,event_count=0,
                            elapsed_s=None,shaft=dict(valid=False,reason='geometry_unavailable'),result=dict(valid=False,spin_rpm=None,
                            omega_rad_s=None,axis_camera=None,failure_reason=str(exc))))+'\n');stream.flush()
                print(name,job['condition'],flush=True)
            del parts
    assert all(old.digest(ROOT/p)==h for p,h in hashes.items())
    print('DONE 88 rows; frozen files unchanged')

if __name__=='__main__':main()
