"""Paired old/new geometry with identical event IDs and unchanged CMAX."""
import importlib.util
import json
from pathlib import Path
from time import perf_counter
import numpy as np
OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('front',OUT/'run.py');front=importlib.util.module_from_spec(s);s.loader.exec_module(front)
base=front.base;old=base.old;p2=base.p2;cmax=base.prev.cmax
if (OUT/'spin_results.json').exists():raise FileExistsError('Preserve results')
rows=json.loads((OUT/'measurements.json').read_text());anchors=json.loads((front.SOURCE/'anchors.json').read_text())
prior=json.loads((OUT.parent/'drill_geometry_diagnosis_20261110/geometry_audit.json').read_text())
files=[OUT/'spin.py',OUT/'spin_protocol.md',OUT/'measurements.json',OUT/'contour.py',base.ROOT/'g1_spin/src/g1_spin/spherical_cmax.py']
hashes={str(p.relative_to(base.ROOT)):old.digest(p) for p in files};old.save(OUT/'spin_freeze.json',hashes)
(OUT/'spin_inputs').mkdir(exist_ok=True);results=[]
for run in [3,5]:
    folder=p2.ARCHIVE/'decoded'/p2.NAMES[run];manifest=json.loads((folder/'manifest.json').read_text())
    parts=[np.load(folder/f,mmap_mode='r') for f in manifest['parts']];origin=int(parts[0]['t'][0]);ots,oxy=p2.observations(run);oxy=old.smooth(ots,oxy)
    for a in [a for a in anchors if a['run']==run]:
        rr=[r for r in rows if r['case']==a['case'] and r['duration_us']==5000 and r['fits']['all']['valid']]
        nts=np.array([r['mid_us'] for r in rr],dtype=np.int64);nxy=np.array([r['fits']['all']['center'] for r in rr]);nxy=old.smooth(nts,nxy)
        end=a['end_us'];begin=end-20000;radius=a['radius']
        try:
            p2.centers(nts,nxy,[begin,end]);p2.centers(ots,oxy,[begin,end])
            raw,ids=old.read_window(parts,origin,begin,end);times=raw['t'].astype(np.int64)-origin;pix=np.column_stack([raw['x'],raw['y']]).astype(float)
            cs={'old_median3':p2.centers(ots,oxy,times),'outer_median3':p2.centers(nts,nxy,times)}
            mask=np.logical_and.reduce([np.linalg.norm(pix-c,axis=1)<=.8*radius for c in cs.values()]);ix=np.flatnonzero(mask);common=ids[ix]
            if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
            np.savez_compressed(OUT/'spin_inputs'/f"{a['case']}.npz",events=raw[ix],ids=ids[ix],common_ids=common,origin_us=origin,radius=radius,**{k:v[ix] for k,v in cs.items()})
            for method,center in cs.items():
                start=perf_counter();rel=pix[ix]-center[ix]
                result=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,radius),(times[ix]-begin)*1e-6,max_rpm=3500)
                results.append(dict(case=a['case'],run=run,method=method,result=result,elapsed_s=perf_counter()-start,input_count=len(ix)))
        except ValueError as e:
            for method in ['old_median3','outer_median3']:results.append(dict(case=a['case'],run=run,method=method,result=dict(valid=False,spin_rpm=None,axis_camera=None,omega_rad_s=None,failure_reason=str(e))))
        print(a['case'],flush=True)
    del parts
old.save(OUT/'spin_results.json',results)
assert len(results)==12 and all(old.digest(base.ROOT/p)==h for p,h in hashes.items())
print('12 paired states complete')
