import json
import numpy as np
from run import OUT,ROOT,FILES,load_file,cmax,score,digest,save
from time import perf_counter
if (OUT/'profile_results.json').exists():raise FileExistsError('Preserve profile')
paths=[OUT/'profile.py',OUT/'profile_protocol.md',OUT/'run.py']+[p for name,p,_,_ in FILES if name in ['disk_motion','drill_1']]
hashes={str(p.relative_to(ROOT)):digest(p) for p in paths};save('profile_freeze.json',hashes)
rows=[]
for name,path,start,duration_us in FILES:
 if name not in ['disk_motion','drill_1']:continue
 xy,t,pol,cs,radius=load_file(path,start);duration=duration_us/1e6;cut=.6*duration
 centers={v:cs+np.column_stack([v*(t-cut),np.zeros(len(t))]) for v in [-80,-40,0,40,80]}
 keep=np.logical_and.reduce([np.linalg.norm(xy-c,axis=1)<.95*radius for c in centers.values()]);train=t[keep]<cut
 for v,c in centers.items():
  rel=(xy[keep]-c[keep])/radius;points=np.column_stack([rel,-np.sqrt(1-np.sum(rel*rel,axis=1))]);result=cmax.estimate(points,t[keep][train],max_rpm=3500) if False else cmax.estimate(points[train],t[keep][train],max_rpm=3500)
  omega=result['omega_rad_s'] if result['valid'] else result['quality'].get('rejected_candidate_omega_rad_s');evaluation=score(points,t[keep],pol[keep],c[keep],radius,np.asarray(omega),duration) if omega is not None else None
  rows.append(dict(case=name,delta_center_vx_px_s=v,events_before=len(t),events_after=int(keep.sum()),result=result,evaluation=evaluation));print(name,v,result['valid'],result['spin_rpm'],evaluation['future_balanced_mean'] if evaluation else None,flush=True)
save('profile_results.json',rows);assert all(digest(ROOT/p)==h for p,h in hashes.items())
