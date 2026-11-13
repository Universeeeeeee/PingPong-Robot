import json
import numpy as np
from scipy.spatial.distance import cdist
from run import OUT,ROOT,FILES,load_file,digest,save
rows=json.loads((OUT/'profile_results.json').read_text());assert len(rows)==10
checks=[]
for row in rows:
 name=row['case'];path,start,du=next((p,s,d) for n,p,s,d in FILES if n==name);xy,t,pol,cs,radius=load_file(path,start);duration=du/1e6;cut=.6*duration
 v=row['delta_center_vx_px_s'];cs=cs.copy();cs[:,0]+=v*(t-cut);assert row['events_before']==row['events_after']==len(t)
 rel=(xy-cs)/radius;pts=np.column_stack([rel,-np.sqrt(1-np.sum(rel*rel,axis=1))]);source=(t<cut)&(t>=cut-min(.005,cut/3));omega=np.array(row['result']['omega_rad_s']);speed=np.linalg.norm(omega);axis=omega/speed;values=[]
 for k in range(4):
  lo=cut+k*duration*.1;hi=cut+(k+1)*duration*.1;mid=(lo+hi)/2;target=(t>=lo)&(t<hi);a=(mid-t[source])*speed
  moved=pts[source]*np.cos(a)[:,None]+np.cross(axis,pts[source])*np.sin(a)[:,None]+np.outer(pts[source]@axis,axis)*(1-np.cos(a))[:,None]
  pred=moved[:,:2]*radius+np.array([np.interp(mid,t,cs[:,j]) for j in range(2)]);visible=moved[:,2]<0;fw=[];bw=[]
  for p in [0,1]:
   sm=pol[source]==p;hm=pol[target]==p;px=pred[sm];seen=visible[sm];actual=xy[target][hm];d1=np.full(len(px),3.);d2=np.full(len(actual),3.)
   if seen.any() and len(actual):
    dist=cdist(px[seen],actual);d1[seen]=np.minimum(dist.min(axis=1),3);d2=np.minimum(dist.min(axis=0),3)
   fw.extend(d1);bw.extend(d2)
  values.append((np.mean(fw)+np.mean(bw))/2)
 error=abs(np.mean(values)-row['evaluation']['future_balanced_mean']);assert error<1e-9;checks.append(dict(case=name,velocity_offset=v,events_unchanged=True,metric_error=error))
for p,h in json.loads((OUT/'profile_freeze.json').read_text()).items():assert digest(ROOT/p)==h
save('profile_verification.json',dict(checks=checks,source_hashes_unchanged=True))
print('10 profile candidates verified against direct raw-pixel brute-force prediction')
