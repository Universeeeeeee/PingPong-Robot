"""Centered 10ms comparison; same geometry as existing 50/250ms experiments."""
import json
import numpy as np
from audit_remote import OUT,save
import spherical_cmax as c
tracks=json.loads((OUT/'event_track.json').read_text());out=[]
for k in [17,19,31,42,45]:
 d=np.load(OUT/f'windows/t{k}.npz');e=d['events'];rel=(e['t'].astype(float)-int(d['origin']))/1e6;e=e[(rel>=k+.120)&(rel<k+.130)];t=(e['t'].astype(float)-int(d['origin']))/1e6
 tt=[];centers=[]
 for row in tracks:
  if k-.3<row['t']<k+.55:
   cc=[q for q in row['candidates'] if 310<q[1]<390 and 60<q[2]<82 and 300<q[0]<700]
   if len(cc)==1:tt.append(row['t']);centers.append(cc[0])
 if len(tt)<3:
  out.append(dict(time=k,error='insufficient_geometry'));continue
 centers=np.array(centers);x=e['x']-np.interp(t,tt,centers[:,0]);y=e['y']-np.interp(t,tt,centers[:,1]);r=float(np.median(centers[:,2]));ix=np.flatnonzero(x*x+y*y<(.8*r)**2);count=len(ix)
 if count<100:
  out.append(dict(time=k,error='insufficient_events',events=count));continue
 ix=ix[np.linspace(0,len(ix)-1,min(4000,len(ix))).astype(int)];p=c.lift(x[ix],y[ix],0,0,r);dt=t[ix]-(k+.125)
 for bound in [1200,3500]:
  free=c.estimate(p,dt,max_rpm=bound)
  free['quality']={key:v for key,v in free['quality'].items() if key not in ['objective_histories','candidate_scan']}
  row=dict(time=k,bound_rpm=bound,start_s=k+.120,end_s=k+.130,events=count,used=len(ix),radius=r,free=free);out.append(row);save('spin_10ms.json',out);print(json.dumps(row),flush=True)
save('spin_10ms.json',out)
