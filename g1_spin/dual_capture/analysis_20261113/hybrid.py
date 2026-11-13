import json
import numpy as np
from audit_remote import OUT,save
from spin import scalar
import spherical_cmax as c
tracks=json.loads((OUT/'event_track.json').read_text());short=json.loads((OUT/'spin_short50ms.json').read_text());out=[]
for row in short:
 k=row['time'];axis=row['free']['axis_camera']
 if axis is None:continue
 d=np.load(OUT/f'windows/t{k}.npz');e=d['events'];t=(e['t'].astype(float)-int(d['origin']))/1e6;tt=[];centers=[]
 for q in tracks:
  if k-.3<q['t']<k+.55:
   cc=[v for v in q['candidates'] if 310<v[1]<390 and 60<v[2]<82 and 300<v[0]<700]
   if len(cc)==1:tt.append(q['t']);centers.append(cc[0])
 centers=np.array(centers);x=e['x']-np.interp(t,tt,centers[:,0]);y=e['y']-np.interp(t,tt,centers[:,1]);r=np.median(centers[:,2]);ix=np.flatnonzero(x*x+y*y<(.8*r)**2);ix=ix[np.linspace(0,len(ix)-1,min(4000,len(ix))).astype(int)]
 result=scalar(c.lift(x[ix],y[ix],0,0,r),t[ix]-(k+.125),axis);out.append(dict(time=k,result=result,axis_source='event-only free 3D 50ms window inside same 250ms window',note='Exploratory reuse of same data; not held-out validation'))
save('hybrid_exploratory.json',out);print(json.dumps(out,indent=2))
