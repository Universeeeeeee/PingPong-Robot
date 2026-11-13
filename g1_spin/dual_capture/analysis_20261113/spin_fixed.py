import json
import numpy as np
from scipy.optimize import minimize_scalar,least_squares
from audit_remote import OUT,save
import spherical_cmax as c

def scalar(p,t,axis):
 axis=np.array(axis,float);axis/=np.linalg.norm(axis)
 def score(rpm):return float(np.var(c.iwe(c.warp(p,t,axis*rpm*2*np.pi/60))/len(p)))
 grid=np.arange(-1200,1201,20);v=[score(x) for x in grid];sol=[]
 for i in range(1,len(grid)-1):
  if v[i]>=v[i-1] and v[i]>=v[i+1]:
   r=minimize_scalar(lambda x:-score(x),bounds=(grid[i-1],grid[i+1]),method='bounded');sol.append((-r.fun,r.x))
 sol.sort(reverse=True);s,r=sol[0];return dict(signed_rpm=r,rpm=abs(r),axis=axis.tolist(),gain=s/max(score(0),1e-30)-1,peaks=[dict(rpm=float(r),score=float(s)) for s,r in sol[:3]])

def main():
 tracks=json.loads((OUT/'event_track.json').read_text());out=[]
 for k in [17,19,31,42,45]:
  d=np.load(OUT/f'windows/t{k}.npz');e=d['events'];rel=(e['t'].astype(float)-int(d['origin']))/1e6;e=e[(rel>=k+.10)&(rel<k+.15)];t=(e['t'].astype(float)-int(d['origin']))/1e6
  tt=[];centers=[]
  for row in tracks:
   if k-.3<row['t']<k+.55:
    cc=[q for q in row['candidates'] if 310<q[1]<390 and 60<q[2]<82 and 300<q[0]<700]
    if len(cc)==1:tt.append(row['t']);centers.append(cc[0])
  if len(tt)<3:out.append(dict(time=k,error='insufficient_geometry'));continue
  centers=np.array(centers);cx=np.full(len(t),np.median(centers[:,0]));cy=np.full(len(t),np.median(centers[:,1]));r=float(np.median(centers[:,2]));x=e['x']-cx;y=e['y']-cy;sel=x*x+y*y<(.8*r)**2;ix=np.flatnonzero(sel);ix=ix[np.linspace(0,len(ix)-1,min(4000,len(ix))).astype(int)];dt=t[ix]-(k+.125);axis=[0,1,0] if k<40 else [.65,.76,0]
  p=c.lift(x[ix],y[ix],0,0,r);one=scalar(p,dt,axis);free=c.estimate(p,dt,max_rpm=1200)
  free['quality']={key:v for key,v in free['quality'].items() if key not in ['objective_histories','candidate_scan']}
  sensitivity=[dict(radius=r+dr,one=scalar(c.lift(x[ix],y[ix],0,0,r+dr),dt,axis)) for dr in [-2,2]]
  result=dict(time=k,events=int(sel.sum()),used=len(ix),geometry=dict(center=np.median(centers[:,:2],axis=0).tolist(),radius=r,track_samples=len(tt),center_range=np.ptp(centers[:,:2],axis=0).tolist()),one=one,free=free,sensitivity=sensitivity)
  out.append(result);save('spin_short_fixed_center.json',out);print(json.dumps(result),flush=True)
 save('spin_short_fixed_center.json',out)
if __name__=='__main__':main()
