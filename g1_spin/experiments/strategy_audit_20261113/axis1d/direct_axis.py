"""Independent image-plane drill-axis annotations for 1D CMAX."""
import json,hashlib,importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3];PHASE=ROOT/'g1_spin/experiments/drill_phase3_20261106';SHAFT=PHASE/'shaft_projection.json'
s=importlib.util.spec_from_file_location('cmax',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py');cmax=importlib.util.module_from_spec(s);s.loader.exec_module(cmax)

def dig(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def scalar(points,t,a,bound=3500):
 a=np.asarray(a,float);a/=np.linalg.norm(a)
 def score(s):return float(np.var(cmax.iwe(cmax.warp(points,t,a*s))/len(points)))
 grid=np.linspace(-bound*2*np.pi/60,bound*2*np.pi/60,281);v=np.array([score(x) for x in grid]);peaks=[]
 for i in range(1,len(grid)-1):
  if v[i]>=v[i-1] and v[i]>=v[i+1]:
   q=minimize_scalar(lambda s:-score(s),bounds=(grid[i-1],grid[i+1]),method='bounded');peaks.append((float(-q.fun),float(q.x)))
 peaks.sort(reverse=True);c,s=peaks[0] if peaks else (score(0),0.)
 return dict(axis_camera=a.tolist(),valid=abs(s)>1e-12,omega_rad_s=(a*s).tolist(),spin_rpm=abs(s)*60/(2*np.pi),signed_rpm=s*60/(2*np.pi),contrast=c,zero_contrast=score(0),relative_gain=c/max(score(0),1e-30)-1,peaks=[dict(signed_rpm=x*60/(2*np.pi),contrast=y) for y,x in peaks[:10]])

def main():
 if (OUT/'direct_results.json').exists():raise FileExistsError('Preserve')
 labels=json.loads(SHAFT.read_text());axes={}
 for run in [1,3,4]:
  z=[x['direction'] for x in labels if x.get('run')==run and x.get('direction')]
  m=np.mean(z,axis=0);m/=np.linalg.norm(m);axes[run]=[float(m[0]),float(m[1]),0.0]
 # image-plane directions are observed shaft-line means; signs are line ambiguity, scalar search keeps both signs.
 axes[2]=[0.0,0.0,-1.0]
 save('direct_axis_definitions.json',dict(axes=axes,run2=dict(status='assumed_optical_axis',axis=[0,0,-1],reason='user confirms shaft approximately normal to image plane; not metric calibration'),frame='event image x right, y down, z=away from camera',source='independent shaft projection for run1/3/4; explicit physical z-axis assumption for run2; no free CMax axis'))
 hashes={str(p.relative_to(ROOT)):dig(p) for p in [OUT/'direct_axis.py',OUT/'direct_axis_definitions.json',PHASE/'results.jsonl',SHAFT,ROOT/'g1_spin/src/g1_spin/spherical_cmax.py']};save('direct_freeze.json',hashes)
 rows=[json.loads(x) for x in (PHASE/'results.jsonl').read_text().splitlines()];out=[]
 for run in [1,2,3,4]:
  for r in [x for x in rows if x['run']==run and x['method']=='G' and x['measurement_valid']]:
   d=np.load(PHASE/'windows'/f"{r['window']}.npz");points=cmax.lift(d['xy'][:,0],d['xy'][:,1],d['center'][0],d['center'][1],float(d['radius_px']));t=(d['time_us']-d['time_us'][0])*1e-6
   x=scalar(points,t,axes[run]);out.append(dict(run=run,window=r['window'],direct_axis=x,free_rpm=r['measurement_rpm'],free_axis=r['measurement_axis']))
 save('direct_results.json',out);assert all(dig(ROOT/p)==h for p,h in hashes.items());save('direct_verification.json',dict(rows=len(out),sources_unchanged=True));print('direct',len(out))
if __name__=='__main__':main()
