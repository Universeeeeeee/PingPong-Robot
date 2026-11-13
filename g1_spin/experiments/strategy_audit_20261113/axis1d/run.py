"""Known-axis scalar CMAX on fixed-axis drill windows; no RPM labels in fitting."""
import json,hashlib,importlib.util
from pathlib import Path
from scipy.optimize import minimize_scalar
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3];PHASE=ROOT/'g1_spin/experiments/drill_phase3_20261106'
s=importlib.util.spec_from_file_location('cmax',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py');cmax=importlib.util.module_from_spec(s);s.loader.exec_module(cmax)

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def fit_axis(points,t,axis,max_rpm):
 axis=np.asarray(axis,float);axis/=np.linalg.norm(axis);bound=max_rpm*2*np.pi/60
 def score(s):return float(np.var(cmax.iwe(cmax.warp(points,t,axis*s))/len(points)))
 grid=np.linspace(-bound,bound,281);values=np.array([score(s) for s in grid]);peaks=[]
 for i in range(1,len(grid)-1):
  if values[i]>=values[i-1] and values[i]>=values[i+1]:
   local=minimize_scalar(lambda x:-score(x),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-5})
   peaks.append((float(-local.fun),float(local.x)))
 peaks.sort(reverse=True);best=peaks[0] if peaks else (float(values.max()),float(grid[values.argmax()]))
 return dict(axis_camera=axis.tolist(),valid=bool(abs(best[1])>1e-12),omega_rad_s=(axis*best[1]).tolist(),spin_rpm=float(abs(best[1])*60/(2*np.pi)),signed_rpm=float(best[1]*60/(2*np.pi)),contrast=float(best[0]),zero_contrast=score(0.),relative_gain=float(best[0]/max(score(0.),1e-30)-1),peaks=[dict(contrast=c,signed_rpm=v*60/(2*np.pi)) for c,v in peaks[:8]],max_rpm=max_rpm,method='scalar_known_axis_cmax')

def main():
 if (OUT/'results.json').exists():raise FileExistsError('Preserve output')
 rows=[json.loads(x) for x in (PHASE/'results.jsonl').read_text().splitlines()];g=[r for r in rows if r['method']=='G' and r['measurement_valid']]
 # Use recordings with physical fixed-axis/pose labels; run1 and fixed run4 are vertical/reverse, run2 faces camera, run3 tilted.
 dirs=json.loads((PHASE/'direction_checks.json').read_text());run1=np.array(dirs['run1_mean_directed_axis']);run1/=np.linalg.norm(run1)
 axis_sets={'run1_vertical_image':[0,-1,0],'run1_free_mean':run1,'run2_data_mean':None,'run3_data_mean':None,'run4_vertical_image':[0,1,0],'run4_free_reverse':-run1}
 for run,key in [(2,'run2_data_mean'),(3,'run3_data_mean')]:
  a=np.array([r['measurement_axis'] for r in g if r['run']==run]);axis_sets[key]=np.mean(a,axis=0);axis_sets[key]/=np.linalg.norm(axis_sets[key])
 selected=[r for r in g if r['run'] in [1,2,3,4]]
 hashes={str(p.relative_to(ROOT)):digest(p) for p in [OUT/'run.py',OUT/'protocol.md',PHASE/'results.jsonl',PHASE/'direction_checks.json',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py']};save('freeze.json',hashes);save('axis_definitions.json',{k:np.asarray(v).tolist() for k,v in axis_sets.items()})
 out=[]
 for r in selected:
  path=PHASE/'windows'/f"{r['window']}.npz";d=np.load(path);points=cmax.lift(d['xy'][:,0],d['xy'][:,1],d['center'][0],d['center'][1],float(d['radius_px']));t=(d['time_us']-d['time_us'][0])*1e-6
  applicable=['run1_vertical_image','run1_free_mean'] if r['run']==1 else ['run4_vertical_image','run4_free_reverse'] if r['run']==4 else ['run2_data_mean'] if r['run']==2 else ['run3_data_mean']
  for key in applicable:
   result=fit_axis(points,t,axis_sets[key],3500)
   out.append(dict(window=r['window'],run=r['run'],condition=r['condition'],axis_source=key,free_rpm=r['measurement_rpm'],free_axis=r['measurement_axis'],**result))
   print(r['window'],key,round(result['signed_rpm'],1),flush=True)
 save('results.json',out);assert all(digest(ROOT/p)==h for p,h in hashes.items());save('verification.json',dict(rows=len(out),source_hashes_unchanged=True,axis_sources='run1/4 physical image-axis hypotheses; run2/3 data-summary directions for diagnostic only'))

if __name__=='__main__':main()
