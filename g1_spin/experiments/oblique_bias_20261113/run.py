import importlib.util,json,hashlib
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];PHASE=ROOT/'g1_spin/experiments/drill_phase3_20261106'
s=importlib.util.spec_from_file_location('axis',ROOT/'g1_spin/experiments/strategy_audit_20261113/axis1d/direct_axis.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
axes=json.loads((ROOT/'g1_spin/experiments/strategy_audit_20261113/axis1d/direct_axis_definitions.json').read_text());axis=np.array(axes['axes']['3'])
shaft={x['window']:x for x in json.loads((PHASE/'shaft_projection.json').read_text())}
rows=[json.loads(x) for x in (PHASE/'results.jsonl').read_text().splitlines()];rows=[x for x in rows if x['run']==3 and x['method']=='G' and x['measurement_valid']]
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 if (OUT/'results.json').exists():raise FileExistsError('Preserve')
 paths=[OUT/'run.py',OUT/'protocol.md',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py',ROOT/'g1_spin/experiments/strategy_audit_20261113/axis1d/direct_axis.py']+[PHASE/'windows'/f"{r['window']}.npz" for r in rows]
 save('freeze.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});out=[]
 for row in rows:
  d=np.load(PHASE/'windows'/f"{row['window']}.npz");xy=d['xy'];center=d['center'];r=float(d['radius_px']);t=(d['time_us']-d['time_us'][0])*1e-6;tm=t-(t[0]+t[-1])/2
  pts=a.cmax.lift(xy[:,0],xy[:,1],*center,r);local=np.r_[shaft[row['window']]['direction'],0]
  jobs=[('reference_first',pts,t,axis),('reference_middle',pts,tm,axis),('per_window_shaft',pts,tm,local)]
  if row['window'] in ['run3_01495000','run3_04365000']:
   for dr in [-2,2,4]:jobs.append((f'radius{dr:+}',a.cmax.lift(xy[:,0],xy[:,1],*center,r+dr),tm,axis))
   for dim in [0,1]:
    for delta in [-2,2]:
     c=center.copy();c[dim]+=delta;jobs.append((f'center{dim}_{delta:+}',a.cmax.lift(xy[:,0],xy[:,1],*c,r),tm,axis))
   for deg in [-15,-5,5,15]:
    theta=np.arctan2(axis[1],axis[0])+np.deg2rad(deg);jobs.append((f'azimuth{deg:+}',pts,tm,[np.cos(theta),np.sin(theta),0]))
   for deg in [-15,15]:
    th=np.deg2rad(deg);jobs.append((f'tilt{deg:+}',pts,tm,[axis[0]*np.cos(th),axis[1]*np.cos(th),np.sin(th)]))
   for polarity in [0,1]:
    mask=d['polarity']==polarity;jobs.append((f'polarity{polarity}',pts[mask],tm[mask],axis))
  for name,points,tt,ax in jobs:
   result=a.scalar(points,tt,ax);free=None
   if name in ['reference_middle','radius+4','center0_+2','center1_+2']:
    free=a.cmax.estimate(points,tt,max_rpm=3500)
   out.append(dict(window=row['window'],variant=name,events=len(points),one=result,free=free,old_free_rpm=row['measurement_rpm']))
  print(row['window'],flush=True)
 save('results.json',out)
 for filename,h in json.loads((OUT/'freeze.json').read_text()).items():assert hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==h
 print('rows',len(out))
if __name__=='__main__':main()
