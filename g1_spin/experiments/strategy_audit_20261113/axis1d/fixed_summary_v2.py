"""Correct run4 fixed-interval aggregation; preserve prior outputs."""
import json,importlib.util
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('direct',P/'direct_axis.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
output=P/'fixed_summary_v2.json'
if output.exists():raise FileExistsError('Preserve corrected results')
previous=json.loads((P/'direct_results.json').read_text())
free=[json.loads(l) for l in (d.PHASE/'results.jsonl').read_text().splitlines()];free={r['window']:r for r in free if r['method']=='G'}
labels=json.loads(d.SHAFT.read_text())
shaft=[r for r in labels if r['run']==4 and r['condition']=='reverse_fixed' and r.get('direction') is not None]
a=np.mean([r['direction'] for r in shaft],axis=0);a/=np.linalg.norm(a);axis=np.r_[a,0.]
corrected=[]
for r in previous:
 if r['run']!=4:corrected.append(r);continue
 if free[r['window']]['condition']!='reverse_fixed':continue
 data=np.load(d.PHASE/'windows'/f"{r['window']}.npz")
 points=d.cmax.lift(data['xy'][:,0],data['xy'][:,1],data['center'][0],data['center'][1],float(data['radius_px']))
 times=(data['time_us']-data['time_us'][0])*1e-6
 corrected.append(dict(**{k:v for k,v in r.items() if k!='direct_axis'},direct_axis=d.scalar(points,times,axis)))
summary=[]
for run in [1,2,3,4]:
 rr=[r for r in corrected if r['run']==run];signed=np.array([r['direct_axis']['signed_rpm'] for r in rr]);rpm=np.array([r['direct_axis']['spin_rpm'] for r in rr]);a=np.array(rr[0]['direct_axis']['axis_camera']);b=np.mean([r['free_axis'] for r in rr],axis=0);b/=np.linalg.norm(b)
 assert np.allclose(abs(signed),rpm,atol=1e-10,rtol=0)
 summary.append(dict(run=run,count=len(rr),reference_axis=a.tolist(),free_mean_directed_axis=b.tolist(),unsigned_axis_angle_deg=float(np.degrees(np.arccos(np.clip(abs(a@b),0,1)))),signed_rpm_median=float(np.median(signed)),rpm_median=float(np.median(rpm)),free_rpm_median=float(np.median([r['free_rpm'] for r in rr])),negative_windows=int((signed<0).sum()),positive_windows=int((signed>0).sum())))
result=dict(note='run4 uses reverse_fixed only for both shaft definition and paired RPM/axis aggregation. Other runs reuse existing results. Selected on prior free-CMAX validity; not all-window acceptance rate. Scalar valid means nonzero candidate, not production quality gate.',run4_shaft_windows=[r['window'] for r in shaft],rows=corrected,summary=summary)
output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,indent=2))
