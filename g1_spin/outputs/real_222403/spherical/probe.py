import sys,json,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'src'))
from g1_spin.events import EventBatch
from g1_spin.spherical_cmax import lift,estimate
root=Path(__file__).resolve().parents[3]/'data/decoded/local_analysis_222403'
m=json.loads((root/'manifest.json').read_text()); e=np.concatenate([np.load(root/f) for f in m['parts']]); t0=int(e['t'][0])
rows=[]
for sec,cx in [(1.0,547),(7.,556),(12.,552)]:
 for duration in [.02,.05,.1]:
  for polarity in [None]:
   a,b=np.searchsorted(e['t'],[t0+round(sec*1e6),t0+round((sec+duration)*1e6)]); s=e[a:b]; cy=224;r=49
   mask=((s['x'].astype(float)-cx)**2+(s['y'].astype(float)-cy)**2 < (r*.8)**2)&(s['y']<240)
   s=s[mask];dt=(s['t'].astype(np.int64)-(t0+round(sec*1e6)))*1e-6
   tick=time.perf_counter(); result=estimate(lift(s['x'],s['y'],cx,cy,r),dt,max_rpm=600)
   result.update(relative_start_s=sec,window_s=duration,roi=[cx,cy,r],runtime_s=time.perf_counter()-tick)
   rows.append(result);print(sec,duration,len(s),result['valid'],result['spin_rpm'],result['axis_camera'],result['failure_reason'],flush=True)
Path(__file__).with_name('pilot.json').write_text(json.dumps(rows,indent=2))
