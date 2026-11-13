"""Predeclared supplementary high-center-velocity interval from prior geometry audit."""
import json,importlib.util
from pathlib import Path
from time import perf_counter
import numpy as np
from run import ROOT,OUT,DATA,MODELS,cmax,center,save

if (OUT/'high_motion_results.json').exists():raise FileExistsError('Preserve supplement')
(OUT/'high_motion_protocol.txt').write_text('补充0.30–0.60s历史几何高平移段：已有中心模型约46px/s。中点0.35/0.45/0.55s×20/50/100ms×三中心=27，中央50ms三法3500RPM=3，共30；独立图案复现source0.45s。固定参数不按本轮结果重拟几何。此段旧名translation只是几何诊断标签，不代表球无自转。\n')
m=json.loads((DATA/'manifest.json').read_text());e=np.concatenate([np.load(DATA/p) for p in m['parts']]);origin=int(e['t'][0]);t=e['t'].astype(np.int64)-origin;models=json.loads(MODELS.read_text())['translation'];rows=[]
for mid in [350000,450000,550000]:
 for ms in [20,50,100]:
  lo,hi=mid-ms*500,mid+ms*500;a,b=np.searchsorted(t,[lo,hi]);raw=e[a:b];tt=t[a:b];pix=np.column_stack([raw['x'],raw['y']]).astype(float)
  cs={'fixed50':np.tile(center(models['50'],mid),(len(raw),1)),'linear50':center(models['50'],tt),'linear100':center(models['100'],tt)}
  mask=(pix[:,0]>=475)&(pix[:,0]<638)&(pix[:,1]>=143)&(pix[:,1]<240)
  for c in cs.values():mask&=np.linalg.norm(pix-c,axis=1)<=.8*48
  ix=np.flatnonzero(mask);common=ix+a
  if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
  name=f'high_motion_{mid:08d}_{ms}ms';np.savez_compressed(OUT/'inputs'/f'{name}.npz',events=raw[ix],ids=ix+a,common_ids=common,origin_us=origin,**{k:v[ix] for k,v in cs.items()})
  jobs=[(method,600) for method in cs]
  if mid==450000 and ms==50:jobs += [(method,3500) for method in cs]
  for method,bound in jobs:
   rel=pix[ix]-cs[method][ix];start=perf_counter();result=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,48),(tt[ix]-lo)*1e-6,max_rpm=bound)
   rows.append(dict(window=name,segment='translation',mid_us=mid,window_ms=ms,method=method,radius=48.,max_rpm=bound,raw_count=len(raw),common_count=len(common),input_count=len(ix),elapsed_s=perf_counter()-start,result=result))
  print(name,flush=True)
save('high_motion_results.json',rows)
s=importlib.util.spec_from_file_location('recurrence',OUT/'recurrence.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r);r.OUT=OUT/'high_motion_recurrence';r.OUT.mkdir(exist_ok=True);r.JOBS={'222403':[.45]};r.main()
