"""Standalone old geometry for the exposure unavailable to the new frontend."""
import importlib.util,json
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('front',OUT/'run.py');f=importlib.util.module_from_spec(s);s.loader.exec_module(f)
if (OUT/'old_standalone.json').exists():raise FileExistsError('Preserve output')
p2=f.base.p2;old=f.old;run=3;end=4335000;begin=end-20000;folder=p2.ARCHIVE/'decoded'/p2.NAMES[run]
m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/x,mmap_mode='r') for x in m['parts']];origin=int(parts[0]['t'][0]);ts,xy=p2.observations(run);xy=old.smooth(ts,xy)
p2.centers(ts,xy,[begin,end]);raw,ids=old.read_window(parts,origin,begin,end);times=raw['t'].astype(np.int64)-origin
centers=p2.centers(ts,xy,times);pixels=np.column_stack([raw['x'],raw['y']]).astype(float);ix=np.flatnonzero(np.linalg.norm(pixels-centers,axis=1)<=.8*59)
if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
rel=pixels[ix]-centers[ix];cmax=f.base.prev.cmax;result=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,59),(times[ix]-begin)*1e-6,max_rpm=3500)
np.savez_compressed(OUT/'spin_inputs/old_standalone.npz',ids=ids[ix],events=raw[ix],centers=centers[ix],origin_us=origin,radius=59)
old.save(OUT/'old_standalone.json',dict(case='run3_04335000',method='old_median3_standalone',result=result,reason_for_supplement='paired preparation rejected both when only new geometry lacked coverage'))
print(result['valid'],result['spin_rpm'])
