"""Recover standalone baseline only where shared 20ms preparation was too strict."""
import importlib.util
import json
from pathlib import Path
from time import perf_counter
import numpy as np
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('trial',OUT/'run.py');t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t)
path=OUT/'baseline_supplement.json'
if path.exists():raise FileExistsError('Preserve supplement')
records=[]
for item in json.loads((OUT/'baseline_eligibility.json').read_text()):
    if not item['baseline_10ms_geometry_available']:continue
    name=item['window'];run=int(name[3]);end=int(name.split('_')[1]);begin=end-10000
    folder=t.p2.ARCHIVE/'decoded'/t.p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text())
    parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']];origin=int(parts[0]['t'][0])
    raw,ids=t.old.read_window(parts,origin,begin,end);ts,xy=t.p2.observations(run)
    center=t.p2.centers(ts,xy,[end-5000])[0];pixels=np.column_stack([raw['x'],raw['y']]).astype(float);radius=t.RADIUS[run]
    selected=np.flatnonzero(np.linalg.norm(pixels-center,axis=1)<=.8*radius)
    if len(selected)>5000:selected=selected[np.linspace(0,len(selected)-1,5000).astype(int)]
    e=raw[selected];rel=pixels[selected]-center
    start=perf_counter();result=t.cmax.estimate(t.cmax.lift(rel[:,0],rel[:,1],0,0,radius),(e['t'].astype(np.int64)-origin-begin)*1e-6,max_rpm=3500)
    records.append(dict(window=name,run=run,begin_us=begin,end_us=end,method='standalone_baseline',result=result,elapsed_s=perf_counter()-start))
    np.savez_compressed(OUT/'inputs'/f'{name}_standalone_baseline.npz',ids=ids[selected],events=e,centers=np.tile(center,(len(e),1)),radius=radius,origin_us=origin)
    del parts
path.write_text(json.dumps(records,indent=2,allow_nan=False))
print(json.dumps([dict(window=r['window'],valid=r['result']['valid'],rpm=r['result']['spin_rpm']) for r in records]))
