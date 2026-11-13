"""Independent boolean-window/explicit NCC check of registered recurrence peaks."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import convolve
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
records=json.loads((OUT/'recurrence.json').read_text())+json.loads((OUT/'high_motion_recurrence/recurrence.json').read_text());checks=[];sources={}
k=np.exp(-np.arange(-2,3)**2/2);k=k/k.sum();kernel=(k[:,None]*k[None,:])[:,:,None]
for name in ['222403','221756']:
 folder=ROOT/f'g1_spin/data/decoded/local_analysis_{name}';manifest=json.loads((folder/'manifest.json').read_text());paths=[folder/'manifest.json']+[folder/f for f in manifest['parts']];e=np.concatenate([np.load(p) for p in paths[1:]]);origin=int(e['t'][0])
 for p in paths:sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
 def image(t):
  mid=round(t*1e6);q=e[(e['t']>=origin+mid-5000)&(e['t']<origin+mid+5000)];im=np.zeros((720,1280,2),np.float32);im[q['y'],q['x'],q['p']]=1;return convolve(im,kernel,mode='mirror')
 for r in [r for r in records if r['recording']==name and not r['background'] and r['highest_peak']]:
  i=r['highest_peak']['index'];src=r['sources'][2];sx,sy=src['template_xy'];dx,dy=src['shifts'][i]
  a=image(src['source_s'])[sy:sy+40,sx:sx+40].astype(float);b=image(src['source_s']+r['lags_s'][i])[sy+dy:sy+dy+40,sx+dx:sx+dx+40].astype(float)
  explicit=float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b)));stored=src['scores'][i];assert abs(explicit-stored)<2e-6,(explicit,stored)
  checks.append(dict(recording=name,anchor_s=r['anchor_s'],source_s=src['source_s'],lag=r['lags_s'][i],explicit_similarity=explicit,opencv_similarity=stored,error=abs(explicit-stored)))
 del e
(OUT/'recurrence_verification.json').write_text(json.dumps(dict(checks=checks,sources_sha256=sources),indent=2))
print('Verified',len(checks),'peak scores with independent window selection, convolution and NCC')
