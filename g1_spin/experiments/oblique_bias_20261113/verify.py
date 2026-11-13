import json,hashlib,importlib.util
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.ndimage import gaussian_filter
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];PHASE=ROOT/'g1_spin/experiments/drill_phase3_20261106'
rows=json.loads((OUT/'results.json').read_text());assert len(rows)==54
checks=[]
for r in rows:
 d=np.load(PHASE/'windows'/f"{r['window']}.npz");xy=d['xy'];center=d['center'].copy();radius=float(d['radius_px']);v=r['variant'];time=(d['time_us']-d['time_us'][0])*1e-6
 if v!='reference_first':time-=(time[0]+time[-1])/2
 if v.startswith('radius'):radius+=int(v[6:])
 if v.startswith('center'):center[int(v[6])]+=int(v.split('_')[1])
 if v.startswith('polarity'):
  m=d['polarity']==int(v[-1]);xy=xy[m];time=time[m]
 p=(xy-center)/radius;p=np.column_stack([p,-np.sqrt(1-np.sum(p*p,axis=1))]);w=np.array(r['one']['omega_rad_s']);moved=Rotation.from_rotvec(-time[:,None]*w).apply(p)
 yy=(np.arcsin(np.clip(moved[:,2],-1,1))/np.pi+.5)*63;xx=((np.arctan2(moved[:,1],moved[:,0])/(2*np.pi)+.5)*128)%128
 iy=np.floor(yy).astype(int);ix=np.floor(xx).astype(int);fy=yy-iy;fx=xx-ix;im=np.zeros((64,128))
 for dy,wy in [(0,1-fy),(1,fy)]:
  for dx,wx in [(0,1-fx),(1,fx)]:np.add.at(im,(np.minimum(iy+dy,63),(ix+dx)%128),wy*wx)
 im=gaussian_filter(im,1,truncate=2,mode=('reflect','wrap'));score=float(np.var(im/len(p)));assert np.isclose(score,r['one']['contrast'],rtol=1e-8,atol=1e-18)
 assert np.isclose(np.linalg.norm(w)*60/(2*np.pi),r['one']['spin_rpm']);checks.append(dict(window=r['window'],variant=v,score_abs_error=abs(score-r['one']['contrast'])))
for name,h in json.loads((OUT/'freeze.json').read_text()).items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
(OUT/'verification.json').write_text(json.dumps(dict(rows=54,independent_scores=checks,frozen_inputs_unchanged=True),indent=2));print('54 scores independently reconstructed')
