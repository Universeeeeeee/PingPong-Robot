import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
D=ROOT/'g1_spin/data/archive_20261105/decoded/DSP03F100025_20261105_201143';m=json.loads((D/'manifest.json').read_text());e=np.concatenate([np.load(D/f) for f in m['parts']]);t=(e['t'].astype(np.int64)-int(e['t'][0]))*1e-6
# Fixed ROI based on previous geometry prior, preserves movement; not tuned to period.
mask=(e['x']>=500)&(e['x']<610)&(e['y']>=365)&(e['y']<475);q=e[mask];tt=t[mask];nb=int(np.ceil(t[-1]/.002))+1
frame=np.floor(tt/.002).astype(int);xx=((q['x'].astype(int)-500)//5);yy=((q['y'].astype(int)-365)//5);feature=q['p'].astype(int)*22*22+yy*22+xx
v=np.bincount(frame*968+feature,minlength=nb*968).reshape(nb,968).astype(float);v/=np.maximum(np.linalg.norm(v,axis=1,keepdims=True),1e-12)
results=[]
for lo,hi in [(1,3),(3,5),(5,7)]:
 a,b=int(lo/.002),min(int(hi/.002),len(v));seq=v[a:b];lags=np.arange(25,min(351,len(seq)//2));scores=[float(np.mean(np.sum(seq[:-k]*seq[k:],axis=1))) for k in lags];indices,_=find_peaks(scores,prominence=.015,distance=12)
 peaks=[dict(lag_s=float(lags[i]*.002),rpm_if_one_repeat=float(60/(lags[i]*.002)),score=scores[i]) for i in indices];results.append(dict(start=lo,end=hi,lags_s=(lags*.002).tolist(),scores=scores,peaks=peaks));print(lo,hi,peaks)
(OUT/'periodicity.json').write_text(json.dumps(dict(roi=[500,365,610,475],frame_ms=2,spatial_cell_px=5,results=results,warning='same sensor, fixed ROI, no unique one-turn guarantee'),indent=2))
