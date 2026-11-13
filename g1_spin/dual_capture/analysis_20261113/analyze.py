import json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from audit_remote import OUT,save

def main():
 d=np.load(OUT/'color_track.npz');tr=d['track'];crops=d['crops'];e=json.loads((OUT/'event_track.json').read_text());e=[dict(t=x['t'],candidates=[c for c in x['candidates'] if 310<c[1]<390 and 60<c[2]<82 and 300<c[0]<700]) for x in e];e=[x for x in e if len(x['candidates'])==1];et=np.array([x['t'] for x in e]);exy=np.array([x['candidates'][0][:2] for x in e]);scans=[]
 for lag in np.arange(-3,3.001,.01):
  y=np.column_stack([np.interp(et+lag,tr[:,0],tr[:,k]) for k in [1,2]]);X=np.c_[exy,np.ones(len(exy))];train=(et<28);test=~train
  init=np.linalg.lstsq(X[train],y[train],rcond=None)[0];fit=least_squares(lambda z:(X[train]@z.reshape(3,2)-y[train]).ravel(),init.ravel(),loss='soft_l1',f_scale=2)
  res=np.linalg.norm(X@fit.x.reshape(3,2)-y,axis=1);scans.append(dict(lag=float(lag),train=float(np.median(res[train])),test=float(np.median(res[test])),affine=fit.x.reshape(3,2).tolist()))
 best=min(scans,key=lambda x:x['train']);save('alignment_preliminary.json',dict(best=best,candidates=sorted(scans,key=lambda x:x['train'])[:10],train_n=int(train.sum()),test_n=int(test.sum()),warning='Uncalibrated spatial affine diagnostic, no clock drift fit; not final synchronization'))
 yy,xx=np.mgrid[:64,:64];mask=(xx-32)**2+(yy-32)**2<22**2
 # Interior RGB patterns, per-pixel temporal mean removed to suppress static lighting.
 features=crops[:,mask,:].astype(float).reshape(len(crops),-1);out=[]
 for begin,end in [(15,21),(16,19),(29,32),(36,39),(40,45),(41,44)]:
  good=(tr[:,0]>=begin)&(tr[:,0]<end);tt=tr[good,0];f=features[good];f=f-f.mean(axis=0);f/=np.maximum(np.linalg.norm(f,axis=1,keepdims=True),1e-9);lags=[];scores=[]
  for lag in range(3,min(100,len(f)//2)):
   valid=np.abs((tt[lag:]-tt[:-lag])-lag/89.9)<.008
   if valid.sum()<30:continue
   lags.append(float(np.median(tt[lag:][valid]-tt[:-lag][valid])));scores.append(float(np.mean(np.sum(f[lag:][valid]*f[:-lag][valid],axis=1))))
  peaks,_=find_peaks(scores,prominence=.03);tops=sorted(peaks,key=lambda i:scores[i],reverse=True)[:5]
  out.append(dict(begin=begin,end=end,frames=len(tt),radius_median=float(np.median(tr[good,3])),peaks=[dict(period_s=lags[i],rpm=60/lags[i],correlation=scores[i]) for i in tops],lags=lags,scores=scores))
 save('rgb_periodicity.json',out);print(json.dumps(dict(alignment=best,rgb=[{k:v for k,v in x.items() if k not in ['lags','scores']} for x in out]),indent=2),flush=True)
if __name__=='__main__':main()
