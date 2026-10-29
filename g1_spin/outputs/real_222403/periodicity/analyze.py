import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks, periodogram
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('g1_spin/data/decoded/local_analysis_222403'); out=Path('g1_spin/outputs/real_222403/periodicity')
meta=json.loads((root/'manifest.json').read_text()); first=meta['statistics']['first_timestamp_us']; dur=meta['statistics']['duration_us']/1e6
bin_us=2000; n=int(np.ceil(dur*1e6/bin_us))+1
# 12x9 bins inside crop, two polarity channels
features=np.zeros((n,9*12*2)); bg=np.zeros(n); roi=np.zeros(n)
for name in meta['parts']:
 a=np.load(root/name); t=((a['t'].astype(np.int64)-first)//bin_us).astype(int)
 mask=(a['x']>=500)&(a['x']<620)&(a['y']>=160)&(a['y']<250)
 b=a[mask]; tb=t[mask]; k=((b['y'].astype(int)-160)//10*12+(b['x'].astype(int)-500)//10)*2+b['p']
 np.add.at(features,(tb,k),1)
 np.add.at(roi,tb,1)
 # equal size background patch immediately left, excluding object
 back=(a['x']>=350)&(a['x']<470)&(a['y']>=160)&(a['y']<250)
 np.add.at(bg,t[back],1)
np.savez_compressed(out/'features.npz',features=features,background=bg,roi=roi,bin_us=bin_us,first_us=first)
# Normalize away global rate per time step, preserve relative spatial features. Median subtract each pixel channel.
features=gaussian_filter1d(features,.7,axis=0)
features=features/np.maximum(features.sum(axis=1,keepdims=True),1)
rows=[]; curves=[]
for center in np.arange(1.0,16.6,.5):
 start=max(0,center-1);end=min(dur,center+1)
 i,j=int(start*500),int(end*500)
 x=features[i:j]; x=x-x.mean(axis=0)
 freqs=np.fft.rfftfreq(len(x),.002)
 spec=np.sum(abs(np.fft.rfft(x*np.hanning(len(x))[:,None],axis=0))**2,axis=1)
 sf=np.fft.rfft(x,n=2*len(x),axis=0)
 ac=np.fft.irfft(np.sum(sf*np.conj(sf),axis=1),n=2*len(x))[:len(x)]
 ac=ac/np.maximum(np.arange(len(x),0,-1),1)
 ac/=max(ac[0],1e-20)
 lo,hi=10,min(len(ac)//2,500)
 peaks,_=find_peaks(ac[lo:hi],prominence=.02)
 peaks=peaks+lo
 peaks=sorted(peaks,key=lambda k:ac[k],reverse=True)
 fmask=(freqs>=1)&(freqs<=40); fidx=np.flatnonzero(fmask)
 pk,_=find_peaks(spec); pk=[int(k) for k in pk if fmask[k]]; pk=sorted(pk,key=lambda k:spec[k],reverse=True)[:6]
 # Report top recurrence peaks, no claim strongest is full revolution.
 rec=[dict(lag_s=float(k*.002),frequency_hz=float(1/(k*.002)),correlation=float(ac[k])) for k in peaks[:6]]
 spectral=[dict(frequency_hz=float(freqs[k]),power=float(spec[k])) for k in pk]
 rows.append(dict(center_s=float(center),start_s=start,end_s=end,recurrence=rec,spectral=spectral,mean_roi_rate=float(roi[i:j].mean()*500),mean_bg_rate=float(bg[i:j].mean()*500)))
 curves.append(ac[:500])
json.dump(rows,open(out/'recurrence.json','w'),indent=2)
fig,axs=plt.subplots(3,1,figsize=(12,10))
times=np.arange(n)*.002
axs[0].plot(times,gaussian_filter1d(roi,25)*500,label='Ball crop');axs[0].plot(times,gaussian_filter1d(bg,25)*500,label='Background patch');axs[0].legend();axs[0].set_ylabel('Events/s')
axs[1].imshow(np.array(curves).T,origin='lower',extent=[1,16.5,0,1],aspect='auto',vmin=0,vmax=.8,cmap='magma');axs[1].set_ylabel('Pattern recurrence lag (s)')
for row in rows:
 for idx,r in enumerate(row['recurrence'][:3]):axs[2].scatter(row['center_s'],r['frequency_hz'],s=30/(idx+1),alpha=.5,c='C0')
axs[2].set_ylim(0,40);axs[2].set_ylabel('Candidate repetition Hz');axs[2].set_xlabel('Seconds since first event')
fig.tight_layout();fig.savefig(out/'recurrence.png',dpi=150)
for row in rows[::4]: print(row['center_s'],[(round(r['frequency_hz'],2),round(r['correlation'],2)) for r in row['recurrence'][:4]],'spec',[(s['frequency_hz']) for s in row['spectral'][:3]])
