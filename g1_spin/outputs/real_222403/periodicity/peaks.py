from pathlib import Path
import numpy as np,json,csv
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks,periodogram
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('g1_spin/outputs/real_222403/periodicity');data=np.load(out/'features.npz');dt=.002
rate=gaussian_filter1d(data['roi'],25)/dt; bg=gaussian_filter1d(data['background'],25)/dt;t=np.arange(len(rate))*dt
pk,props=find_peaks(rate,prominence=16000,distance=90)
period=np.diff(pk)*dt; mids=(pk[:-1]+pk[1:])*.5*dt
rows=[dict(time_s=float(tm),previous_peak_s=float(a*dt),next_peak_s=float(b*dt),period_s=float(p),candidate_hz=float(1/p),candidate_rpm=float(60/p)) for tm,a,b,p in zip(mids,pk[:-1],pk[1:],period)]
with (out/'cycle_candidates.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
# Independent low-frequency rate power ball vs equally sized background; low-frequency candidate power not simply illumination.
comparisons=[]
for start,end in [(0,4),(4,8),(8,12),(12,16)]:
 a,b=int(start/dt),int(end/dt)
 fs,s=periodogram(data['roi'][a:b],fs=1/dt);_,bs=periodogram(data['background'][a:b],fs=1/dt)
 band=(fs>=1)&(fs<=10)
 k=np.flatnonzero(band)[np.argmax(s[band])]
 corr=float(np.corrcoef(rate[a:b],bg[a:b])[0,1])
 comparisons.append(dict(start_s=start,end_s=end,strongest_roi_rate_hz=float(fs[k]),roi_spectral_power=float(s[k]),background_spectral_power=float(bs[k]),smoothed_rate_correlation=corr))
json.dump(dict(engineering_parameters={'smooth_sigma_ms':50,'peak_prominence_events_per_s':16000,'min_peak_distance_ms':180},peak_times_s=(pk*dt).tolist(),background_comparisons=comparisons,caveat='Peak repetition may be one revolution or harmonics. No independent RPM truth; no 3D axis.'),open(out/'cycle_method.json','w'),indent=2)
fig,ax=plt.subplots(2,1,figsize=(12,7));ax[0].plot(t,rate,label='Ball crop, 50ms Gaussian');ax[0].scatter(pk*dt,rate[pk],s=20,color='C1',label='Repeat candidates');ax[0].plot(t,bg,label='Equal-area background');ax[0].set_ylabel('Events/s');ax[0].legend()
ax[1].plot(mids,60/period,'o-');ax[1].set_ylabel('Conditional RPM (one repeat/rev)');ax[1].set_xlabel('Seconds from first event');ax[1].set_ylim(0,250);ax[1].grid(alpha=.2)
fig.tight_layout();fig.savefig(out/'candidate_slowdown.png',dpi=150)
print('peaks',pk*dt)
print('candidates',[(round(r['time_s'],2),round(r['candidate_rpm'],1)) for r in rows])
print(comparisons)
