from pathlib import Path
import numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('g1_spin/data/decoded/local_analysis_222403');out=Path('g1_spin/outputs/real_222403/periodicity');meta=json.load(open(root/'manifest.json'));first=meta['statistics']['first_timestamp_us']
times=np.array([.306,.446,.586,6.08,6.296,6.512,11.006,11.416,11.826]);imgs=np.zeros((9,90,120))
for part in meta['parts']:
 a=np.load(root/part);t=(a['t'].astype(np.int64)-first)/1e6
 for i,tm in enumerate(times):
  m=(abs(t-tm)<.01)&(a['x']>=500)&(a['x']<620)&(a['y']>=160)&(a['y']<250)
  np.add.at(imgs[i],(a['y'][m].astype(int)-160,a['x'][m].astype(int)-500),1)
fig,axs=plt.subplots(3,3,figsize=(10,8));titles=['Rate peak','Half-period later','Next rate peak']
for i,ax in enumerate(axs.flat):
 ax.imshow(np.log1p(imgs[i]),cmap='magma',vmin=0,vmax=2.8);ax.set_title(f'{titles[i%3]} / {times[i]:.3f} s');ax.axis('off')
fig.suptitle('20 ms event crops: recurrent patterns vs half-cycle (not intensity frames)');fig.tight_layout();fig.savefig(out/'repeated_patterns.png',dpi=180)
