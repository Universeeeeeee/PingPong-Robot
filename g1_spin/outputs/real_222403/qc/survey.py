import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/Users/vae/Projects/PingPongRobot/g1_spin/data/decoded/local_analysis_222403')
out=Path(__file__).parent
m=json.loads((root/'manifest.json').read_text())
e=np.concatenate([np.load(root/f) for f in m['parts']])
t0=int(e['t'][0]); ts=e['t'].astype(np.int64)-t0
h,w=m['height'],m['width']
counts=np.bincount(e['y'].astype(int)*w+e['x'],minlength=h*w).reshape(h,w)
plt.figure(figsize=(16,9));plt.imshow(np.log1p(counts),cmap='magma');plt.colorbar(label='log(1+events)');plt.title('Entire 17.58 s accumulation (sensor pixels)');plt.savefig(out/'full_accumulation.png',dpi=150);plt.close()
times=np.arange(.25,17.5,1)
fig,axs=plt.subplots(5,4,figsize=(20,15))
for ax,t in zip(axs.flat,times):
 a,b=np.searchsorted(ts,[int(t*1e6),int((t+.1)*1e6)])
 s=e[a:b];c=np.bincount(s['y'].astype(int)*w+s['x'],minlength=h*w).reshape(h,w)
 ax.imshow(np.log1p(c),cmap='magma',vmin=0,vmax=3);ax.set_title(f'{t:.2f}–{t+.1:.2f} s, N={len(s)}');ax.set_xticks([]);ax.set_yticks([])
for ax in list(axs.flat)[len(times):]:ax.axis('off')
fig.tight_layout();fig.savefig(out/'contact_sheet.png',dpi=130);plt.close(fig)
bins=np.bincount(ts//10000)
np.savetxt(out/'event_rate_10ms.csv',np.column_stack((np.arange(len(bins))*.01,bins)),delimiter=',',header='relative_seconds,events_10ms',comments='')
plt.figure(figsize=(14,3));plt.plot(np.arange(len(bins))*.01,bins/0.01,lw=.5);plt.xlabel('Time from first CD event (s)');plt.ylabel('Global events / s');plt.tight_layout();plt.savefig(out/'event_rate.png',dpi=150);plt.close()
flat=counts.ravel();idx=np.argsort(flat)[-30:][::-1]
qc=dict(t0_us=t0,duration_s=float(ts[-1]/1e6),events=len(e),on_fraction=float(e['p'].mean()),top30_pixel_fraction=float(flat[idx].sum()/len(e)),top_pixels=[dict(x=int(i%w),y=int(i//w),count=int(flat[i])) for i in idx])
(out/'summary.json').write_text(json.dumps(qc,indent=2));print(json.dumps(qc,indent=2))
