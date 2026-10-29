from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path('/Users/vae/Projects/PingPongRobot/g1_spin')
for name in ['221544','221756']:
 root=base/'data/decoded'/('local_analysis_'+name); out=Path(__file__).parent/name;out.mkdir(exist_ok=True)
 m=json.loads((root/'manifest.json').read_text());st=m['statistics'];t0=st['first_timestamp_us'];duration=st['duration_us']/1e6
 times=np.linspace(.3,duration-.3,12);starts=t0+(times*1e6).astype(np.int64);ends=starts+30000
 h,w=m['height'],m['width'];hist=np.zeros((h,w),np.int64);frames=np.zeros((12,h,w),np.int32);rates=np.zeros(int(duration*100)+1,np.int64)
 for part in m['parts']:
  e=np.load(root/part,mmap_mode='r');ts=e['t'];idx=e['y'].astype(np.int64)*w+e['x'];hist+=np.bincount(idx,minlength=h*w).reshape(h,w)
  bins=((ts-t0)//10000).astype(np.int64);rates+=np.bincount(bins,minlength=len(rates))[:len(rates)]
  for i,(s,end) in enumerate(zip(starts,ends)):
   if ts[-1]<s or ts[0]>=end:continue
   a,b=np.searchsorted(ts,[s,end]);ev=e[a:b];frames[i]+=np.bincount(ev['y'].astype(np.int64)*w+ev['x'],minlength=h*w).reshape(h,w)
 fig,axs=plt.subplots(3,4,figsize=(20,9))
 for ax,t,c in zip(axs.flat,times,frames):
  ax.imshow(np.log1p(c),cmap='magma',vmin=0,vmax=3);ax.set_title(f'{t:.2f}s +30ms N={c.sum()}');ax.set_xticks([]);ax.set_yticks([])
 fig.tight_layout();fig.savefig(out/'contact_sheet.png',dpi=130);plt.close(fig)
 plt.figure(figsize=(14,8));plt.imshow(np.log1p(hist),cmap='magma');plt.colorbar();plt.title(f'{name}: {duration:.2f}s full accumulation');plt.savefig(out/'accumulation.png',dpi=140);plt.close()
 plt.figure(figsize=(12,3));plt.plot(np.arange(len(rates))*.01,rates*100,lw=.4);plt.xlabel('Relative time /s');plt.ylabel('events/s');plt.tight_layout();plt.savefig(out/'event_rate.png',dpi=150);plt.close()
 np.savez_compressed(out/'survey_arrays.npz',times=times,frames=frames,histogram=hist,rates=rates)
 (out/'summary.json').write_text(json.dumps(st,indent=2));print(name,st,flush=True)
