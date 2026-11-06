from pathlib import Path
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).parent
for name in ['221544','221756']:
 d=np.load(base/name/'survey_arrays.npz');frames=d['frames'];times=d['times'];hist=d['histogram'];
 fig,axs=plt.subplots(3,4,figsize=(16,10))
 rows=[]
 for ax,t,c in zip(axs.flat,times,frames):
  ax.imshow(np.log1p(c[130:280,450:680]),origin='upper',extent=[450,680,280,130],cmap='magma',vmin=0,vmax=4)
  ax.set_title(f'{t:.2f}s');ax.grid(alpha=.2)
  roi=c[150:250,490:650];rows.append(dict(time_s=float(t),ball_box_events=int(roi.sum()),global_events=int(c.sum()),ball_box_fraction=float(roi.sum()/c.sum())))
 fig.tight_layout();fig.savefig(base/name/'ball_crops.png',dpi=140);plt.close(fig)
 (base/name/'roi_statistics.json').write_text(json.dumps(rows,indent=2))
 print(name,rows)
