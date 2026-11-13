import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
OUT=Path(__file__).resolve().parent
rows=json.loads((OUT/'prediction/summary.json').read_text());profile=json.loads((OUT/'prediction/profile_results.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(16,5))
rr=[r for r in rows if not r['case'].startswith('synthetic')];labels=[r['case'].replace('disk_','D:').replace('drill_','R:') for r in rr]
axes[0].barh(labels,[r['future_gain_vs_zero']*100 for r in rr],color=['#237d9b' if r['native_valid'] else '#aaaaaa' for r in rr]);axes[0].set(xlabel='Future pixel-distance reduction vs zero (%)',title='Conditional prediction evidence\nGray: estimator rejected candidate');axes[0].invert_yaxis()
for name,label,color in [('disk_motion','Moving disk','#c87b25'),('drill_1','Fixed-axis drill','#237d9b')]:
 r=[x for x in profile if x['case']==name];r.sort(key=lambda x:x['delta_center_vx_px_s']);x=[x['delta_center_vx_px_s'] for x in r]
 axes[1].plot(x,[x['result']['spin_rpm'] for x in r],'o-',label=label,color=color)
 vals=np.array([x['evaluation']['future_balanced_mean'] for x in r]);nominal=vals[2]
 axes[2].plot(x,100*(vals/nominal-1),'o-',label=label,color=color)
axes[1].set(xlabel='Assumed center-velocity offset (px/s)',ylabel='Candidate RPM',title='Different RPM explanations')
axes[2].set(xlabel='Assumed center-velocity offset (px/s)',ylabel='Change from nominal prediction distance (%)',title='Almost unchanged future prediction');axes[2].axhline(0,color='gray',lw=.7)
for ax in axes:ax.grid(alpha=.2)
axes[1].legend();axes[2].legend();fig.tight_layout();fig.savefig(OUT/'strategy_evidence.png',dpi=150);plt.close(fig)
