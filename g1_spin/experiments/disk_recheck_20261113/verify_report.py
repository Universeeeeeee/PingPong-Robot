"""Independent raw-mask/geometry/metric reconstruction and reporting."""
import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];DATA=ROOT/'g1_spin/data/decoded/local_analysis_222403'
rows=[json.loads(x) for x in (OUT/'results.jsonl').read_text().splitlines()]
if (OUT/'high_motion_results.json').exists():rows+=json.loads((OUT/'high_motion_results.json').read_text())
models=json.loads((ROOT/'g1_spin/experiments/real_center_sensitivity_20261102/models.json').read_text());m=json.loads((DATA/'manifest.json').read_text());e=np.concatenate([np.load(DATA/p) for p in m['parts']]);origin=int(e['t'][0]);checks=[]
for path in sorted((OUT/'inputs').glob('*.npz')):
 rr=[r for r in rows if r['window']==path.stem];assert rr
 first=rr[0];mid=first['mid_us'];ms=first['window_ms'];begin,end=mid-ms*500,mid+ms*500
 rawids=np.flatnonzero((e['t']>=origin+begin)&(e['t']<origin+end));raw=e[rawids];t=raw['t'].astype(np.int64)-origin;pix=np.column_stack([raw['x'],raw['y']]).astype(float)
 cs={}
 for method in ['fixed50','linear50','linear100']:
  model=models[first['segment']]['100' if method=='linear100' else '50'];query=np.full(len(raw),mid) if method=='fixed50' else t
  cs[method]=np.asarray(model['center_anchor_px'])+np.outer(query/1e6-model['anchor_s'],model['velocity_px_s'])
 minr=min(r['radius'] for r in rr);keep=(pix[:,0]>=475)&(pix[:,0]<638)&(pix[:,1]>=143)&(pix[:,1]<240)
 for c in cs.values():keep&=np.sum((pix-c)**2,axis=1)<=(.8*minr)**2
 ix=np.flatnonzero(keep);data=np.load(path);assert np.array_equal(rawids[ix],data['common_ids'])
 if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
 assert np.array_equal(rawids[ix],data['ids']) and np.array_equal(raw[ix],data['events'])
 for method in cs:assert np.allclose(cs[method][ix],data[method],atol=1e-10,rtol=0)
 for r in rr:
  result=r['result'];assert r['input_count']==len(ix) and r['common_count']==keep.sum()
  if result['valid']:
   w=np.array(result['omega_rad_s']);assert np.isclose(np.linalg.norm(w)*60/(2*np.pi),result['spin_rpm']);assert np.allclose(w/np.linalg.norm(w),result['axis_camera'])
 checks.append(dict(window=path.stem,events=len(ix),raw_ids_centers_mask_verified=True))
for p,h in json.loads((OUT/'freeze.json').read_text()).items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
summary=[];bound=[];radius=[]
for segment in ['translation','early_turn','middle','late']:
 for ms in [20,50,100]:
  for method in ['fixed50','linear50','linear100']:
   rr=[r for r in rows if r['segment']==segment and r['window_ms']==ms and r['method']==method and r['radius']==48 and r['max_rpm']==600]
   good=[r['result']['spin_rpm'] for r in rr if r['result']['valid']]
   if rr:summary.append(dict(segment=segment,window_ms=ms,method=method,valid=len(good),planned=len(rr),median=float(np.median(good)) if good else None,minimum=min(good) if good else None,maximum=max(good) if good else None))
 for r in [r for r in rows if r['segment']==segment and r['max_rpm']==3500]:
  base=next(x for x in rows if x['window']==r['window'] and x['method']==r['method'] and x['radius']==48 and x['max_rpm']==600)
  bound.append(dict(segment=segment,method=r['method'],rpm600=base['result']['spin_rpm'],rpm3500=r['result']['spin_rpm'],delta=r['result']['spin_rpm']-base['result']['spin_rpm'] if r['result']['valid'] and base['result']['valid'] else None))
 for r in [r for r in rows if r['segment']==segment and r['radius']!=48]:radius.append(dict(segment=segment,radius=r['radius'],rpm=r['result']['spin_rpm']))
rec=json.loads((OUT/'recurrence.json').read_text())
if (OUT/'high_motion_recurrence/recurrence.json').exists():rec+=json.loads((OUT/'high_motion_recurrence/recurrence.json').read_text())
# Retain all recurrence peaks; never find a peak nearest CMAX.
rec_summary=[dict(recording=r['recording'],source_anchor_s=r['anchor_s'],background=r['background'],highest_peak=r['highest_peak'],all_peaks=r['peaks']) for r in rec]
(OUT/'summary.json').write_text(json.dumps(dict(cmax=summary,search_bound=bound,radius=radius,recurrence=rec_summary),indent=2))
(OUT/'verification.json').write_text(json.dumps(dict(rows=len(rows),inputs=len(checks),checks=checks,protected_hashes_unchanged=True),indent=2))
fig,axes=plt.subplots(1,2,figsize=(13,4))
for method,color in [('fixed50','#777777'),('linear50','#168f9d'),('linear100','#b56f25')]:
 for ms,marker in [(20,'x'),(50,'o'),(100,'s')]:
  rr=[r for r in rows if r['method']==method and r['max_rpm']==600 and r['radius']==48 and r['result']['valid']]
  rr=[r for r in rr if r['window_ms']==ms];axes[0].scatter([r['mid_us']/1e6 for r in rr],[r['result']['spin_rpm'] for r in rr],color=color,marker=marker,s=24,label=f'{method} {ms}ms')
axes[0].set(xlabel='Seconds since first event',ylabel='Candidate RPM',title='Eccentric disk: slowdown survives center hypotheses');axes[0].legend(fontsize=6,ncol=3);axes[0].grid(alpha=.2)
for seg,color in [('early_turn','#168f9d'),('middle','#b56f25'),('late','#8058aa')]:
 rr=[r for r in rows if r['segment']==seg and r['method']=='linear50' and r['window_ms']==50 and r['max_rpm']==600 and (r['radius']!=48 or r['mid_us']==round({'early_turn':1.75,'middle':6.35,'late':11.35}[seg]*1e6))];rr.sort(key=lambda r:r['radius']);axes[1].plot([r['radius'] for r in rr],[r['result']['spin_rpm'] for r in rr],'o-',label=seg,color=color)
axes[1].set(xlabel='Assumed radius (px)',ylabel='Candidate RPM',title='Absolute scale remains geometry-sensitive');axes[1].legend();axes[1].grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'speed_summary.png',dpi=150);plt.close(fig)
print(json.dumps(dict(rows=len(rows),input_windows=len(checks),max_search_bound_delta=max(abs(x['delta']) for x in bound if x['delta'] is not None),high_motion=[x for x in summary if x['segment']=='translation']),indent=2))
