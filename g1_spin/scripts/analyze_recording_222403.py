"""Reproduce selected real-clip analysis; no truth, SDK or Ubuntu required.

Circle annotations come from the saved upper-arc fits. Per-clip constant ROI
is an approximation, not compensation for the observed eccentric motion.
"""
import csv,json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
from PIL import Image
from g1_spin.spherical_cmax import lift,estimate,iwe,warp

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/decoded/local_analysis_222403'
OUT=ROOT/'outputs/real_222403/selected'
OUT.mkdir(parents=True,exist_ok=True)
meta=json.loads((DATA/'manifest.json').read_text())
events=np.concatenate([np.load(DATA/p) for p in meta['parts']])
t0=int(events['t'][0])
fits=np.genfromtxt(ROOT/'outputs/real_222403/qc/circle_fits.csv',delimiter=',',names=True)
period=np.genfromtxt(ROOT/'outputs/real_222403/periodicity/cycle_candidates.csv',delimiter=',',names=True)

def save_json(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False))

def sample(start,duration,roi):
    a,b=np.searchsorted(events['t'],[t0+round(start*1e6),t0+round((start+duration)*1e6)])
    raw=events[a:b]
    cx,cy,r=roi
    inside=(raw['x'].astype(float)-cx)**2+(raw['y'].astype(float)-cy)**2<(r*.8)**2
    # Avoid the visibly occluded lower part, not a guessed texture classifier.
    s=raw[inside & (raw['y']<240)]
    dt=(s['t'].astype(np.int64)-(t0+round(start*1e6)))*1e-6
    p=lift(s['x'],s['y'],cx,cy,r)
    result=estimate(p,dt,max_rpm=600)
    result.update(relative_start_s=start,duration_s=duration,roi=list(roi),roi_events=len(s),
                  physics_validated=False,validation_status='provisional_eccentric_motion_unmodeled')
    return result,p,dt

def visualize_clip(folder,start,end,roi):
    frames=[];cx,cy,r=roi
    for sec in np.arange(start,end,.005):
        a,b=np.searchsorted(events['t'],[t0+round(sec*1e6),t0+round((sec+.005)*1e6)])
        s=events[a:b];s=s[(s['x']>=485)&(s['x']<625)&(s['y']>=155)&(s['y']<265)]
        im=np.full((110,140,3),15,np.uint8)
        for pol,color in [(0,(255,150,70)),(1,(80,200,255))]:
            q=s[s['p']==pol];im[q['y']-155,q['x']-485]=color
        cv2.circle(im,(round(cx)-485,round(cy)-155),round(r),(140,140,140),1)
        cv2.circle(im,(round(cx)-485,round(cy)-155),round(.8*r),(80,220,100),1)
        cv2.line(im,(0,85),(139,85),(160,160,160),1)
        im=cv2.resize(im,(560,440),interpolation=cv2.INTER_NEAREST)
        cv2.putText(im,f'{sec:.3f}s | 5ms events | 0.25x',(8,22),0,.5,(230,230,230),1)
        frames.append(Image.fromarray(cv2.cvtColor(im,cv2.COLOR_BGR2RGB)))
    frames[0].save(folder/'events.gif',save_all=True,append_images=frames[1:],duration=20,loop=0)
    frames[len(frames)//2].save(folder/'roi.png')

summaries=[]
for label,start,end in [('early',1.6,1.9),('middle',6.2,6.5),('late',11.2,11.5)]:
    folder=OUT/label;folder.mkdir(exist_ok=True)
    f=fits[(fits['relative_s']>=start-1e-6)&(fits['relative_s']<end-1e-6)]
    roi=tuple(float(np.median(f[k])) for k in ['cx','cy','radius'])
    drift=float(np.hypot(np.ptp(f['cx']),np.ptp(f['cy'])))
    results=[];best=None
    for sec in np.arange(start,end-.001,.05):
        result,p,dt=sample(float(sec),.05,roi)
        results.append(result)
        if best is None or (result['valid'] and not best[0]['valid']):best=(result,p,dt)
        print(label,round(sec,2),result['spin_rpm'],result['failure_reason'],flush=True)
    # Sensitivity is a robustness check, not a confidence interval.
    sensitivity=[]
    for offset in [(0,0,0),(-2,0,0),(2,0,0),(0,-2,0),(0,2,0),(0,0,-2),(0,0,2)]:
        perturbed=tuple(v+d for v,d in zip(roi,offset))
        q,_,_=sample(start+.1,.05,perturbed)
        sensitivity.append(dict(offset=list(offset),valid=q['valid'],spin_rpm=q['spin_rpm'],axis_camera=q['axis_camera'],failure_reason=q['failure_reason']))
    valid=[r for r in results if r['valid']]
    speeds=[r['spin_rpm'] for r in valid]
    sv=[r['spin_rpm'] for r in sensitivity if r['valid']]
    axes=np.array([r['axis_camera'] for r in valid])
    axis_center=np.mean(axes,axis=0) if len(axes) else None
    if axis_center is not None:axis_center/=np.linalg.norm(axis_center)
    axis_spread=float(np.max(np.degrees(np.arccos(np.clip(axes@axis_center,-1,1))))) if len(axes) else None
    mid=(start+end)/2
    summary=dict(clip=label,start_s=start,end_s=end,roi=roi,excluded_y_ge=240,inner_radius_ratio=.8,
      fitted_center_span_px=drift,arc_rms_px=float(np.median(f['rms'])),
      computational_valid_fraction=len(valid)/len(results),candidate_rpm_median=float(np.median(speeds)) if speeds else None,
      candidate_rpm_min=min(speeds) if speeds else None,candidate_rpm_max=max(speeds) if speeds else None,
      period_candidate_rpm=float(np.interp(mid,period['time_s'],period['candidate_rpm'])),
      candidate_axis_camera=axis_center.tolist() if axis_center is not None else None,
      axis_max_deviation_deg=axis_spread,
      roi_sensitivity_rpm_range=[min(sv),max(sv)] if sv else None,
      physics_validated=False,physical_failure_reason='eccentric_motion_and_unverified_ROI_geometry',
      caution='Algorithmic convergence is not validated physical spin; periodicity is a conditional cross-check, not truth.',
      windows=results,roi_sensitivity=sensitivity)
    save_json(folder/'result.json',summary)
    with (folder/'windows.csv').open('w') as stream:
        keys=['relative_start_s','duration_s','roi_events','valid','spin_rpm','omega_rad_s','axis_camera','failure_reason','physics_validated']
        writer=csv.DictWriter(stream,fieldnames=keys,extrasaction='ignore');writer.writeheader();writer.writerows(results)
    result,p,dt=best;dt=dt-(dt[0]+dt[-1])/2 if len(dt) else dt
    after=iwe(warp(p,dt,result['omega_rad_s'])) if result['valid'] else iwe(p)
    before=iwe(p);vmax=max(before.max(),after.max())
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for ax,img,title in zip(axs,[before,after],['Uncompensated','Candidate compensation (not validated)']):
        ax.imshow(img,origin='lower',vmin=0,vmax=vmax,cmap='magma');ax.set_title(title);ax.set_xlabel('Longitude bin');ax.set_ylabel('Latitude bin')
    fig.tight_layout();fig.savefig(folder/'spherical_before_after.png',dpi=140);plt.close(fig)
    visualize_clip(folder,start,end,roi)
    summaries.append(summary)
save_json(OUT/'summary.json',summaries)
fig,ax=plt.subplots(figsize=(11,4))
ax.plot(period['time_s'],period['candidate_rpm'],'o-',ms=3,label='Pattern repetition (one repeat/revolution assumed)',color='#cf853c')
for s in summaries:
    ax.scatter([r['relative_start_s']+.025 for r in s['windows'] if r['valid']],
               [r['spin_rpm'] for r in s['windows'] if r['valid']],s=24,color='#168f9d')
    ax.axvspan(s['start_s'],s['end_s'],alpha=.12,color='#168f9d')
ax.scatter([],[],color='#168f9d',label='s-CMax candidates; eccentric motion unmodeled')
ax.axvspan(14.5,17.579495,color='gray',alpha=.13,label='Insufficient reliable repetition evidence')
ax.set(xlabel='Seconds from first event',ylabel='Conditional candidate RPM',title='Real recording: slowdown evidence, no independent ground truth')
ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'candidate_comparison.png',dpi=160);plt.close(fig)
print('Analysis saved',OUT,flush=True)
