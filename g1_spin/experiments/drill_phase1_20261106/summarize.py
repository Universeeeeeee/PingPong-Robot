"""Evidence, diagnostics and plots; does not modify or re-run the estimator."""
import csv
import json
from collections import Counter
from pathlib import Path

import cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from explore import OUT, load, select, picture

COLORS={'B0':'#888888','B1':'#007a9a','B2':'#d96719'}

def vector(row,candidate=False):
    if row['valid']:return np.array([row['wx'],row['wy'],row['wz']])
    if candidate and row['result']:
        a=row['result']['quality'].get('rejected_candidate_omega_rad_s')
        return None if a is None else np.array(a)
    return None

def angle(a,b,unoriented=False):
    if a is None or b is None or min(np.linalg.norm(a),np.linalg.norm(b))<1e-9:return None
    a=a/np.linalg.norm(a);b=b/np.linalg.norm(b)
    value=np.dot(a,b)
    return float(np.rad2deg(np.arccos(np.clip(abs(value) if unoriented else value,-1,1))))

def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def main():
    rows=[json.loads(l) for l in (OUT/'results.jsonl').read_text().splitlines()]
    audit=json.loads((OUT/'window_audit.json').read_text())
    summary=[]
    for state in ['rotation','static','translation']:
        for part in ['all','development','evaluation']:
            for v in COLORS:
                rr=[r for r in rows if r['state']==state and r['variant']==v and (part=='all' or r['partition']==part)]
                good=[r for r in rr if r['valid']];candidate=[r['candidate_rpm'] for r in rr if r['candidate_rpm'] is not None]
                timings=[r['runtime_s'] for r in rr if r['runtime_s']>0]
                axes=np.array([vector(r)/np.linalg.norm(vector(r)) for r in good])
                mean=axes.mean(axis=0) if len(axes) else None
                deviations=[angle(a,mean) for a in axes] if len(axes) else []
                summary.append(dict(state=state,partition=part,variant=v,planned_windows=len(rr),
                    geometry_available=sum(r['geometry_status']!='unavailable' for r in rr),
                    optimizer_supported=sum(r['optimizer_count']>=100 for r in rr),valid_nonzero=len(good),
                    failures=dict(Counter(r['failure_reason'] for r in rr if r['failure_reason'])),
                    valid_rpm_median=float(np.median([r['rpm'] for r in good])) if good else None,
                    candidate_rpm_median=float(np.median(candidate)) if candidate else None,
                    candidate_rpm_p95=float(np.percentile(candidate,95)) if candidate else None,
                    axis_spread_median_deg=float(np.median(deviations)) if deviations else None,
                    axis_spread_max_deg=max(deviations) if deviations else None,
                    runtime_p50_s=float(np.median(timings)) if timings else None,
                    runtime_p95_s=float(np.percentile(timings,95)) if timings else None))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    sensitivity=[]
    for r in rows:
        if r['variant'] in COLORS:continue
        base=next(q for q in rows if q['window']==r['window'] and q['variant']=='B2')
        sensitivity.append(dict(window=r['window'],state=r['state'],variant=r['variant'],valid=r['valid'],
            failure_reason=r['failure_reason'],base_valid=base['valid'],
            directed_axis_difference_deg=angle(vector(r),vector(base)),
            candidate_vector_difference_rad_s=float(np.linalg.norm(vector(r,True)-vector(base,True)))
                if vector(r,True) is not None and vector(base,True) is not None else None,
            valid_rpm_difference=r['rpm']-base['rpm'] if r['valid'] and base['valid'] else None,
            candidate_rpm=r['candidate_rpm']))
    write_csv(OUT/'sensitivity.csv',sensitivity)
    ambiguous=[]
    for r in rows:
        if r['variant'] not in COLORS or not r['optimizer_candidates']:continue
        candidates=r['optimizer_candidates'];best=max(c['contrast'] for c in candidates)
        close=[c for c in candidates if c['contrast']>=.98*best and np.linalg.norm(c['omega'])*.01>=.03]
        angles=[angle(np.array(a['omega']),np.array(b['omega'])) for i,a in enumerate(close) for b in close[i+1:]]
        ambiguous.append(dict(window=r['window'],state=r['state'],variant=r['variant'],
            returned_valid=r['valid'],total_solutions=len(candidates),near_score_nontrivial_solutions=len(close),
            max_pairwise_directed_axis_deg=max(angles) if angles else None,
            note='2% score and 0.03rad displacement diagnostic, not a new acceptance gate'))
    write_csv(OUT/'candidate_ambiguity.csv',ambiguous)
    # Independent visible shaft projection diagnostic, never fed into spin fitting.
    e,t0=load('rotation');shaft=[]
    for w in [w for w in audit if w['state']=='rotation' and w['geometry_status']!='unavailable']:
        d=np.load(OUT/'windows'/f'{w["window"]}.npz');center=d['B1'][0]
        raw=select(e,t0,w['start_us']/1e6,10)
        q=raw[(abs(raw['x'].astype(float)-center[0])<40)&
              (raw['y']>center[1]+75)&(raw['y']<center[1]+230)]
        pts=[]
        for low in np.arange(center[1]+75,center[1]+230,12):
            a=q[(q['y']>=low)&(q['y']<low+12)]
            if len(a)>=3:pts.append([np.median(a['x']),np.median(a['y'])])
        pts=np.array(pts);record=dict(window=w['window'],shaft_status='insufficient_support',slope=None,rms_px=None,
                                    y_span_px=None,B0_projection_angle_deg=None,B1_projection_angle_deg=None,B2_projection_angle_deg=None)
        if len(pts)>=6 and np.ptp(pts[:,1])>=70:
            yy=pts[:,1]-center[1];xx=pts[:,0]-center[0]
            fit=least_squares(lambda p:p[0]*yy+p[1]-xx,[0,0],loss='soft_l1')
            residual=fit.x[0]*yy+fit.x[1]-xx;rms=float(np.sqrt(np.mean(residual**2)))
            record.update(slope=float(fit.x[0]),rms_px=rms,y_span_px=float(np.ptp(pts[:,1])),
                          shaft_status='accepted_projection_candidate' if rms<=3 else 'large_line_residual')
            if rms<=3:
                for v in COLORS:
                    r=next(r for r in rows if r['window']==w['window'] and r['variant']==v)
                    a=vector(r)
                    if a is not None and np.linalg.norm(a[:2])/np.linalg.norm(a)>.2:
                        record[v+'_projection_angle_deg']=angle(a[:2],np.array([fit.x[0],1]),True)
        shaft.append(record)
    write_csv(OUT/'shaft_projection.csv',shaft)
    # Main outcome plot: no interpolation through absent or invalid windows.
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    for v,color in COLORS.items():
        rr=[r for r in rows if r['state']=='rotation' and r['variant']==v and r['valid']]
        axes[0,0].scatter([r['start_us']/1e6 for r in rr],[r['rpm'] for r in rr],s=28,c=color,label=v)
        ss=[r for r in shaft if r[v+'_projection_angle_deg'] is not None]
        axes[0,1].scatter([int(r['window'].split('_')[-1])/1e6 for r in ss],
                          [r[v+'_projection_angle_deg'] for r in ss],s=28,c=color,label=v)
        rr=[r for r in rows if r['state']=='translation' and r['variant']==v and r['candidate_rpm'] is not None]
        axes[1,0].scatter([r['start_us']/1e6 for r in rr],[r['candidate_rpm'] for r in rr],s=25,c=color,marker='x',label=v)
        good=[r for r in rr if r['valid']]
        axes[1,0].scatter([r['start_us']/1e6 for r in good],[r['rpm'] for r in good],s=110,facecolors='none',edgecolors=color)
        axes[1,1].scatter([r['start_us']/1e6 for r in rr],[100*r['gain'] for r in rr],s=25,c=color,label=v)
    axes[0,0].set(title='Rotation: accepted candidate speed (no RPM truth)',ylabel='RPM',xlabel='First recording time (s)')
    axes[0,1].set(title='Rotation: angle to observed shaft projection',ylabel='Unoriented 2D angle (deg), not 3D error',xlabel='First recording time (s)')
    axes[1,0].set(title='Translation: optimizer candidates; circles = accepted',ylabel='Candidate RPM, rejected values are NOT measurements',xlabel='Fifth recording time (s)')
    axes[1,1].set(title='Translation: gain over zero-spin score',ylabel='Contrast gain (%)',xlabel='Fifth recording time (s)')
    axes[1,1].axhline(8,color='black',ls='--',lw=1,label='Existing 8% gate')
    for ax in axes.flat:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Same-event B0/B1/B2 comparison | unvalidated center hypotheses | no spin smoothing')
    fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for v,color in COLORS.items():
        rr=[r for r in rows if r['state']=='rotation' and r['variant']==v and r['valid']]
        for dim,ax in zip(['wx','wy','wz'],axes):ax.scatter([r['start_us']/1e6 for r in rr],[r[dim] for r in rr],c=color,label=v,s=25);ax.set_ylabel(dim+' (rad/s)')
    axes[-1].set_xlabel('First recording time (s)')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Rotation angular velocity components; missing windows are not interpolated')
    fig.tight_layout();fig.savefig(OUT/'omega_components.png',dpi=140);plt.close(fig)
    # Same raw events, before/after translational alignment only, no spin warping.
    names=['rotation_01425000','rotation_03665000','translation_04165000','translation_06425000']
    fig,axes=plt.subplots(4,3,figsize=(12,13))
    for i,name in enumerate(names):
        d=np.load(OUT/'windows'/f'{name}.npz');tt=(d['time_us']-d['time_us'][0])/1000
        for j,v in enumerate(COLORS):
            xy=d['xy']-d[v];ax=axes[i,j]
            ax.scatter(xy[:,0],xy[:,1],c=tt,cmap='viridis',s=.5,vmin=0,vmax=10)
            ax.set(xlim=(-50,50),ylim=(50,-50),aspect='equal',title=f'{name}\n{v}',xlabel='center-relative x',ylabel='center-relative y')
    fig.suptitle('Identical event IDs, color=time (0-10ms); center subtraction only, no rotation compensation')
    fig.tight_layout();fig.savefig(OUT/'event_alignment.png',dpi=140);plt.close(fig)
    # Geometry overlay videos: raw 5ms samples every20ms, not a spin truth visualization.
    video_checks=[]
    for key in ['rotation','control']:
        e,t0=load(key);duration=(int(e['t'][-1])-t0)*1e-6
        gg=list(csv.DictReader((OUT/f'{key}_geometry.csv').open()));times=np.array([float(r['time_s']) for r in gg])
        path=OUT/f'{key}_geometry_preview.mp4';writer=cv2.VideoWriter(str(path),cv2.VideoWriter_fourcc(*'mp4v'),50,(1280,720));assert writer.isOpened();n=0
        for s in np.arange(.005,duration-.005,.02):
            im=cv2.cvtColor(picture(select(e,t0,s-.0025,5)),cv2.COLOR_RGB2BGR)
            g=gg[int(np.argmin(abs(times-s)))];valid=g['geometry_valid']=='True'
            if valid:cv2.circle(im,(round(float(g['cx'])),round(float(g['cy']))),60,(60,230,60),1)
            status='geometry candidate, not truth' if valid else 'NO GEOMETRY (not filled)'
            cv2.putText(im,f'{key} {s:.3f}s | 5ms / 20ms sampled | {status}',(10,25),0,.6,(230,230,230),1)
            writer.write(im);n+=1
        writer.release();cap=cv2.VideoCapture(str(path));count=0
        while True:
            ok,frame=cap.read()
            if not ok:break
            assert frame.shape[:2]==(720,1280);count+=1
        cap.release();assert count==n
        video_checks.append(dict(file=path.name,frames=n,redecoded=count,fps=50,sample_step_ms=20,window_ms=5))
    (OUT/'video_verification.json').write_text(json.dumps(video_checks,indent=2))
    print('Summary, candidate audit, sensitivity, projection diagnostics, plots and verified videos written')

if __name__=='__main__':main()
