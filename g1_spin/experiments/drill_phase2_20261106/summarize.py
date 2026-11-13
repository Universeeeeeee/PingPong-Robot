"""Summarize the fixed 46 windows; never rerun or select the estimator."""
import csv
import json
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import OUT,OLD,CONFIG,load,select,picture,table,write_json,ALPHAS
from diagnostics import shaft

def angle(a,b,unoriented=False):
    if a is None or b is None or min(np.linalg.norm(a),np.linalg.norm(b))<1e-9:return None
    dot=np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)
    return float(np.rad2deg(np.arccos(np.clip(abs(dot) if unoriented else dot,-1,1))))

def metric(values,method='median'):
    values=[v for v in values if v is not None]
    return None if not values else float(np.median(values) if method=='median' else max(values))

def main():
    rows=[json.loads(l) for l in (OUT/'results.jsonl').read_text().splitlines()]
    audit=json.loads((OUT/'window_audit.json').read_text());summary=[];paired=[];projection=[];ambiguity=[]
    geometry_support=[]
    for w in audit:
        g=dict(window=w['window'],run=w['run'],geometry_status=w['geometry_status'],radius_px=w['radius_px'],
            B1_center=None,B2_first_center=None,B2_last_center=None,B0_center=None,
            used_geometry_observations=0,min_support_bins=None,max_ring_rms_px=None,geometry_failure=w['failure'])
        if w['geometry_status']!='unavailable':
            d=np.load(OUT/'windows'/f'{w["window"]}.npz');b1='A0' if w['kind']=='ablation' else 'B1';b2='A1' if w['kind']=='ablation' else 'B2'
            g.update(B1_center=d[b1][0].tolist(),B2_first_center=d[b2][0].tolist(),B2_last_center=d[b2][-1].tolist(),
                B0_center=d['B0'][0].tolist() if 'B0' in d.files else None)
            if w['run'] in [1,5]:
                source=list(csv.DictReader((OLD/('rotation_geometry.csv' if w['run']==1 else 'control_geometry.csv')).open()))
            else:source=[r for r in json.loads((OUT/'geometry.json').read_text()) if r['run']==w['run']]
            source=[r for r in source if str(r['geometry_valid'])=='True'];source.sort(key=lambda r:float(r['time_s']))
            ts=np.array([round(float(r['time_s'])*1e6) for r in source]);query=np.r_[d['time_us'],w['reference_us'],w['start_us'],w['end_us']]
            if w['anchor_us'] is not None:query=np.r_[query,w['anchor_us']]
            right=np.searchsorted(ts,query,side='left');used=set(right.tolist());used.update((right[ts[right]!=query]-1).tolist())
            supports=[source[i] for i in used]
            g.update(used_geometry_observations=len(supports),min_support_bins=min(int(r['angular_bins']) for r in supports),
                max_ring_rms_px=max(float(r['ring_rms_px']) for r in supports))
        geometry_support.append(g)
    table(OUT/'window_geometry_support.csv',geometry_support);write_json(OUT/'window_geometry_support.json',geometry_support)
    for run in [2,3,4]:
        for part in ['all','development','evaluation']:
            for v in ['B1','B2']:
                rr=[r for r in rows if r['run']==run and r['variant']==v and (part=='all' or r['partition']==part)]
                valid=[r for r in rr if r['valid']]
                summary.append(dict(run=run,partition=part,variant=v,planned=len(rr),
                    geometry_available=sum(r['geometry_status']!='unavailable' for r in rr),valid=len(valid),
                    failures=dict(Counter(r['failure_reason'] for r in rr if r['failure_reason'])),
                    rpm_median=metric([r['rpm'] for r in valid]),abs_axis_z_median=metric([abs(r['axis'][2]) for r in valid]),
                    mean_directed_axis=np.mean([r['axis'] for r in valid],axis=0).tolist() if valid else None,
                    runtime_median_s=metric([r['runtime_s'] for r in rr if r['runtime_s']>0])))
        e,t0=load(run);fig,axes=plt.subplots(3,4,figsize=(15,11))
        for ax,w in zip(axes.flat,[w for w in audit if w['run']==run]):
            q,_=select(e,t0,w['start_us'],10000);b1=next(r for r in rows if r['window']==w['window'] and r['variant']=='B1');b2=next(r for r in rows if r['window']==w['window'] and r['variant']=='B2')
            pr=dict(window=w['window'],run=run,shaft_status='geometry_unavailable',points=0,span_px=None,rms_px=None,
                direction=None,B1_projection_angle_deg=None,B2_projection_angle_deg=None)
            center=np.array(CONFIG[run]['seed'])
            if w['geometry_status']!='unavailable':
                d=np.load(OUT/'windows'/f'{w["window"]}.npz');center=d['B1'][0]
                if run==2:pr['shaft_status']='no_reliable_visible_shaft'
                else:
                    s=shaft(q,center,w['radius_px'],[np.cos(np.deg2rad(35)),np.sin(np.deg2rad(35))] if run==3 else [0,1])
                    pr.update(shaft_status=s['status'],**{k:s[k] for k in ['points','span_px','rms_px','direction']})
                    if s['status']=='accepted_projection_candidate':
                        for r in [b1,b2]:
                            if r['valid'] and np.linalg.norm(r['axis'][:2])>.2:
                                pr[r['variant']+'_projection_angle_deg']=angle(np.array(r['axis'][:2]),np.array(s['direction']),True)
            projection.append(pr)
            box=np.r_[np.floor(center-[90,90]).astype(int),np.ceil(center+[150 if run==3 else 90,150]).astype(int)]
            ax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]])
            if w['geometry_status']!='unavailable':
                ax.add_patch(plt.Circle(center,w['radius_px'],fill=False,color='lime',lw=.8))
                ax.plot(d['B2'][[0,-1],0],d['B2'][[0,-1],1],'.-',color='magenta',lw=.8,ms=3)
                if pr['shaft_status']=='accepted_projection_candidate':
                    direction=np.array(pr['direction']);line=center+np.array([w['radius_px']+15,w['radius_px']+150])[:,None]*direction
                    ax.plot(line[:,0],line[:,1],color='white',lw=.8)
            ax.set_title(f'{w["start_us"]/1e6:.3f}s {w["partition"][:3]}\n'+('B1 circle / B2 trace' if w['geometry_status']!='unavailable' else 'NO GEOMETRY'),fontsize=9)
            ax.tick_params(labelsize=7)
            paired.append(dict(window=w['window'],run=run,partition=w['partition'],geometry_status=w['geometry_status'],
                B1_valid=b1['valid'],B2_valid=b2['valid'],axis_difference_deg=angle(b1['axis'],b2['axis']) if b1['valid'] and b2['valid'] else None,
                rpm_difference=b2['rpm']-b1['rpm'] if b1['valid'] and b2['valid'] else None,
                candidate_axis_difference_deg=angle(b1['axis'],b2['axis']),
                B1_failure=b1['failure_reason'],B2_failure=b2['failure_reason']))
        fig.suptitle(f'run{run}: all 12 fixed raw 10ms windows | geometry hypotheses only; white=external shaft projection candidate')
        fig.tight_layout();fig.savefig(OUT/f'run{run}_geometry_evidence.png',dpi=150);plt.close(fig)
    for r in rows:
        candidates=r['optimizer_candidates'];best=max([c['contrast'] for c in candidates],default=0.)
        close=[c for c in candidates if c['contrast']>=best*.98 and np.linalg.norm(c['omega'])*.01>=.03]
        angles=[angle(a['omega'],b['omega']) for i,a in enumerate(close) for b in close[i+1:]]
        ambiguity.append(dict(window=r['window'],variant=r['variant'],valid=r['valid'],
            candidates=len(candidates),near_score_nontrivial_solutions=len(close),max_pairwise_axis_deg=max(angles) if angles else None))
    table(OUT/'paired_comparison.csv',paired);table(OUT/'shaft_projection.csv',projection);table(OUT/'candidate_ambiguity.csv',ambiguity)
    write_json(OUT/'summary.json',summary)
    # Ablation diagnostics include rejected candidates and do not choose alpha.
    predictions=json.loads((OUT/'velocity_predictions.json').read_text());ablation=[]
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,p in zip(axes.flat,predictions):
        rr=[r for r in rows if r['window']==p['window']];by={r['variant']:r for r in rr};base=np.array(by['A0']['omega']);actual=np.array(by['A1']['omega'])-base
        out=dict(**p,actual_delta_omega=actual.tolist(),x_near_prediction_residual=float(actual[0]-p['near_axis_delta_omega'][0]),
            full_prediction_residual_norm=float(np.linalg.norm(actual-np.array(p['full_sphere_linear_delta_omega']))),
            A0_A1_valid=[by['A0']['valid'],by['A1']['valid']],axis_A0_A1_difference_deg=angle(by['A0']['axis'],by['A1']['axis']),
            alpha_valid={v:r['valid'] for v,r in by.items()},alpha_gain={v:r['gain'] for v,r in by.items()})
        ablation.append(out);aa=list(ALPHAS.values())
        ax.plot(aa,[by[v]['omega'][0] for v in ALPHAS],'.-',label='optimizer candidate wx')
        ax.plot(aa,[base[0]+a*p['near_axis_delta_omega'][0] for a in aa],'--',label='near-axis prediction')
        ax.plot(aa,[base[0]+a*p['full_sphere_linear_delta_omega'][0] for a in aa],':',label='sphere displacement prediction')
        ax.set(title=f'run{p["run"]} {p["start_us"]/1e6:.3f}s',xlabel='velocity factor alpha',ylabel='wx (rad/s)');ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=7);fig.suptitle('Velocity-hypothesis ablation | includes rejected optimizer candidates | no alpha selected')
    fig.tight_layout();fig.savefig(OUT/'velocity_ablation.png',dpi=150);plt.close(fig)
    write_json(OUT/'velocity_ablation.json',ablation);table(OUT/'velocity_ablation.csv',ablation)
    # Directed reversal uses the old fixed set; the recovered supplement stays separate.
    oldrows=[json.loads(l) for l in (OLD/'results.jsonl').read_text().splitlines()]
    ref=[r for r in oldrows if r['state']=='rotation' and r['variant']=='B1' and r['valid']]
    axis1=np.array([r['result']['axis_camera'] for r in ref]).mean(axis=0);axis1/=np.linalg.norm(axis1)
    reversal=[]
    for r in [r for r in rows if r['run']==4 and r['valid']]:
        reversal.append(dict(window=r['window'],variant=r['variant'],wy=r['omega'][1],
            directed_angle_to_run1_mean_deg=angle(r['axis'],axis1),angle_to_negated_run1_mean_deg=angle(r['axis'],-axis1),
            wy_reversed_relative_to_all_12_run1_B1=all(r['omega'][1]*a['wy']<0 for a in ref)))
    table(OUT/'reversal.csv',reversal)
    write_json(OUT/'reversal.json',dict(run1_reference_windows=len(ref),run1_mean_directed_axis=axis1.tolist(),rows=reversal))
    # Availability and three components; absent/invalid samples remain visibly absent.
    fig,axes=plt.subplots(4,3,figsize=(15,10),sharex='col')
    for col,run in enumerate([2,3,4]):
        for v,color,offset in [('B1','#007a9a',0),('B2','#d96719',.12)]:
            rr=[r for r in rows if r['run']==run and r['variant']==v]
            for r in rr:
                status=0 if r['geometry_status']=='unavailable' else (2 if r['valid'] else 1)
                axes[0,col].scatter(r['start_us']/1e6,status+offset,s=25,c=color,label=v if r==rr[0] else None)
            good=[r for r in rr if r['valid']]
            for dim in range(3):axes[dim+1,col].scatter([r['start_us']/1e6 for r in good],[r['omega'][dim] for r in good],s=25,c=color,label=v)
        axes[0,col].set(title=f'run{run}',yticks=[0,1,2],yticklabels=['geometry missing','rejected','accepted']);axes[-1,col].set_xlabel('Seconds from first CD event')
        for dim in range(3):axes[dim+1,col].set_ylabel(['wx','wy','wz'][dim]+' rad/s')
        for ax in axes[:,col]:ax.grid(alpha=.2)
    axes[0,0].legend();fig.suptitle('Fixed extension windows | accepted candidates only in component plots | no gap filling')
    fig.tight_layout();fig.savefig(OUT/'extension_results.png',dpi=150);plt.close(fig)
    print('Wrote summaries, all-window geometry evidence, projections, reversal, ablation and plots')

if __name__=='__main__':main()
