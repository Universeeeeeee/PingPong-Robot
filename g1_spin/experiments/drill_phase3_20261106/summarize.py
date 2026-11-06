"""Report fixed 324 states, failure-only flow costs, raw recurrence ambiguity."""
import json
from collections import Counter
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import OUT,load,select,picture,write_json,table

def angle(a,b,unoriented=False):
    if a is None or b is None or min(np.linalg.norm(a),np.linalg.norm(b))<1e-9:return None
    value=np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(abs(value) if unoriented else value,-1,1))))

def median(values):
    values=[v for v in values if v is not None];return float(np.median(values)) if values else None

def shaft(raw,center,radius,direction):
    direction=np.array(direction,float);direction/=np.linalg.norm(direction);perp=np.array([-direction[1],direction[0]])
    xy=np.column_stack([raw['x'],raw['y']]).astype(float)-center;along=xy@direction;across=xy@perp;pts=[]
    for low in np.arange(radius+15,radius+170,12):
        keep=(along>=low)&(along<low+12)&(abs(across)<25)
        if keep.sum()>=3:pts.append([np.median(along[keep]),np.median(across[keep])])
    record=dict(status='insufficient_shaft_support',points=len(pts),span_px=None,rms_px=None,direction=None)
    if len(pts)<6:return record
    pts=np.array(pts);span=float(np.ptp(pts[:,0]));record['span_px']=span
    if span<70:return record
    result=least_squares(lambda p:p[0]*pts[:,0]+p[1]-pts[:,1],[0.,0.],loss='soft_l1')
    rms=float(np.sqrt(np.mean((result.x[0]*pts[:,0]+result.x[1]-pts[:,1])**2)))
    projected=direction+result.x[0]*perp;projected/=np.linalg.norm(projected)
    record.update(status='accepted_projection_candidate' if rms<=3 else 'large_shaft_line_residual',rms_px=rms,direction=projected.tolist())
    return record

def grouped(rows,key):
    output=[]
    for value in dict.fromkeys(r[key] for r in rows):
        for method in ['G','A','S']:
            rr=[r for r in rows if r[key]==value and r['method']==method];eligible=[r for r in rr if r['geometry_status']=='available' and r['event_count']>=100]
            valid=[r for r in rr if r['measurement_valid']]
            output.append({key:value,'method':method,'planned':len(rr),'geometry_available':sum(r['geometry_status']=='available' for r in rr),
                'input_eligible':len(eligible),'initialization_available':sum(r['init_available'] is True for r in rr) if method!='G' else None,
                'search_executed':sum(r['search'] is not None and r['search']['nm_calls']>0 for r in rr),
                'common_nonzero_evidence':sum(r['common_quality']['passed'] for r in rr),'measurement_valid':len(valid),
                'native_statuses':dict(Counter(r['native_status'] for r in rr)),'failures':dict(Counter(r['failure_reason'] for r in rr if r['failure_reason'])),
                'measurement_rpm_median':median([r['measurement_rpm'] for r in valid]),
                'candidate_rpm_median_including_rejected':median([r['candidate_rpm'] for r in rr]),
                'pipeline_s_median_on_input_eligible':median([r['pipeline_s'] for r in eligible]),
                'pipeline_s_median_on_measurements':median([r['pipeline_s'] for r in valid]),
                'cost_note':'G search cost' if method=='G' else 'Failure-only plane extraction+LS cost; no real A/S searches or speed comparison'})
    return output

def main():
    rows=[json.loads(l) for l in (OUT/'results.jsonl').read_text().splitlines()];audit=json.loads((OUT/'window_audit.json').read_text())
    for key in ['run','condition','geometry_source','partition']:
        summary=grouped(rows,key);write_json(OUT/f'summary_{key}.json',summary);table(OUT/f'summary_{key}.csv',summary)
    optimizers=[]
    for r in rows:
        s=r['search'];optimizers.append(dict(window=r['window'],method=r['method'],failure_reason=r['failure_reason'],
            raw_optimizer_candidate_omega=None if s is None else s['optimizer_candidate_omega'],raw_optimizer_candidate_contrast=None if s is None else s['optimizer_candidate_contrast'],
            zero_contrast=None if s is None else s['zero_contrast'],zero_won=r['zero_won'],selected_candidate_omega=r['candidate_omega'],
            selected_candidate_rpm=r['candidate_rpm'],measurement_valid=r['measurement_valid'],measurement_rpm=r['measurement_rpm'],
            nm_calls=0 if s is None else s['nm_calls'],objective_calls=0 if s is None else s['objective_calls'],contrast_calls=0 if s is None else s['contrast_calls']))
    table(OUT/'optimizer_candidates.csv',optimizers)
    failures=Counter();seed_counts=0
    for r in rows:
        if r['method']=='A' and r['flow_metadata']:
            failures.update(r['flow_metadata']['failures']);seed_counts+=r['flow_metadata']['seeds_attempted']
    eligible=[w for w in audit if w['geometry_status']=='available' and w['event_count']>=100]
    flow_summary=dict(planned=108,input_eligible=len(eligible),flow_seeds=seed_counts,seed_exit_reasons=dict(failures),
        eligible_not_downsampled=sum(w['common_count']==w['event_count'] for w in eligible),eligible_downsampled=sum(w['common_count']>w['event_count'] for w in eligible),
        A_init_available=0,S_init_available=0,real_A_searches=0,real_S_searches=0,
        caveat='All local-plane sets were empty. Real-data A/S direction, RPM, final accuracy and successful-pipeline speed are not identifiable.')
    write_json(OUT/'flow_diagnosis.json',flow_summary)
    # Independent fixed-period peaks remain unchanged; only now compare with G.
    period=[]
    for p in json.loads((OUT/'periodicity.json').read_text()):
        rr=[r for r in rows if r['method']=='G' and r['run']==p['run'] and r['start_us']>=round(p['start_s']*1e6) and r['end_us']<=round(p['end_s']*1e6) and r['measurement_valid']]
        gmed=median([r['measurement_rpm'] for r in rr]);chosen=p['selected_by_similarity_only']
        period.append(dict(attempt=p['attempt'],selected_peak_unchanged=chosen,all_peaks=p['all_peaks'],half_double_peak_pairs=p['half_double_peak_pairs'],
            G_valid_windows=len(rr),G_window_ids=[r['window'] for r in rr],G_candidate_rpm_median=gmed,
            primary_peak_rpm_one_repeat_per_rotation=None if chosen is None else chosen['rpm_if_one_repeat_per_rotation'],
            primary_to_G_rpm_ratio=None if chosen is None or gmed is None else chosen['rpm_if_one_repeat_per_rotation']/gmed,
            statement='No closest-CMax peak selection. Half/double periods retained; no independent RPM accuracy validation.'))
    write_json(OUT/'periodicity_comparison.json',period);table(OUT/'periodicity_comparison.csv',[{k:v for k,v in p.items() if k not in ['all_peaks','half_double_peak_pairs']} for p in period])
    # Two-dimensional shaft projection only; no hidden initialization input.
    projection=[];cfg=json.loads((OUT/'config.json').read_text())
    for run in [1,3,4]:
        e,t0=load(run)
        for w in [w for w in audit if w['run']==run]:
            r=next(r for r in rows if r['window']==w['window'] and r['method']=='G')
            pr=dict(window=w['window'],run=run,condition=w['condition'],status='geometry_unavailable',points=0,span_px=None,rms_px=None,
                direction=None,G_projection_angle_deg=None,measurement_valid=r['measurement_valid'])
            if w['geometry_status']=='available':
                d=np.load(OUT/w['input_file']);q,_=select(e,t0,w['start_us'],10000)
                direction=[0,1] if run in [1,4] else [np.cos(np.deg2rad(35)),np.sin(np.deg2rad(35))]
                if w['geometry_source']=='phase3_new_block_candidate':
                    c=next(c for c in cfg if c['run']==run and c['block_s']==w['block_s']);theta=np.deg2rad(np.mean(c['sector_deg']));direction=[np.cos(theta),np.sin(theta)]
                pr.update(shaft(q,d['center'],w['radius_px'],direction))
                if pr['status']=='accepted_projection_candidate' and r['measurement_valid'] and np.linalg.norm(r['measurement_axis'][:2])>.2:
                    pr['G_projection_angle_deg']=angle(np.array(r['measurement_axis'][:2]),np.array(pr['direction']),True)
            projection.append(pr)
    table(OUT/'shaft_projection.csv',projection);write_json(OUT/'shaft_projection.json',projection)
    ref=[r for r in rows if r['run']==1 and r['method']=='G' and r['measurement_valid']]
    mean=np.mean([r['measurement_axis'] for r in ref],axis=0);mean/=np.linalg.norm(mean)
    reverse=[]
    for r in [r for r in rows if r['condition']=='reverse_fixed' and r['method']=='G' and r['measurement_valid']]:
        reverse.append(dict(window=r['window'],wy=r['measurement_omega'][1],directed_angle_to_run1_mean_deg=angle(r['measurement_axis'],mean),
            all_run1_wy_opposite=all(r['measurement_omega'][1]*a['measurement_omega'][1]<0 for a in ref)))
    write_json(OUT/'direction_checks.json',dict(run1_reference_windows=len(ref),run1_mean_directed_axis=mean.tolist(),
        run2_abs_z_median=median([abs(r['measurement_axis'][2]) for r in rows if r['run']==2 and r['method']=='G' and r['measurement_valid']]),
        run3_abs_z_median=median([abs(r['measurement_axis'][2]) for r in rows if r['run']==3 and r['method']=='G' and r['measurement_valid']]),
        run4_fixed_reverse=reverse,caveat='Qualitative and 2D projection support only; A/S have no real measurements to compare.'))
    # All 24 new raw windows, circles only when geometry exists.
    for run in [4,5]:
        e,t0=load(run);ww=[w for w in audit if w['run']==run and w['geometry_source']=='phase3_new_block_candidate']
        fig,axes=plt.subplots(len(ww)//3,3,figsize=(12,3.8*(len(ww)//3)))
        for ax,w in zip(axes.flat,ww):
            q,_=select(e,t0,w['start_us'],10000);c=next(c for c in cfg if c['run']==run and c['block_s']==w['block_s'])
            center=np.array(w['center']) if w['center'] is not None else np.array(c['seed']);box=np.r_[np.floor(center-110).astype(int),np.ceil(center+110).astype(int)]
            ax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]])
            if w['geometry_status']=='available':ax.add_patch(plt.Circle(center,w['radius_px'],fill=False,color='lime',lw=.7))
            r=next(r for r in rows if r['window']==w['window'] and r['method']=='G')
            ax.set_title(f'{w["start_us"]/1e6:.3f}s r{w["radius_px"]:g} | '+(r['failure_reason'] or 'G candidate'),fontsize=9)
            ax.tick_params(labelsize=7)
        fig.suptitle(f'run{run} all new fixed windows | raw10ms and B1 geometry hypotheses | no spin warping')
        fig.tight_layout(rect=[0,0,1,.965]);fig.savefig(OUT/f'run{run}_new_windows.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(13,9),sharex=True)
    for run,color in [(1,'#007a9a'),(2,'#6f4a9c'),(3,'#318447'),(4,'#d56e23')]:
        rr=[r for r in rows if r['method']=='G' and r['run']==run and r['measurement_valid']]
        for dim,ax in enumerate(axes):ax.scatter([r['start_us']/1e6 for r in rr],[r['measurement_omega'][dim] for r in rr],s=23,c=color,label=f'run{run}')
    for dim,ax in enumerate(axes):ax.set_ylabel(['wx','wy','wz'][dim]+' rad/s');ax.grid(alpha=.2);ax.legend(ncol=4)
    axes[-1].set_xlabel('Per-recording time from first CD event (s)')
    fig.suptitle('G accepted angular-motion candidates only | A/S initialization coverage=0 | no gap filling')
    fig.tight_layout();fig.savefig(OUT/'G_components.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,5));conditions=list(dict.fromkeys(r['condition'] for r in rows));xx=np.arange(len(conditions))
    for m,offset,color in [('G',-.22,'#007a9a'),('A',0,'#d56e23'),('S',.22,'#6f4a9c')]:
        vals=[sum(r['measurement_valid'] for r in rows if r['condition']==c and r['method']==m) for c in conditions]
        ax.bar(xx+offset,vals,.22,label=m,color=color)
    ax.set_xticks(xx,conditions,rotation=20,ha='right');ax.set_ylabel('Accepted nonzero candidate windows');ax.legend();ax.grid(axis='y',alpha=.2)
    ax.set_title('All fixed input windows retained; A/S=0 due to empty local-plane sets, no fallback')
    fig.tight_layout();fig.savefig(OUT/'coverage.png',dpi=150);plt.close(fig)
    print('Wrote all state/source/condition summaries, original optimizer candidates, recurrence comparison, projections and raw evidence')

if __name__=='__main__':main()
