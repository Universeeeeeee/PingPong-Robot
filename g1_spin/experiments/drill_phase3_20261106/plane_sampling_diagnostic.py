"""Same three windows/physical seeds, full inner support versus 5000-event input."""
import importlib.util
import json
import sys
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import OUT,ROOT,P1,load,digest,write_json,table

def audit_seeds(xy,times,polarity,seed_indices,physical_ids,center,radius,window,support):
    tree=cKDTree(xy);rows=[]
    for seed,physical in zip(seed_indices,physical_ids):
        indices=np.asarray(tree.query_ball_point(xy[seed],4),int)
        indices=indices[(abs(times[indices]-times[seed])<=.002)&(polarity[indices]==polarity[seed])]
        row=dict(window=window,support=support,seed_index=int(seed),physical_event_id=int(physical),neighbors=len(indices),
            spatial_ratio=None,time_std_s=None,plane_rms_s=None,rms_over_std=None,threshold_s=None,gradient_norm_s_per_px=None,status='too_few_neighbors')
        if len(indices)>=10:
            positions=xy[indices];c=positions.mean(axis=0);sp=positions-c;singular=np.linalg.svd(sp,compute_uv=False)
            ratio=float(singular[1]/singular[0]) if singular[0]>1e-6 else 0.
            dt=times[indices]-times[seed];std=float(np.std(dt));design=np.column_stack([sp,np.ones(len(indices))])
            coeff,_,rank,_=np.linalg.lstsq(design,dt,rcond=None);rms=float(np.sqrt(np.mean((design@coeff-dt)**2)));threshold=max(20e-6,.15*std)
            gradient=float(np.linalg.norm(coeff[:2]));local=(c-center)/radius
            status='accepted'
            if ratio<.08:status='spatial_plane_degenerate'
            elif std<1e-6:status='insufficient_local_time_span'
            elif rank<3 or rms>threshold:status='plane_residual_or_rank'
            elif np.sum(local*local)>=.9**2:status='outside_flow_inner_sphere'
            elif gradient<1e-8:status='vanishing_gradient'
            row.update(spatial_ratio=ratio,time_std_s=std,plane_rms_s=rms,rms_over_std=rms/std if std else None,
                threshold_s=threshold,gradient_norm_s_per_px=gradient,status=status)
        rows.append(row)
    return rows

def main():
    sys.path.insert(0,str(ROOT/'g1_spin/src'));from g1_spin.events import EventBatch
    spec=importlib.util.spec_from_file_location('frozen_flow_readonly',P1/'frozen/optical_flow.py');frozen=importlib.util.module_from_spec(spec);spec.loader.exec_module(frozen)
    assert digest(ROOT/'g1_spin/src/g1_spin/optical_flow.py')==digest(P1/'frozen/optical_flow.py')
    audit=json.loads((OUT/'window_audit.json').read_text());records=[];summaries=[];fig,axes=plt.subplots(2,3,figsize=(15,8))
    names=['run1_01495000','run2_03425000','run4_07495000']
    for col,name in enumerate(names):
        w=next(w for w in audit if w['window']==name);d=np.load(OUT/w['input_file']);e,origin=load(w['run'])
        full=e[d['common_ids']];full_xy=np.column_stack([full['x'],full['y']]).astype(float);full_t=full['t'].astype(np.int64)-origin
        seed=np.linspace(0,len(d['ids'])-1,min(300,len(d['ids'])),dtype=int);physical=d['ids'][seed]
        fullseed=np.searchsorted(d['common_ids'],physical);assert np.array_equal(d['common_ids'][fullseed],physical)
        native=[]
        for row,support,xy,times,pol,seed_indices,ids in [(0,'selected_max5000',d['xy'],d['time_us'],d['polarity'],seed,d['ids']),
                (1,'full_same_inner_support',full_xy,full_t,full['p'],fullseed,d['common_ids'])]:
            seconds=(times-int(times[0]))*1e-6
            metrics=audit_seeds(xy,seconds,pol,seed_indices,physical,d['center'],w['radius_px'],name,support);records+=metrics
            result=frozen.estimate_flow(EventBatch(xy[:,0],xy[:,1],times,pol),*d['center'],w['radius_px'])
            counts=dict(Counter(m['status'] for m in metrics));ratio=[m['rms_over_std'] for m in metrics if m['rms_over_std'] is not None]
            native.append(dict(support=support,events=len(xy),same_physical_seeds=len(seed_indices),exit_reasons=counts,
                independent_accepted=sum(m['status']=='accepted' for m in metrics),native_result=result,
                median_neighbors=float(np.median([m['neighbors'] for m in metrics])),median_rms_over_std=float(np.median(ratio)) if ratio else None,
                timestamp_dtype=str(times.dtype),timestamp_span_us=int(np.ptp(times)),seconds_span=float(np.ptp(seconds)),
                note='Native full-support interface uses its own 300 uniform seeds; independent comparison reuses exact physical seed IDs from selected input.'))
            ax=axes[row,col];ax.hist(ratio,bins=np.linspace(0,1,26));ax.axvline(.15,color='red',ls='--');ax.set(title=f'{name}\n{support}',xlabel='plane RMS / local time std',ylabel='Seed count')
        summaries.append(dict(window=name,radius_px=w['radius_px'],center=d['center'].tolist(),native_and_independent=native))
        del e
    table(OUT/'plane_seed_diagnostics.csv',records);write_json(OUT/'plane_sampling_diagnostic.json',dict(status='completed',frozen_flow_sha256=digest(P1/'frozen/optical_flow.py'),
        representative_windows=summaries,seed_rows=len(records),same_thresholds=True,main_324_results_unchanged=True,
        caveat='Diagnostic only, no replacement inputs or threshold changes. Three windows cannot establish a general cause of missing plane support.'))
    fig.suptitle('Same 3 real windows, same 300 physical seeds, unchanged neighborhood/plane gates | red=.15 (20us floor also applies)')
    fig.tight_layout();fig.savefig(OUT/'plane_residual_diagnostics.png',dpi=150);plt.close(fig)
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
