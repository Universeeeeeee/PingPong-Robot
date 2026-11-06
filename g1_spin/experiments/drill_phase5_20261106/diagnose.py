"""Fixed raw patches and a nonparametric single-valued time-map lower bound."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];P4=OUT.parent/'drill_phase4_20261106';FRONT=P4/'real_frontend'
CASES=['run1_01495000','run2_02095000','run5_04095000'];CONTROL='run5_01525000'
NAMES={1:'DSP03F100025_20261105_201058',2:'DSP03F100025_20261105_201123',5:'DSP03F100025_20261105_215030'}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False))
def table(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def dist(values):
    a=np.asarray(values,float);return dict(n=len(a),min=float(a.min()) if len(a) else None,median=float(np.median(a)) if len(a) else None,
        p95=float(np.percentile(a,95)) if len(a) else None,max=float(a.max()) if len(a) else None)

def decomposition(xy,time,polarity=None):
    if not len(time):return dict(events=0,groups=0,status='empty')
    keys=xy if polarity is None else np.column_stack([xy,polarity]);unique,inverse,count=np.unique(keys,axis=0,return_inverse=True,return_counts=True)
    means=np.bincount(inverse,weights=time)/count;centered=time-means[inverse]
    within=float(np.mean(centered*centered));between=float(np.sum(count*(means-time.mean())**2)/len(time));total=float(np.var(time))
    gate=max(20e-6,.15*np.sqrt(total));lower=float(np.sqrt(within))
    return dict(events=len(time),groups=len(unique),single_event_groups=int((count==1).sum()),multiple_event_groups=int((count>1).sum()),
        event_fraction_in_multiple_groups=float(count[count>1].sum()/len(time)),total_variance_s2=total,within_variance_s2=within,between_variance_s2=between,
        variance_identity_error=abs(total-within-between),irreducible_RMS_us=lower*1e6,old_gate_us=float(gate*1e6),lower_over_gate=lower/gate,
        impossible_under_old_gate=bool(lower>gate),maximum_single_value_R2=None if total<=0 else 1-within/total,
        gate_branch='20us_floor' if .15*np.sqrt(total)<20e-6 else '0.15std',
        group_mean_times_s=means.tolist(),group_keys=unique.tolist(),group_counts=count.tolist())

def selected_patches(xy,times):
    mask=times<.006;cells=np.floor(xy[mask]/12).astype(int);unique,count=np.unique(cells,axis=0,return_counts=True)
    candidates=sorted(range(len(unique)),key=lambda i:(-count[i],int(unique[i,0]),int(unique[i,1])));chosen=[]
    for i in candidates:
        if count[i]<20:continue
        cell=unique[i];center=cell*12+6
        if any(np.linalg.norm(center-np.array(r['center']))<16 for r in chosen):continue
        chosen.append(dict(patch=len(chosen),cell=cell.tolist(),center=center.tolist(),box=[int(cell[0]*12),int(cell[1]*12),int(cell[0]*12+12),int(cell[1]*12+12)],selection_events_first6ms=int(count[i])))
        if len(chosen)==6:break
    return chosen

def load(run):
    folder=ROOT/'g1_spin/data/archive_20261105/decoded'/NAMES[run];manifest=json.loads((folder/'manifest.json').read_text())
    events=np.concatenate([np.load(folder/p) for p in manifest['parts']]);return events,int(events['t'][0]),folder,manifest

def per_pixel(xy,t,pol):
    rows=[];intervals=[]
    if not len(t):return rows,intervals
    keys=np.column_stack([xy,pol]);unique,inv=np.unique(keys,axis=0,return_inverse=True)
    for k,key in enumerate(unique):
        times=np.sort(t[inv==k]);gap=np.diff(times);intervals.extend(gap.tolist());clusters=int(1+np.count_nonzero(gap>50e-6))
        rows.append(dict(x=int(key[0]),y=int(key[1]),polarity=int(key[2]),count=len(times),first_us=float(times[0]*1e6),
            median_us=float(np.median(times)*1e6),latest_us=float(times[-1]*1e6),time_std_us=float(np.std(times)*1e6),
            span_us=float(np.ptp(times)*1e6),clusters_gap50us=clusters,interval_us=(gap*1e6).tolist(),
            fraction_intervals_le50us=float(np.mean(gap<=50e-6)) if len(gap) else None))
    return rows,intervals

def metrics(xy,t,pol,outside_counts):
    before=t<.006;pixel,intervals=per_pixel(xy[before],t[before],pol[before]);allpixel,allintervals=per_pixel(xy,t,pol)
    counts,_=np.histogram(t,bins=np.linspace(0,.01,201));corr=None
    if np.std(counts)>0 and np.std(outside_counts)>0:corr=float(np.corrcoef(counts,outside_counts)[0,1])
    overlap=None;nearest=[]
    if len(xy):
        off=np.unique(xy[(pol==0)&before],axis=0);on=np.unique(xy[(pol==1)&before],axis=0)
        if len(off) and len(on):
            nearest=np.r_[cKDTree(off).query(on)[0],cKDTree(on).query(off)[0]].tolist()
            overlap=len(set(map(tuple,off))&set(map(tuple,on)))/len(set(map(tuple,off))|set(map(tuple,on)))
    return dict(events_full10ms=len(t),events_first6ms=int(before.sum()),events_last4ms=int((~before).sum()),
        training_only_decomposition=decomposition(xy[before],t[before],pol[before]),
        training_only_per_polarity={str(p):decomposition(xy[before&(pol==p)],t[before&(pol==p)]) for p in [0,1]},
        training_only_pixel_records=pixel,full10ms_pixel_records=allpixel,training_same_pixel_intervals_us=dist(np.array(intervals)*1e6),
        full10ms_same_pixel_intervals_us=dist(np.array(allintervals)*1e6),
        training_fraction_intervals_le50us=float(np.mean(np.array(intervals)<=50e-6)) if intervals else None,
        training_pixel_polarity_overlap_Jaccard=overlap,training_opposite_polarity_nearest_pixel_distance=dist(nearest),
        full10ms_rate_bin50us=counts.tolist(),outside_rate_correlation_full10ms=corr,
        analysis_rule='Only first6ms enters selection/variance; full10ms and rate correlation are descriptive observation, not training or future validation.')

def plot_patch(path,xy,t,pol,box,pixels):
    if all(path.with_name(path.name+suffix).exists() for suffix in ['_xyt.png','_slices.png','_maps.png','_raster_intervals.png']):return
    # Separate polarities and time units; every raw event is retained in saved NPZ.
    fig=plt.figure(figsize=(13,6));axes=[fig.add_subplot(1,2,i+1,projection='3d') for i in range(2)]
    for p,ax in enumerate(axes):
        keep=pol==p;ax.scatter(xy[keep,0],xy[keep,1],t[keep]*1000,s=4,c=t[keep]*1000,cmap='viridis',vmin=0,vmax=10)
        ax.set(xlim=(box[0],box[2]),ylim=(box[1],box[3]),zlim=(0,10),xlabel='x pixel',ylabel='y pixel',zlabel='time ms',title=f'polarity{p} all raw events')
    fig.tight_layout();fig.savefig(path.with_name(path.name+'_xyt.png'),dpi=130);plt.close(fig)
    fig,axes=plt.subplots(4,5,figsize=(13,10))
    for k,ax in enumerate(axes.flat):
        keep=(t>=k*.0005)&(t<(k+1)*.0005)
        for p,color in [(0,'#249bc5'),(1,'#e39424')]:
            m=keep&(pol==p);ax.scatter(xy[m,0],xy[m,1],s=9,c=color)
        ax.set(xlim=(box[0]-.5,box[2]-.5),ylim=(box[3]-.5,box[1]-.5),aspect='equal',title=f'{k*.5:g}–{(k+1)*.5:g}ms n{int(keep.sum())}');ax.tick_params(labelsize=7)
    fig.suptitle('Raw event0.5ms slices; t<6ms selection/measurement, later slices observation only');fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(path.with_name(path.name+'_slices.png'),dpi=130);plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(14,7))
    for p in [0,1]:
        for col,key in enumerate(['count','first_us','median_us','latest_us']):
            arr=np.full((12,12),np.nan)
            for r in pixels:
                if r['polarity']==p:arr[r['y']-box[1],r['x']-box[0]]=r[key] if key=='count' else r[key]/1000
            im=axes[p,col].imshow(arr,extent=[box[0]-.5,box[2]-.5,box[3]-.5,box[1]-.5],vmin=0,vmax=None if key=='count' else 10,cmap='viridis')
            axes[p,col].set_title(f'pol{p} {key if key=="count" else key.replace("_us"," ms")}');fig.colorbar(im,ax=axes[p,col],fraction=.05)
    fig.suptitle('Full10ms per-pixel maps; missing pixels remain NaN');fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(path.with_name(path.name+'_maps.png'),dpi=130);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,8))
    for p in [0,1]:
        keep=pol==p;rows=(xy[keep,1]-box[1])*12+xy[keep,0]-box[0]
        axes[p,0].scatter(t[keep]*1000,rows,s=5);axes[p,0].axvline(6,color='red',ls='--');axes[p,0].set(xlim=(0,10),xlabel='time ms',ylabel='pixel local index',title=f'pol{p} event raster')
        intervals=[]
        for r in pixels:
            if r['polarity']==p:intervals+=r['interval_us']
        if intervals:axes[p,1].hist(intervals,bins=np.r_[np.linspace(0,200,21),np.linspace(400,10000,25)])
        axes[p,1].axvline(50,color='red',ls='--');axes[p,1].set(xlabel='same pixel interval us',ylabel='count',title=f'pol{p} intervals, no inlier selection')
    fig.tight_layout();fig.savefig(path.with_name(path.name+'_raster_intervals.png'),dpi=130);plt.close(fig)

def old_neighborhood_bound(manifest):
    rows=[];summaries=[]
    for case in CASES:
        r=next(r for r in manifest if r['case']==case);data=np.load(FRONT/r['input_file']);lookup={int(i):k for k,i in enumerate(data['ids'])}
        detail=json.loads(gzip.decompress((FRONT/'development'/f'{case}_B1_F0.json.gz').read_bytes()));local=[]
        for seed in detail['extraction_metadata']['seed_diagnostics']:
            ids=seed['fit_physical_ids']
            if len(ids)<10:continue
            ix=np.array([lookup[i] for i in ids]);m=decomposition(data['xy'][ix],data['times'][ix],data['polarity'][ix])
            row=dict(case=case,seed_physical_id=seed['seed_physical_id'],neighbors=len(ids),**{k:v for k,v in m.items() if not k.startswith('group_')})
            rows.append(row);local.append(row)
        summaries.append(dict(case=case,eligible_neighbors=len(local),single_value_lower_bound_exceeds_old_gate=sum(r['impossible_under_old_gate'] for r in local),
            lower_over_gate=dist([r['lower_over_gate'] for r in local]),interpretation='Any single-valued pixel/polarity timestamp function fails these neighborhoods, not merely the linear optimizer. F1 is outside this argument.'))
    table(OUT/'old_F0_neighbor_lower_bounds.csv',rows);write(OUT/'old_F0_neighbor_lower_bound_summary.json',summaries)
    return summaries

def main():
    if (OUT/'patch_metrics.json').exists():raise FileExistsError('Preserve fixed raw structure inspection')
    (OUT/'patches').mkdir(exist_ok=True);manifest=json.loads((FRONT/'input_manifest.json').read_text());protected={}
    paths=[FRONT/'input_manifest.json',FRONT/'frontends.py',FRONT/'pipeline.py',FRONT/'prediction.py',FRONT/'parameter_freeze.json',FRONT/'selection_development.json']
    for p in paths:protected[str(p.relative_to(ROOT))]=digest(p)
    for case in CASES+[CONTROL]:
        r=next(r for r in manifest if r['case']==case);p=FRONT/r['input_file'];protected[str(p.relative_to(ROOT))]=digest(p)
        if case in CASES:
            p=FRONT/'development'/f'{case}_B1_F0.json.gz';protected[str(p.relative_to(ROOT))]=digest(p)
    selections=[];allmetrics=[];outside=[]
    for run in [1,2,5]:
        events,origin,folder,source_manifest=load(run)
        for name in source_manifest['parts']:
            p=folder/name;protected[str(p.relative_to(ROOT))]=digest(p)
        protected[str((folder/'manifest.json').relative_to(ROOT))]=digest(folder/'manifest.json')
        for case in [c for c in CASES+[CONTROL] if c.startswith(f'run{run}_')]:
            r=next(r for r in manifest if r['case']==case);data=np.load(FRONT/r['input_file']);chosen=selected_patches(data['xy'],data['times'])
            begin=r['start_us'];lo,hi=np.searchsorted(events['t'],[origin+begin,origin+begin+10000]);raw=events[lo:hi];xy=np.column_stack([raw['x'],raw['y']]).astype(int);time=(raw['t'].astype(np.int64)-origin-begin)*1e-6
            outer=np.linalg.norm(xy-np.array(r['center']),axis=1)>1.5*r['radius_px'];ocounts,_=np.histogram(time[outer],bins=np.linspace(0,.01,201))
            sphere_count,_=np.histogram(time[np.linalg.norm(xy-np.array(r['center']),axis=1)<=r['radius_px']],bins=np.linspace(0,.01,201))
            outside.append(dict(case=case,raw_window_events=len(raw),outside_1p5r_events=int(outer.sum()),outside_rate_bin50us=ocounts.tolist(),
                within_r_rate_bin50us=sphere_count.tolist(),sphere_outside_correlation=float(np.corrcoef(ocounts,sphere_count)[0,1]) if np.std(ocounts)>0 and np.std(sphere_count)>0 else None,
                caveat='Outside includes shaft, hand and background motion; event-rate correlation alone cannot identify illumination or synchronization bursts.'))
            selection=dict(case=case,condition=r['condition'],center=r['center'],radius_px=r['radius_px'],selection_support='stage4 complete common_ids first6ms only',
                selected_count=len(chosen),maximum_allowed=6,insufficient_patch_slots=6-len(chosen),patches=chosen,common_events=len(data['ids']),
                first6ms_common_events=int((data['times']<.006).sum()),source_file=r['input_file'],source_sha256=digest(FRONT/r['input_file']))
            selections.append(selection)
            for patch in chosen:
                box=patch['box'];keep=(xy[:,0]>=box[0])&(xy[:,0]<box[2])&(xy[:,1]>=box[1])&(xy[:,1]<box[3]);ids=np.flatnonzero(keep)+lo
                px=xy[keep];tt=time[keep];pol=raw['p'][keep];prefix=f'{case}_patch{patch["patch"]}'
                np.savez_compressed(OUT/'patches'/f'{prefix}.npz',ids=ids,xy=px,time_s=tt,polarity=pol,in_stage4_common=np.isin(ids,data['ids']),origin_us=origin,start_us=begin)
                result=metrics(px,tt,pol,ocounts);result.update(case=case,patch=patch,events_file=f'patches/{prefix}.npz')
                allmetrics.append(result);plot_patch(OUT/'patches'/prefix,px,tt,pol,box,result['full10ms_pixel_records'])
            print('raw structure',case,'selected',len(chosen),'events',[(m['patch']['patch'],m['events_full10ms']) for m in allmetrics if m['case']==case],flush=True)
        del events
    write(OUT/'patch_selection.json',selections);write(OUT/'patch_metrics.json',allmetrics);write(OUT/'outside_activity.json',outside);write(OUT/'source_hashes.json',protected)
    old=old_neighborhood_bound(manifest);print('old F0 lower bounds',[(r['case'],r['single_value_lower_bound_exceeds_old_gate'],r['eligible_neighbors'],r['lower_over_gate']['median']) for r in old],flush=True)

if __name__=='__main__':main()
