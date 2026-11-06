"""B1 input interventions only; B0 local appearance templates and all gates unchanged."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

OUT=Path(__file__).resolve().parent
B0=OUT.parent/'marker_b0_20261101'
sys.path.insert(0,str(B0))
spec=importlib.util.spec_from_file_location('b0_frozen_run',B0/'run.py')
b0=importlib.util.module_from_spec(spec);spec.loader.exec_module(b0)
sys.path.pop(0)
ACTIVE=[13,14,15]
VARIANTS=['occlusion','clone_plus8','outside_context_delete']


def digest(path):
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(path.iterdir()) if p.is_file()}


def write_csv(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def read_b0():
    result={}
    for r in csv.DictReader((B0/'observations.csv').open()):
        if r['duration_ms']!='5':continue
        for k in ['frame','anchor_x','anchor_y']:
            r[k]=float(r[k]) if r[k] else np.nan
        result[r['segment'],r['feature'],int(r['frame'])]=r
    return result


def forward(frames,images,seed,identity):
    """Reference frame10 to frame19 only. Reverse check consumes only prior/current frames."""
    center=np.array(seed,dtype=int);previous=10;active=True;rows=[]
    for k in range(10,20):
        row=dict(**identity,frame=k,start_s=frames[k]['start_s'],end_s=frames[k]['end_s'],
                 midpoint_s=(frames[k]['start_s']+frames[k]['end_s'])/2,anchor_x=np.nan,anchor_y=np.nan,
                 position_valid=False,max_frame_consumed=k)
        if k==10:
            enough=len(b0.unique_xy(b0.local_events(frames[k]['events'],center)))>=12
            row.update(status='manual_seed' if enough else 'manual_seed_insufficient_pixels',anchor_x=int(center[0]),anchor_y=int(center[1]),position_valid=enough)
            active=enough;rows.append(row);continue
        if not active:
            row.update(status='not_attempted_after_loss');rows.append(row);continue
        diag=b0.match(images[previous],images[k],center,center);reasons=[]
        if not diag['match_available']:
            reasons=['image_boundary'];proposal=center
        else:
            proposal=np.array([diag['proposed_x'],diag['proposed_y']])
            source_count=len(b0.unique_xy(b0.local_events(frames[previous]['events'],center)))
            target_count=len(b0.unique_xy(b0.local_events(frames[k]['events'],proposal)))
            reverse=b0.match(images[k],images[previous],proposal,center)
            fb=float(np.hypot(reverse['proposed_x']-center[0],reverse['proposed_y']-center[1])) if reverse['match_available'] else np.inf
            diag.update(source_distinct_pixels=source_count,target_distinct_pixels=target_count,forward_backward_px=fb,
                        reverse_score=reverse.get('score',np.nan),source_frame=previous,direction=1)
            if min(source_count,target_count)<12:reasons.append('insufficient_pixels')
            if diag['score']<.35:reasons.append('score')
            if diag['score_margin']<.025:reasons.append('ambiguous_peak')
            if diag['search_boundary']:reasons.append('search_boundary')
            if fb>2:reasons.append('forward_backward')
        row.update(diag)
        if reasons:
            active=False;row['status']='rejected_'+'+'.join(reasons)
        else:
            row.update(status='accepted',position_valid=True,anchor_x=int(proposal[0]),anchor_y=int(proposal[1]))
            center=proposal;previous=k
        rows.append(row)
    return rows


def keep_square(events,center,half):
    return (abs(events['x'].astype(float)-center[0])<=half)&(abs(events['y'].astype(float)-center[1])<=half)


def occupancy(events):
    out=np.zeros((2,b0.H,b0.W),dtype=np.uint8)
    e=events[(events['x']>=b0.X0)&(events['x']<b0.X0+b0.W)&(events['y']>=b0.Y0)&(events['y']<b0.Y0+b0.H)]
    out[(e['p']>0).astype(int),e['y'].astype(int)-b0.Y0,e['x'].astype(int)-b0.X0]=1
    return out


def modify(frames,segment,variant,baseline):
    output=[];stats=[];archive={}
    for k,frame in enumerate(frames):
        original=frame['events'];modified=original.copy();center=np.array([baseline[segment,'F1',k]['anchor_x'],baseline[segment,'F1',k]['anchor_y']])
        added=original[:0];removed=original[:0];source=original[:0];protected=np.zeros(len(original),dtype=bool)
        if variant in ['occlusion','clone_plus8'] and k in ACTIVE:
            chosen=keep_square(original,center,13)
            if variant=='occlusion':
                removed=original[chosen].copy();modified=original[~chosen].copy()
                assert np.array_equal(modified,original[~chosen])
            else:
                source=original[chosen].copy();added=source.copy();added['x']+=8
                assert np.array_equal(added['t'],source['t']) and np.array_equal(added['p'],source['p'])
                for key in original.dtype.names:
                    if key!='x':assert np.array_equal(added[key],source[key])
                assert np.array_equal(added['x'].astype(int)-source['x'].astype(int),np.full(len(source),8))
                modified=np.concatenate([original.copy(),added]);modified=modified[np.argsort(modified['t'],kind='stable')]
        elif variant=='outside_context_delete' and k>=10:
            protected=keep_square(original,center,30)
            if k>10:
                prior=baseline[segment,'F1',k-1];protected|=keep_square(original,[prior['anchor_x'],prior['anchor_y']],30)
            modified=original[protected].copy();removed=original[~protected].copy()
            assert np.array_equal(modified,original[protected])
        output.append(dict(events=modified,start_s=frame['start_s'],end_s=frame['end_s']))
        if k>=10:
            before,after=occupancy(original),occupancy(modified)
            stats.append(dict(segment=segment,variant=variant,frame=k,intervention_active=bool(k in ACTIVE) if variant!='outside_context_delete' else True,
                              baseline_anchor_x=center[0],baseline_anchor_y=center[1],original_events=len(original),modified_events=len(modified),
                              added_events=len(added),removed_events=len(removed),occupancy_changed_polarity_pixels=int(np.count_nonzero(before!=after)),
                              occupancy_added_polarity_pixels=int(np.count_nonzero((after>0)&(before==0))),
                              occupancy_removed_polarity_pixels=int(np.count_nonzero((before>0)&(after==0))),
                              source_patch_events=len(source),source_patch_distinct_pixels=len(b0.unique_xy(source)),
                              original_payload_preserved_outside_mask=True))
            archive[f'frame_{k:02d}']=modified;archive[f'added_{k:02d}']=added;archive[f'removed_{k:02d}']=removed
            if len(source):archive[f'clone_source_{k:02d}']=source
    np.savez_compressed(OUT/f'{segment}_{variant}_events.npz',**archive)
    return output,stats


def add_diagnostics(rows,baseline,variant):
    for r in rows:
        b=baseline[r['segment'],'F1',r['frame']];base=np.array([b['anchor_x'],b['anchor_y']],dtype=float)
        r.update(baseline_anchor_x=base[0],baseline_anchor_y=base[1],copy_anchor_x=base[0]+8 if r['frame'] in ACTIVE and variant=='clone_plus8' else np.nan,
                 copy_anchor_y=base[1] if r['frame'] in ACTIVE and variant=='clone_plus8' else np.nan)
        if r['position_valid']:
            p=np.array([r['anchor_x'],r['anchor_y']]);distance=float(np.linalg.norm(p-base));r['accepted_distance_from_baseline']=distance
            if variant=='clone_plus8' and r['frame'] in ACTIVE:
                dc=float(np.linalg.norm(p-(base+[8,0])));r['accepted_distance_from_copy_anchor']=dc
                r['copy_anchor_attraction']=bool(dc<=3 and distance>=5 and dc+3<=distance)
        if 'proposed_x' in r:
            p=np.array([r['proposed_x'],r['proposed_y']]);r['proposal_distance_from_baseline']=float(np.linalg.norm(p-base))
            if variant=='clone_plus8' and r['frame'] in ACTIVE:r['proposal_distance_from_copy_anchor']=float(np.linalg.norm(p-(base+[8,0])))


def independent_blurred(events):
    """Independent bincount accumulation, distinct occupancy, explicit three-pixel finite kernel."""
    use=events[(events['x']>=b0.X0)&(events['x']<b0.X0+b0.W)&(events['y']>=b0.Y0)&(events['y']<b0.Y0+b0.H)]
    images=[]
    for pol in [False,True]:
        e=use[(use['p']>0)==pol]
        counts=np.bincount((e['y'].astype(int)-b0.Y0)*b0.W+e['x'].astype(int)-b0.X0,minlength=b0.H*b0.W).reshape(b0.H,b0.W)
        images.append(gaussian_filter((counts>0).astype(float),.8,radius=3))
    return np.array(images)


def context_validation(original,modified,baseline_rows,segment):
    checks=[];original_images=[independent_blurred(f['events']) for f in original];modified_images=[independent_blurred(f['events']) for f in modified]
    for r in baseline_rows:
        if r['status']!='accepted':continue
        k=r['frame'];previous=int(r['source_frame']);source=next(x for x in baseline_rows if x['frame']==previous)
        anchor=[source['anchor_x'],source['anchor_y']];proposal=[r['proposed_x'],r['proposed_y']]
        # Forward template/search and reverse template/search: exactly what b0.match reads.
        accesses=[(previous,anchor,13,'forward_template'),(k,anchor,25,'forward_search'),
                  (k,proposal,13,'reverse_template'),(previous,anchor,25,'reverse_search')]
        for frame,center,half,name in accesses:
            a=b0.patch(original_images[frame],center,half);b=b0.patch(modified_images[frame],center,half)
            assert np.array_equal(a,b)
            checks.append(dict(segment=segment,frame=frame,matching_target_frame=k,access=name,halfwidth=half,
                               blurred_arrays_bitwise_equal=True,sha256=hashlib.sha256(a.tobytes()).hexdigest()))
    return checks


def main():
    previous_hashes=digest(B0);oldrows=read_b0();events,manifest=b0.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0])
    baselines=[];observations=[];interventions=[];contexts=[];cache={};isolation=[]
    for segment,(start,end) in b0.raw.SEGMENTS.items():
        frames=[]
        for k in range(20):
            lo=t0+round(start*1e6)+k*5000;hi=lo+5000;a,b=np.searchsorted(times,[lo,hi]);selected=events[a:b]
            assert np.array_equal(selected,events[(times>=lo)&(times<hi)])
            frames.append(dict(events=selected,start_s=(lo-t0)/1e6,end_s=(hi-t0)/1e6))
        images=[b0.image(f['events']) for f in frames];cache[segment,'baseline']=frames
        segment_baselines={}
        for feature,seed in b0.SEEDS[segment].items():
            rr=forward(frames,images,seed,dict(segment=segment,feature=feature,variant='baseline'))
            repeated=forward(frames,images,seed,dict(segment=segment,feature=feature,variant='baseline'))
            assert json.dumps(rr,sort_keys=True)==json.dumps(repeated,sort_keys=True)
            for r in rr:
                prior=oldrows[segment,feature,r['frame']]
                assert r['status']==prior['status']
                for key in ['anchor_x','anchor_y','proposed_x','proposed_y','score','second_score','score_margin','forward_backward_px','reverse_score','source_distinct_pixels','target_distinct_pixels']:
                    if key in r:
                        x,y=float(r[key]),float(prior[key]);assert (np.isnan(x) and np.isnan(y)) or x==y,(segment,feature,r['frame'],key,x,y)
            isolation.append(dict(segment=segment,feature=feature,identical_when_called_independently=True,n_rows=len(rr)))
            baselines.extend(rr);segment_baselines[feature]=rr
        for variant in VARIANTS:
            altered,ss=modify(frames,segment,variant,oldrows);cache[segment,variant]=altered;interventions.extend(ss)
            altered_images=[b0.image(f['events']) for f in altered]
            rr=forward(altered,altered_images,b0.SEEDS[segment]['F1'],dict(segment=segment,feature='F1',variant=variant))
            add_diagnostics(rr,oldrows,variant);observations.extend(rr)
            if variant=='outside_context_delete':
                contexts.extend(context_validation(frames,altered,segment_baselines['F1'],segment))
                for a,b in zip(segment_baselines['F1'],rr):
                    for key in ['status','anchor_x','anchor_y','score','second_score','score_margin','forward_backward_px']:
                        if key in a:
                            x,y=a[key],b[key];assert x==y or (isinstance(x,float) and np.isnan(x) and np.isnan(y))
            print(segment,variant,[(r['frame'],r['status'],r.get('accepted_distance_from_baseline')) for r in rr],flush=True)
    write_csv('baseline.csv',baselines);write_csv('observations.csv',observations);write_csv('interventions.csv',interventions)
    write_csv('context_checks.csv',contexts);write_csv('independent_feature_checks.csv',isolation)
    assert previous_hashes==digest(B0)
    metadata=dict(source=str(b0.raw.DATA),source_t0_us=t0,source_events=len(events),source_dtype=str(events.dtype),
                  manifest_sha256=hashlib.sha256((b0.raw.DATA/'manifest.json').read_bytes()).hexdigest(),
                  input_parts_sha256={f:hashlib.sha256((b0.raw.DATA/f).read_bytes()).hexdigest() for f in manifest['parts']},
                  baseline_source_files_sha256=previous_hashes,seeds=b0.SEEDS,segments=b0.raw.SEGMENTS,
                  normal_forward_frames=list(range(10,20)),active_pressure_frames=ACTIVE,clone_offset_xy=[8,0],
                  removal_patch_halfwidth=13,outside_protected_halfwidth=30,blur_kernel_radius_px=3,
                  matching_image_roi_inclusive=[450,670,120,280],gate_source=str(B0/'protocol.md'),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  protocol_sha256=hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest(),
                  baseline_rows=len(baselines),intervention_observation_rows=len(observations),
                  intervention_archives_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*_events.npz')},
                  caveats=['No human logo labels or physical Marker IDs','B0 already uses no 3D layout, but DOES use local appearance templates',
                           'Oracle baseline paths construct synthetic inputs only, not normal initialization/control','No real layout change recordings',
                           'No fixed-template or subpixel method changes in this ablation'])
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    plot(cache,baselines+observations)


def plot(cache,rows):
    for variant in ['baseline','occlusion','clone_plus8']:
        fig,axes=plt.subplots(3,4,figsize=(16,11))
        for i,segment in enumerate(b0.raw.SEGMENTS):
            for j,k in enumerate([12,13,15,16]):
                ax=axes[i,j];frames=cache[segment,variant];b0.raw.draw_raw(ax,b0.raw.crop(frames[k]['events']),2)
                r=next(r for r in rows if r['segment']==segment and r['feature']=='F1' and r['variant']==variant and r['frame']==k)
                base=next(r for r in rows if r['segment']==segment and r['feature']=='F1' and r['variant']=='baseline' and r['frame']==k)
                p=np.array([base['anchor_x'],base['anchor_y']])
                ax.plot(*p,'+',color='#444444',ms=5,label='baseline oracle anchor')
                if variant=='occlusion' and k in ACTIVE:
                    ax.add_patch(plt.Rectangle(p-13,26,26,fill=False,ls=':',color='#cc4444',lw=1))
                if variant=='clone_plus8' and k in ACTIVE:
                    ax.plot(*(p+[8,0]),'+',color='#00aa99',ms=5)
                    ax.add_patch(plt.Rectangle(p+[8,0]-13,26,26,fill=False,ls=':',color='#00aa99',lw=1))
                if r['position_valid']:
                    q=np.array([r['anchor_x'],r['anchor_y']]);ax.plot(*q,'x',color='#b03090',ms=6)
                    ax.add_patch(plt.Rectangle(q-13,26,26,fill=False,color='#b03090',lw=.7))
                elif 'proposed_x' in r:ax.plot(r['proposed_x'],r['proposed_y'],'x',color='#dd2222',ms=6)
                ax.set_title(f'{segment} f{k} [{frames[k]["start_s"]:.3f},{frames[k]["end_s"]:.3f})s\n{r["status"]}',fontsize=8)
                ax.tick_params(labelsize=6)
                ax.set_xlabel(('pressure ACTIVE' if k in ACTIVE and variant!='baseline' else 'original events')+'; no recovery after rejection',fontsize=7)
        fig.suptitle(f'{variant}: raw events + F1 patch anchor; gray + is offline baseline reference\n'
                     'Magenta x = accepted; red x = rejected proposal; no physical identities or logo labels',fontsize=12)
        fig.subplots_adjust(left=.04,right=.995,bottom=.055,top=.90,wspace=.18,hspace=.5)
        fig.savefig(OUT/f'{variant}_review.png',dpi=180);plt.close(fig)


if __name__=='__main__':main()
