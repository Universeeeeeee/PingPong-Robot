"""B2 raw temporal evidence cards; all morphology descriptors remain nonsemantic."""
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
from scipy.ndimage import label

OUT=Path(__file__).resolve().parent
B0=OUT.parent/'marker_b0_20261101'
sys.path.insert(0,str(B0))
spec=importlib.util.spec_from_file_location('b0_reference',B0/'run.py');b0=importlib.util.module_from_spec(spec);spec.loader.exec_module(b0)
sys.path.pop(0)
REVIEW_FRAMES=[0,5,10,15,19]


def digest(path):return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in path.iterdir() if p.is_file()}


def read_rows():
    rows=list(csv.DictReader((B0/'observations.csv').open()))
    for r in rows:
        for k,v in r.items():
            if v in ['True','False']:r[k]=v=='True'
            elif not v:r[k]=np.nan
            else:
                try:r[k]=float(v)
                except ValueError:pass
    return rows


def write_csv(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def morphology(events,center):
    local=b0.local_events(events,center);pts=b0.unique_xy(local)
    occupancy=np.zeros((27,27),dtype=bool)
    if len(pts):occupancy[(pts[:,1]-center[1]+13).astype(int),(pts[:,0]-center[0]+13).astype(int)]=True
    labels,n=label(occupancy,np.ones((3,3),dtype=int));sizes=np.bincount(labels.ravel())[1:]
    border=np.zeros((27,27),dtype=bool);border[:2]=True;border[-2:]=True;border[:,:2]=True;border[:,-2:]=True
    result=dict(raw_patch_events=len(local),distinct_pixels=len(pts),components_ge3_pixels=int(np.sum(sizes>=3)),
                largest_component_pixel_fraction=float(np.max(sizes)/len(pts)) if len(pts) else np.nan,
                boundary_2px_pixel_fraction=float(np.sum(occupancy&border)/len(pts)) if len(pts) else np.nan)
    centroids=[]
    for pol,name in [(True,'on'),(False,'off')]:
        e=local[(local['p']>0)==pol];p=b0.unique_xy(e)
        result[name+'_events']=len(e);result[name+'_pixels']=len(p)
        centroids.append(np.mean(p,axis=0) if len(p) else np.full(2,np.nan))
    result['on_off_centroid_distance']=float(np.linalg.norm(centroids[0]-centroids[1]))
    return result


def main():
    old_hashes=digest(B0);rows=read_rows();lookup={(r['segment'],int(r['duration_ms']),r['feature'],int(r['frame'])):r for r in rows}
    events,manifest=b0.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);cache={};window_checks=[]
    for segment,(start,end) in b0.raw.SEGMENTS.items():
        for duration in [2,5]:
            for frame in range(100//duration):
                lo=t0+round(start*1e6)+frame*duration*1000;hi=lo+duration*1000;a,b=np.searchsorted(times,[lo,hi])
                selected=events[a:b];direct=events[(times>=lo)&(times<hi)];assert np.array_equal(selected,direct)
                cache[segment,duration,frame]=(selected,(lo-t0)/1e6,(hi-t0)/1e6)
                window_checks.append(dict(segment=segment,duration_ms=duration,frame=frame,start_absolute_us=lo,end_absolute_us=hi,
                                          events_full_sensor=len(selected),independent_mask_matches=True))
    morphology_rows=[]
    for r in rows:
        identity={k:r[k] for k in ['segment','duration_ms','feature','frame','start_s','end_s','midpoint_s','status','position_valid','anchor_x','anchor_y']}
        selected=cache[r['segment'],int(r['duration_ms']),int(r['frame'])][0]
        diag=morphology(selected,[r['anchor_x'],r['anchor_y']]) if r['position_valid'] else {}
        for k in ['circle_valid','ellipse_valid','circle_reason','ellipse_reason','circle_occupied_30deg_bins','ellipse_occupied_30deg_bins',
                  'circle_max_empty_angle_deg','ellipse_max_empty_angle_deg','circle_median_abs_residual','ellipse_median_abs_residual']:
            diag[k]=r.get(k,np.nan)
        morphology_rows.append(dict(**identity,**diag))
    write_csv('morphology.csv',morphology_rows);write_csv('window_checks.csv',window_checks)
    card_records=[]
    for segment in b0.raw.SEGMENTS:
        for feature,seed in b0.SEEDS[segment].items():
            fig,axes=plt.subplots(2,5,figsize=(17,8))
            for j,frame in enumerate(REVIEW_FRAMES):
                r5=lookup[segment,5,feature,frame]
                center=np.array([r5['anchor_x'],r5['anchor_y']]) if r5['position_valid'] else np.array(seed,dtype=float)
                mid=r5['midpoint_s'];frames2=[r for r in rows if r['segment']==segment and r['duration_ms']==2 and r['feature']==feature]
                r2=min(frames2,key=lambda r:abs(r['midpoint_s']-mid));frame2=int(r2['frame'])
                for i,(duration,k,r) in enumerate([(5,frame,r5),(2,frame2,r2)]):
                    ax=axes[i,j];selected,lo,hi=cache[segment,duration,k]
                    visible=selected[(abs(selected['x'].astype(float)-center[0])<=20)&(abs(selected['y'].astype(float)-center[1])<=20)]
                    b0.raw.draw_raw(ax,visible,4)
                    ax.set(xlim=(center[0]-20.5,center[0]+20.5),ylim=(center[1]+20.5,center[1]-20.5),aspect='equal')
                    ax.add_patch(plt.Rectangle(center-13,26,26,fill=False,color='#b03090' if i==0 else '#777777',lw=.8,linestyle='-' if r5['position_valid'] else ':'))
                    if r['position_valid']:ax.plot(r['anchor_x'],r['anchor_y'],'x',color='#b03090' if i==0 else '#008b77',ms=5)
                    association='valid anchor' if r['position_valid'] else 'NO matched anchor'
                    crop_name='5ms valid-anchor crop' if r5['position_valid'] else 'FIXED SEED crop; not a track'
                    circle='Y' if r.get('circle_valid') is True else 'N';ellipse='Y' if r.get('ellipse_valid') is True else 'N'
                    ax.set_title(f'{duration}ms f{k:02d} [{lo:.3f},{hi:.3f})s\n{association}; circle/ellipse={circle}/{ellipse}',fontsize=8)
                    ax.tick_params(labelsize=6)
                    if i==1:ax.set_xlabel(f'{crop_name}\n2ms−5ms mid={1000*(r2["midpoint_s"]-mid):+.1f}ms',fontsize=7)
                    card_records.append(dict(segment=segment,feature=feature,review_5ms_frame=frame,raw_duration_ms=duration,raw_frame=k,
                                             start_s=lo,end_s=hi,raw_midpoint_s=(lo+hi)/2,crop_center_x=center[0],crop_center_y=center[1],
                                             crop_source='valid_5ms_anchor' if r5['position_valid'] else 'fixed_manual_seed_NO_association',
                                             independent_duration_anchor_valid=r['position_valid'],independent_duration_status=r['status'],
                                             midpoint_difference_2_minus_5_ms=1000*(r2['midpoint_s']-mid),
                                             human_semantics='UNREVIEWED; no automatic logo/Marker label'))
            fig.suptitle(f'{segment} {feature}: raw temporal morphology, ON orange / OFF blue; same-pixel dual polarity split\n'
                         '27px box is matching support; 41px context reveals truncation/neighbors. Circle failure does NOT imply logo.',fontsize=11)
            fig.subplots_adjust(left=.04,right=.995,bottom=.1,top=.87,wspace=.22,hspace=.35)
            fig.savefig(OUT/f'card_{segment}_{feature}.png',dpi=180);plt.close(fig)
    write_csv('card_windows.csv',card_records)
    # Full-context midpoints prevent local cropping from silently removing neighboring explanations.
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,(segment,(start,end)) in zip(axes,b0.raw.SEGMENTS.items()):
        selected,lo,hi=cache[segment,5,10];b0.raw.draw_raw(ax,b0.raw.crop(selected),3)
        for feature,seed in b0.SEEDS[segment].items():
            ax.plot(*seed,'x',color='#b03090');ax.text(seed[0]+2,seed[1]-2,feature,fontsize=9,color='#b03090')
            ax.add_patch(plt.Rectangle(np.array(seed)-13,26,26,fill=False,color='#b03090',lw=.8))
        ax.set_title(f'{segment} [{lo:.3f},{hi:.3f})s\nAll F labels are manual patch names',fontsize=10)
    fig.suptitle('Reference full ROI: raw events, no semantic labels or model circles',fontsize=12)
    fig.subplots_adjust(left=.04,right=.995,bottom=.1,top=.82,wspace=.18);fig.savefig(OUT/'reference_context.png',dpi=180);plt.close(fig)
    assert digest(B0)==old_hashes
    meta=dict(source=str(b0.raw.DATA),source_t0_us=t0,segments=b0.raw.SEGMENTS,seeds=b0.SEEDS,review_5ms_frames=REVIEW_FRAMES,
              raw_window_count=len(window_checks),morphology_rows=len(morphology_rows),card_panels=len(card_records),
              morphology_is_descriptive_not_semantic_classifier=True,invalid_positions_not_filled=True,
              previous_B0_sha256=old_hashes,manifest_sha256=hashlib.sha256((b0.raw.DATA/'manifest.json').read_bytes()).hexdigest(),
              input_parts_sha256={f:hashlib.sha256((b0.raw.DATA/f).read_bytes()).hexdigest() for f in manifest['parts']},
              review_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              protocol_sha256=hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest())
    (OUT/'review_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Created nine raw evolution cards (90 panels), reference context, 630 morphology rows and 210 independent window checks.')


if __name__=='__main__':main()
