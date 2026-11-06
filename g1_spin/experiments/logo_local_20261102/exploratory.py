"""Post-survey candidate: raw previews first; unmodified B0 forward matching only when requested."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import review

OUT=review.OUT
CONFIG=json.loads((OUT/'exploratory_protocol.json').read_text())
SEEDS=CONFIG['seeds']
START,END=CONFIG['segment_s']


def windows():
    events,manifest=review.b0.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);frames={};checks=[]
    for duration in [2,5]:
        frames[duration]=[]
        for k in range(100//duration):
            lo=t0+round(START*1e6)+k*duration*1000;hi=lo+duration*1000;a,b=np.searchsorted(times,[lo,hi])
            selected=events[a:b];direct=events[(times>=lo)&(times<hi)];assert np.array_equal(selected,direct)
            frames[duration].append(dict(events=selected,start_s=(lo-t0)/1e6,end_s=(hi-t0)/1e6))
            checks.append(dict(duration_ms=duration,frame=k,start_s=(lo-t0)/1e6,end_s=(hi-t0)/1e6,events_sensor=len(selected),independent_absolute_mask_matches=True))
    return frames,checks


def preview(frames,checks):
    for duration in [5,2]:
        n=len(frames[duration]);cols=5 if duration==5 else 10;fig,axes=plt.subplots(n//cols,cols,figsize=(20,12))
        for k,ax in enumerate(axes.flat):
            f=frames[duration][k];review.b0.raw.draw_raw(ax,review.b0.raw.crop(f['events']),1.8 if duration==5 else 1.2)
            ax.set_title(f'f{k:02d} [{f["start_s"]:.3f},{f["end_s"]:.3f})s',fontsize=8 if duration==5 else 6)
            ax.tick_params(labelsize=6 if duration==5 else 4)
        fig.suptitle(f'Post-survey exploratory 11.70–11.80s: ALL consecutive non-overlapping {duration}ms raw windows\n'
                     'No automatic logo labels; chosen after fixed full-recording survey, not an independent test',fontsize=12)
        fig.subplots_adjust(left=.03,right=.995,bottom=.04,top=.9,hspace=.48,wspace=.2)
        fig.savefig(OUT/f'exploratory_raw_{duration}ms.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(11,11))
    for i,(feature,center) in enumerate(SEEDS.items()):
        center=np.array(center)
        for j,k in enumerate([9,10,11]):
            ax=axes[i,j];f=frames[5][k];e=f['events'];local=e[(abs(e['x'].astype(float)-center[0])<=21)&(abs(e['y'].astype(float)-center[1])<=21)]
            review.b0.raw.draw_raw(ax,local,4);ax.set(xlim=(center[0]-21.5,center[0]+21.5),ylim=(center[1]+21.5,center[1]-21.5),aspect='equal')
            ax.plot(*center,'x',color='#b03090');ax.add_patch(plt.Rectangle(center-13,26,26,fill=False,color='#b03090',lw=.8))
            ax.set_title(f'{feature} fixed seed {tuple(center)}\n[{f["start_s"]:.3f},{f["end_s"]:.3f})s',fontsize=9)
    fig.suptitle('Post-survey seeds BEFORE tracking: raw ON/OFF at fixed coordinates, 27px patch boxes\n'
                 'L = irregular-logo candidate only; A1/A2 = round/partial-arc candidates; no verified physical IDs',fontsize=11)
    fig.subplots_adjust(left=.06,right=.98,bottom=.04,top=.9,hspace=.32,wspace=.22)
    fig.savefig(OUT/'exploratory_seed_review.png',dpi=180);plt.close(fig)
    review.write_csv('exploratory_window_checks.csv',checks)
    print('Wrote post-survey raw 5/2ms full sequences and three-seed raw detail card; seeds unchanged.')


def track(frames):
    # Preserve the B0 control flow, generalized only in the number of forward blocks.
    rows=[]
    for duration in [5,2]:
        active=frames[duration]
        images=[review.b0.image(f['events']) for f in active]
        for feature,seed in SEEDS.items():
            rr=forward_general(active,images,seed,duration,feature)
            for r in rr:
                if r['position_valid']:
                    r.update(review.b0.geometry(active[int(r['frame'])]['events'],[r['anchor_x'],r['anchor_y']]))
                    r.update(review.morphology(active[int(r['frame'])]['events'],[r['anchor_x'],r['anchor_y']]))
            rows.extend(rr)
            print(duration,feature,[(r['frame'],r['status']) for r in rr],flush=True)
    review.write_csv('exploratory_observations.csv',rows)
    metadata=dict(selection=CONFIG['selection'],normal_forward_evaluation_ms=50,seeds=SEEDS,durations_ms=[2,5],
                  seed_windows={'5ms':[11.750,11.755],'2ms':[11.750,11.752]},
                  gates='Exact unmodified B0 matching/rejection gates; geometry diagnostic only',
                  configuration_sha256=hashlib.sha256((OUT/'exploratory_protocol.json').read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  B0_run_sha256=hashlib.sha256((review.B0/'run.py').read_bytes()).hexdigest())
    (OUT/'exploratory_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    plot_tracks(frames,rows)


def forward_general(frames,images,seed,duration,feature):
    reference=len(frames)//2;center=np.array(seed,dtype=int);previous=reference;active=True;rows=[]
    for k in range(reference,len(frames)):
        row=dict(segment='post_survey_11_75',duration_ms=duration,feature=feature,frame=k,reference_frame=reference,
                 start_s=frames[k]['start_s'],end_s=frames[k]['end_s'],midpoint_s=(frames[k]['start_s']+frames[k]['end_s'])/2,
                 anchor_x=np.nan,anchor_y=np.nan,position_valid=False,max_frame_consumed=k)
        if k==reference:
            enough=len(review.b0.unique_xy(review.b0.local_events(frames[k]['events'],center)))>=12
            row.update(status='manual_seed' if enough else 'manual_seed_insufficient_pixels',anchor_x=int(center[0]),anchor_y=int(center[1]),position_valid=enough)
            active=enough;rows.append(row);continue
        if not active:row.update(status='not_attempted_after_loss');rows.append(row);continue
        d=review.b0.match(images[previous],images[k],center,center);reasons=[]
        if not d['match_available']:reasons=['image_boundary'];proposal=center
        else:
            proposal=np.array([d['proposed_x'],d['proposed_y']]);source_count=len(review.b0.unique_xy(review.b0.local_events(frames[previous]['events'],center)))
            target_count=len(review.b0.unique_xy(review.b0.local_events(frames[k]['events'],proposal)))
            rev=review.b0.match(images[k],images[previous],proposal,center)
            fb=float(np.hypot(rev['proposed_x']-center[0],rev['proposed_y']-center[1])) if rev['match_available'] else np.inf
            d.update(source_distinct_pixels=source_count,target_distinct_pixels=target_count,forward_backward_px=fb,
                     reverse_score=rev.get('score',np.nan),source_frame=previous,direction=1)
            if min(source_count,target_count)<12:reasons.append('insufficient_pixels')
            if d['score']<.35:reasons.append('score')
            if d['score_margin']<.025:reasons.append('ambiguous_peak')
            if d['search_boundary']:reasons.append('search_boundary')
            if fb>2:reasons.append('forward_backward')
        row.update(d)
        if reasons:active=False;row['status']='rejected_'+'+'.join(reasons)
        else:row.update(status='accepted',position_valid=True,anchor_x=int(proposal[0]),anchor_y=int(proposal[1]));center=proposal;previous=k
        rows.append(row)
    return rows


def plot_tracks(frames,rows):
    fig,axes=plt.subplots(2,5,figsize=(18,9));colors={'L':'#b03090','A1':'#008b45','A2':'#7755cc'}
    for i,duration in enumerate([5,2]):
        chosen=[10,12,14,16,19] if duration==5 else [25,30,35,40,49]
        for ax,k in zip(axes[i],chosen):
            f=frames[duration][k];review.b0.raw.draw_raw(ax,review.b0.raw.crop(f['events']),2)
            missing=[]
            for feature in SEEDS:
                r=next(r for r in rows if r['duration_ms']==duration and r['feature']==feature and r['frame']==k)
                if r['position_valid']:
                    p=np.array([r['anchor_x'],r['anchor_y']]);ax.plot(*p,'x',color=colors[feature]);ax.text(p[0]+1,p[1]-2,feature,color=colors[feature],fontsize=8)
                    ax.add_patch(plt.Rectangle(p-13,26,26,fill=False,color=colors[feature],lw=.7))
                else:missing.append(feature+' lost')
            ax.set_title(f'{duration}ms f{k} [{f["start_s"]:.3f},{f["end_s"]:.3f})s',fontsize=9)
            ax.set_xlabel('; '.join(missing),fontsize=8,color='#aa3333')
    fig.suptitle('Post-survey exploratory forward-only observations: no auxiliary fills or identity claims',fontsize=12)
    fig.subplots_adjust(left=.035,right=.995,bottom=.06,top=.9,wspace=.18,hspace=.45)
    fig.savefig(OUT/'exploratory_tracking_review.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--track',action='store_true');args=parser.parse_args()
    frames,checks=windows()
    if args.track:track(frames)
    else:preview(frames,checks)
