"""B0 local feature observations and short associations; fixed protocol.md, no pose model."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares
from scipy.signal import correlate2d

import continuous_preview as raw
OUT=Path(__file__).resolve().parent
SEEDS={'early_turn': {'F1':[565,208],'F2':[595,191],'F3':[525,210]},
       'middle': {'F1':[540,223],'F2':[581,222],'F3':[550,181]},
       'late': {'F1':[533,220],'F2':[574,193],'F3':[579,236]}}
HALF,SEARCH=13,12
X0,Y0,W,H=450,120,221,161


def write_csv(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def image(events):
    use=events[(events['x']>=X0)&(events['x']<X0+W)&(events['y']>=Y0)&(events['y']<Y0+H)]
    out=np.zeros((2,H,W),dtype=float)
    for k,on in enumerate([False,True]):
        e=use[(use['p']>0)==on]
        out[k,e['y'].astype(int)-Y0,e['x'].astype(int)-X0]=1.
        out[k]=gaussian_filter(out[k],.8)
    return out


def patch(im,center,pad=HALF):
    x,y=np.rint(center).astype(int)-[X0,Y0]
    if x-pad<0 or x+pad>=W or y-pad<0 or y+pad>=H:return None
    return im[:,y-pad:y+pad+1,x-pad:x+pad+1]


def local_events(events,center):
    x,y=center
    return events[(abs(events['x'].astype(float)-x)<=HALF)&(abs(events['y'].astype(float)-y)<=HALF)]


def unique_xy(events):
    return np.unique(np.column_stack((events['x'],events['y'])).astype(float),axis=0)


def cosine(a,b):
    if a is None or b is None:return np.nan
    den=np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.sum(a*b)/den) if den>0 else np.nan


def match(source,target,source_center,expected_center):
    template=patch(source,source_center)
    search=patch(target,expected_center,HALF+SEARCH)
    if template is None or search is None:return dict(match_available=False)
    numerator=sum(correlate2d(search[k],template[k],mode='valid') for k in range(2))
    energy=sum(correlate2d(search[k]**2,np.ones((27,27)),mode='valid') for k in range(2))
    den=np.linalg.norm(template)*np.sqrt(np.maximum(0,energy))
    scores=np.divide(numerator,den,out=np.zeros_like(numerator),where=den>1e-12)
    iy,ix=np.unravel_index(np.argmax(scores),scores.shape)
    other=scores.copy();other[max(0,iy-3):iy+4,max(0,ix-3):ix+4]=-np.inf
    sy,sx=np.unravel_index(np.argmax(other),other.shape)
    best=np.rint(expected_center).astype(int)+[ix-SEARCH,iy-SEARCH]
    second=np.rint(expected_center).astype(int)+[sx-SEARCH,sy-SEARCH]
    return dict(match_available=True,proposed_x=int(best[0]),proposed_y=int(best[1]),
                score=float(scores[iy,ix]),second_score=float(other[sy,sx]),
                score_margin=float(scores[iy,ix]-other[sy,sx]),second_x=int(second[0]),second_y=int(second[1]),
                search_boundary=bool(ix in [0,2*SEARCH] or iy in [0,2*SEARCH]))


def model(points,anchor,ellipse=False):
    prefix='ellipse' if ellipse else 'circle'
    output={prefix+'_valid':False,prefix+'_reason':'insufficient_pixels'}
    if len(points)<15:return output
    anchor=np.asarray(anchor,dtype=float)
    if ellipse:
        cov=np.cov(points.T);vals,vec=np.linalg.eigh(cov);theta=float(np.arctan2(vec[1,-1],vec[0,-1]))
        theta=(theta+np.pi/2)%np.pi-np.pi/2
        axes=np.clip(np.sqrt(np.maximum(vals[::-1],1)*2),3.1,17.9)
        initial=np.r_[anchor,axes,theta];lower=np.r_[anchor-8,[3,3,-np.pi]];upper=np.r_[anchor+8,[18,18,np.pi]]
        def residual(p):
            d=points-p[:2];c,s=np.cos(p[4]),np.sin(p[4]);u=(c*d[:,0]+s*d[:,1])/p[2];v=(-s*d[:,0]+c*d[:,1])/p[3]
            return (np.hypot(u,v)-1)*(p[2]+p[3])/2
    else:
        initial=np.r_[anchor,np.clip(np.median(np.linalg.norm(points-anchor,axis=1)),3.1,17.9)]
        lower=np.r_[anchor-8,3];upper=np.r_[anchor+8,18]
        def residual(p):return np.linalg.norm(points-p[:2],axis=1)-p[2]
    result=least_squares(residual,initial,bounds=(lower,upper),loss='soft_l1',f_scale=1.,max_nfev=200)
    p=result.x;err=residual(p)
    d=points-p[:2]
    if ellipse:
        c,s=np.cos(p[4]),np.sin(p[4]);angles=np.arctan2((-s*d[:,0]+c*d[:,1])/p[3],(c*d[:,0]+s*d[:,1])/p[2])
        axis_ratio=max(p[2:4])/min(p[2:4])
    else:angles=np.arctan2(d[:,1],d[:,0]);axis_ratio=1.
    bins=np.floor((angles+np.pi)/(2*np.pi)*12).astype(int)%12
    occupied=np.bincount(bins,minlength=12)>0
    empty=''.join('0' if b else '1' for b in np.r_[occupied,occupied])
    gap=min(12,max(map(len,empty.split('0'))))*30
    reasons=[]
    if not result.success:reasons.append('optimizer')
    if np.any(result.active_mask):reasons.append('bounds')
    if occupied.sum()<8:reasons.append('angular_bins')
    if gap>90:reasons.append('angular_gap')
    if np.median(abs(err))>1.5:reasons.append('median_residual')
    if np.sqrt(np.mean(err**2))>2.5:reasons.append('rms')
    if axis_ratio>2:reasons.append('axis_ratio')
    output.update({prefix+'_valid':not reasons,prefix+'_reason':'valid' if not reasons else '+'.join(reasons),
                   prefix+'_raw_cx':p[0],prefix+'_raw_cy':p[1],prefix+'_cx':p[0] if not reasons else np.nan,
                   prefix+'_cy':p[1] if not reasons else np.nan,prefix+'_radius_or_a':p[2],
                   prefix+'_b':p[3] if ellipse else p[2],prefix+'_theta':p[4] if ellipse else 0.,
                   prefix+'_occupied_30deg_bins':int(occupied.sum()),prefix+'_max_empty_angle_deg':gap,
                   prefix+'_median_abs_residual':float(np.median(abs(err))),prefix+'_rms':float(np.sqrt(np.mean(err**2))),
                   prefix+'_axis_ratio':axis_ratio})
    return output


def geometry(events,center):
    use=local_events(events,center);points=unique_xy(use)
    out=dict(local_events=len(use),local_distinct_pixels=len(points),
             local_on_events=int(np.sum(use['p']>0)),local_off_events=int(np.sum(use['p']<=0)))
    for name,e in [('on',use[use['p']>0]),('off',use[use['p']<=0]),('combined',use)]:
        pts=unique_xy(e)
        out[name+'_distinct_pixels']=len(pts)
        out[name+'_centroid_x']=float(np.mean(pts[:,0])) if len(pts) else np.nan
        out[name+'_centroid_y']=float(np.mean(pts[:,1])) if len(pts) else np.nan
    return dict(**out,**model(points,center),**model(points,center,True))


def track(frames,images,seed,identity):
    reference=len(frames)//2
    rowmap={}
    def base(k):
        return dict(**identity,frame=k,reference_frame=reference,start_s=frames[k]['start_s'],end_s=frames[k]['end_s'],
                    midpoint_s=(frames[k]['start_s']+frames[k]['end_s'])/2,
                    anchor_x=np.nan,anchor_y=np.nan,position_valid=False)
    center=np.array(seed,dtype=int)
    enough=len(unique_xy(local_events(frames[reference]['events'],center)))>=12
    rowmap[reference]=dict(base(reference),status='manual_seed' if enough else 'manual_seed_insufficient_pixels',
                           anchor_x=center[0],anchor_y=center[1],position_valid=enough,**geometry(frames[reference]['events'],center))
    for direction in [-1,1]:
        previous=reference;center=np.array(seed,dtype=int);active=enough
        for k in range(reference+direction,len(frames) if direction>0 else -1,direction):
            row=base(k)
            if not active:
                rowmap[k]=dict(row,status='not_attempted_after_loss');continue
            matchdiag=match(images[previous],images[k],center,center)
            reasons=[]
            if not matchdiag['match_available']:
                reasons=['image_boundary'];proposal=center
            else:
                proposal=np.array([matchdiag['proposed_x'],matchdiag['proposed_y']])
                source_count=len(unique_xy(local_events(frames[previous]['events'],center)))
                target_count=len(unique_xy(local_events(frames[k]['events'],proposal)))
                reverse=match(images[k],images[previous],proposal,center)
                fb=float(np.hypot(reverse['proposed_x']-center[0],reverse['proposed_y']-center[1])) if reverse['match_available'] else np.inf
                matchdiag.update(source_distinct_pixels=source_count,target_distinct_pixels=target_count,forward_backward_px=fb,
                                 reverse_score=reverse.get('score',np.nan),source_frame=previous,direction=direction)
                if min(source_count,target_count)<12:reasons.append('insufficient_pixels')
                if matchdiag['score']<.35:reasons.append('score')
                if matchdiag['score_margin']<.025:reasons.append('ambiguous_peak')
                if matchdiag['search_boundary']:reasons.append('search_boundary')
                if fb>2:reasons.append('forward_backward')
            row.update(matchdiag)
            if reasons:
                active=False;row['status']='rejected_'+'+'.join(reasons)
            else:
                row.update(status='accepted',position_valid=True,anchor_x=int(proposal[0]),anchor_y=int(proposal[1]))
                row.update(geometry(frames[k]['events'],proposal))
                next_k=k+direction
                if 0<=next_k<len(frames):
                    predicted=proposal+(proposal-center)
                    future=patch(images[next_k],predicted)
                    row.update(future_frame=next_k,future_predicted_x=int(predicted[0]),future_predicted_y=int(predicted[1]),
                               future_predicted_patch_score=cosine(patch(images[k],proposal),future),
                               future_predicted_distinct_pixels=len(unique_xy(local_events(frames[next_k]['events'],predicted))))
                center=proposal;previous=k
            rowmap[k]=row
    return [rowmap[k] for k in range(len(frames))]


def semantic_checks():
    a=np.zeros((2,H,W));b=np.zeros_like(a)
    yy,xx=np.mgrid[:9,:9];pattern=((xx*7+yy*11)%5<2).astype(float)
    a[0,70:79,90:99]=pattern;a[1,73:82,93:102]=pattern[::-1]
    b[:,3:,4:]=a[:,:-3,:-4]
    result=match(a,b,[X0+96,Y0+76],[X0+96,Y0+76])
    assert result['proposed_x']==X0+100 and result['proposed_y']==Y0+79
    # A compact footprint keeps both copies entirely isolated inside each 27px test patch.
    compact=np.zeros_like(a);compact[0,75:78,95:98]=[[1,0,1],[0,1,0],[1,1,0]]
    compact[1,75:78,95:98]=[[0,1,0],[1,0,1],[0,0,1]]
    template=patch(compact,[X0+96,Y0+76]);repeated=np.zeros_like(a)
    for dx in [-8,8]:
        cx,cy=96+dx,76;repeated[:,cy-13:cy+14,cx-13:cx+14]+=template
    ambiguous=match(compact,repeated,[X0+96,Y0+76],[X0+96,Y0+76])
    assert ambiguous['score_margin']<.025
    empty=match(a,np.zeros_like(a),[X0+96,Y0+76],[X0+96,Y0+76])
    assert empty['score']<.35 and empty['score_margin']<.025
    dtype=np.dtype([('x','u2'),('y','u2'),('p','i2'),('t','i8')])
    sample=np.array([(550+dx,200+dy,1,0) for dx in range(-2,2) for dy in range(-2,2)],dtype=dtype)
    frames=[dict(events=sample if k!=3 else sample[:0],start_s=k*.002,end_s=(k+1)*.002) for k in range(5)]
    tracked=track(frames,[image(f['events']) for f in frames],[550,200],dict(segment='synthetic',duration_ms=2,feature='test'))
    assert tracked[3]['status'].startswith('rejected_') and tracked[4]['status']=='not_attempted_after_loss'
    return dict(integer_translation_recovered=[4,3],duplicated_patch_rejected_by_margin=True,
                empty_patch_rejected=True,reappearance_after_loss_not_reconnected=True)


def main():
    checks=semantic_checks();events,manifest=raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);rows=[];cache={}
    seeds=dict(candidates=SEEDS,meaning='Human-reviewed local patch anchors, not geometric centers or physical Marker identities',
               basis={'early_turn':{'F1':'central split-polarity structure','F2':'upper-right mixed structure, overlap stress','F3':'left separate small structure'},
                      'middle':{'F1':'lower-left structure','F2':'lower-right structure','F3':'upper thin arc stress'},
                      'late':{'F1':'left lower structure','F2':'upper-right structure','F3':'lower near-occlusion stress'}},
               patch_size=27,review='Parent reviewed all three complete 5ms contact sheets and froze all nine proposals before tracking')
    (OUT/'seeds.json').write_text(json.dumps(seeds,indent=2)+'\n')
    for segment,(start,end) in raw.SEGMENTS.items():
        for duration in [2,5]:
            frames=[]
            for k in range(100//duration):
                lo=t0+round(start*1e6)+k*duration*1000;hi=lo+duration*1000;a,b=np.searchsorted(times,[lo,hi])
                frames.append(dict(events=events[a:b],start_s=(lo-t0)/1e6,end_s=(hi-t0)/1e6))
            images=[image(f['events']) for f in frames];cache[segment,duration]=frames
            for feature,seed in SEEDS[segment].items():
                rr=track(frames,images,seed,dict(segment=segment,duration_ms=duration,feature=feature));rows.extend(rr)
                print(segment,duration,feature,{s:sum(r['status']==s for r in rr) for s in ['accepted','manual_seed','not_attempted_after_loss']},flush=True)
    write_csv('observations.csv',rows)
    metadata=dict(source=str(raw.DATA),source_t0_us=t0,manifest_sha256=hashlib.sha256((raw.DATA/'manifest.json').read_bytes()).hexdigest(),
                  input_parts_sha256={f:hashlib.sha256((raw.DATA/f).read_bytes()).hexdigest() for f in manifest['parts']},
                  segments=raw.SEGMENTS,durations_ms=[2,5],expected_observation_rows=630,actual_observation_rows=len(rows),
                  patch_size=27,search_halfwidth_px=12,occupancy_blur_sigma=.8,min_distinct_pixels=12,
                  min_cosine_score=.35,min_second_peak_margin=.025,second_peak_exclusion_chebyshev_px=3,
                  max_forward_backward_px=2,stop_after_rejection=True,geometry_is_diagnostic_only=True,
                  layout_note='10-point Motive template found separately, not linked to current approximately 14-marker ball',
                  semantic_checks=checks,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  protocol_sha256=hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest(),
                  seeds_sha256=hashlib.sha256((OUT/'seeds.json').read_bytes()).hexdigest())
    assert len(rows)==630
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    plot(rows,cache)


def plot(rows,cache):
    colors={'F1':'#b03090','F2':'#008b45','F3':'#7b55cc'}
    fig,axes=plt.subplots(3,2,figsize=(12,11))
    for i,segment in enumerate(raw.SEGMENTS):
        for j,duration in enumerate([2,5]):
            ax=axes[i,j];frames=cache[segment,duration]
            selected=frames[len(frames)//2]['events'];raw.draw_raw(ax,raw.crop(selected),2)
            for feature,seed in SEEDS[segment].items():
                q=[r for r in rows if r['segment']==segment and r['duration_ms']==duration and r['feature']==feature]
                x=[r['anchor_x'] for r in q];y=[r['anchor_y'] for r in q]
                ax.plot(x,y,'o-',color=colors[feature],ms=2,lw=.7,label=feature)
                ax.plot(*seed,'x',color=colors[feature],ms=6);ax.text(seed[0]+2,seed[1]-2,feature,color=colors[feature],fontsize=9)
                ax.add_patch(plt.Rectangle(np.asarray(seed)-13,26,26,fill=False,edgecolor=colors[feature],linestyle=':',lw=.7))
            ax.set_title(f'{segment} {duration}ms: accepted patch-anchor trajectories\nbackground = reference raw block; boxes = fixed seed patches',fontsize=10)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('F1/F2/F3 are local feature names, not physical Marker IDs or geometric centers',fontsize=12)
    fig.subplots_adjust(left=.05,right=.98,bottom=.05,top=.92,wspace=.18,hspace=.4)
    fig.savefig(OUT/'trajectories.png',dpi=180);plt.close(fig)


if __name__=='__main__':main()
