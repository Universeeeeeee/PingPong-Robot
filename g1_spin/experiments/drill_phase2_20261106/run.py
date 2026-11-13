"""Stage 2: independent output, frozen estimator, no physical truth claimed."""
import argparse
import csv
import hashlib
import importlib.util
import json
import platform
from pathlib import Path
from time import perf_counter

import cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
OLD=OUT.parent/'drill_phase1_20261106'
ARCHIVE=ROOT/'g1_spin/data/archive_20261105'
NAMES={1:'DSP03F100025_20261105_201058',2:'DSP03F100025_20261105_201123',
       3:'DSP03F100025_20261105_201143',4:'DSP03F100025_20261105_203829',5:'DSP03F100025_20261105_215030'}
BLOCKS={2:([.8,2.0],[3.4,5.2]),3:([1.4,2.6],[4.2,6.2]),4:([.8,1.6],[2.4,3.2])}
CONFIG={2:dict(radius=63.,seed=[625.,495.],sector=None),
        3:dict(radius=59.,seed=[555.,430.],sector=[5,65]),
        4:dict(radius=60.,seed=[598.,355.],sector=[60,120])}
ANCHORS={2:[.9,3.5],3:[1.5,4.3],4:[.9,2.5]}
SPEED={1:[1.425,3.665,4.365,1.495,3.595,5.695],5:[4.095,5.095]}
ALPHAS={'A-1':-1.,'A0':0.,'A0.5':.5,'A1':1.}
spec=importlib.util.spec_from_file_location('phase2_frozen_cmax',OLD/'frozen/spherical_cmax.py')
cmax=importlib.util.module_from_spec(spec);spec.loader.exec_module(cmax)

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def write_json(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False))

def table(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load(run):
    folder=ARCHIVE/'decoded'/NAMES[run]
    m=json.loads((folder/'manifest.json').read_text())
    e=np.concatenate([np.load(folder/p,allow_pickle=False) for p in m['parts']])
    return e,int(e['t'][0])

def select(e,t0,start_us,duration_us):
    a,b=np.searchsorted(e['t'],[t0+int(start_us),t0+int(start_us)+int(duration_us)])
    return e[a:b],int(a)

def picture(q,box):
    x0,y0,x1,y1=map(int,box)
    q=q[(q['x']>=x0)&(q['x']<x1)&(q['y']>=y0)&(q['y']<y1)]
    im=np.full((y1-y0,x1-x0,3),18,np.uint8)
    for p,color in [(0,[75,175,255]),(1,[255,175,65])]:
        a=q[q['p']==p];im[a['y']-y0,a['x']-x0]=color
    return im

def centers(ts,xy,query):
    """Validate only each requested interpolation bracket; exact knots are valid."""
    query=np.asarray(query,dtype=np.int64)
    if not len(ts) or np.any(query<ts[0]) or np.any(query>ts[-1]):
        raise ValueError('geometry_outside_coverage')
    right=np.searchsorted(ts,query,side='left')
    exact=ts[right]==query
    left=np.maximum(right-1,0)
    if np.any((~exact)&(ts[right]-ts[left]>30000)):
        raise ValueError('geometry_gap_over_30ms')
    return np.column_stack([np.interp(query-ts[0],ts-ts[0],xy[:,i]) for i in range(2)])

def observations(run):
    if run in [1,5]:
        rows=list(csv.DictReader((OLD/('rotation_geometry.csv' if run==1 else 'control_geometry.csv')).open()))
    else:rows=[r for r in json.loads((OUT/'geometry.json').read_text()) if r['run']==run]
    good=[r for r in rows if str(r['geometry_valid'])=='True']
    ts=np.array([round(float(r['time_s'])*1e6) for r in good],np.int64)
    xy=np.array([[float(r['cx']),float(r['cy'])] for r in good])
    order=np.argsort(ts)
    return ts[order],xy[order]

def fit(q,seed,radius,sector):
    x0,y0=np.floor(np.asarray(seed)-105).astype(int)
    q=q[(q['x']>=x0)&(q['x']<x0+210)&(q['y']>=y0)&(q['y']<y0+210)]
    im=np.zeros((210,210),np.uint8);im[q['y']-y0,q['x']-x0]=255
    circles=cv2.HoughCircles(cv2.GaussianBlur(im,(5,5),1),cv2.HOUGH_GRADIENT,1,50,
                            param1=80,param2=20,minRadius=round(radius-7),maxRadius=round(radius+7))
    if circles is None:return None
    xy=np.unique(np.column_stack((q['x'],q['y'])),axis=0).astype(float);candidates=[]
    for cc in circles[0]:
        c=cc[:2]+[x0,y0];delta=xy-c;dist=np.linalg.norm(delta,axis=1)
        theta=np.arctan2(delta[:,1],delta[:,0]);bins=np.floor((theta+np.pi)/(2*np.pi)*36).astype(int)%36
        near=abs(dist-radius)<4
        if sector:near&=~((theta>np.deg2rad(sector[0]))&(theta<np.deg2rad(sector[1])))
        ids=[j for j in range(36) if np.count_nonzero(near&(bins==j))>=2]
        pts=np.array([np.median(xy[near&(bins==j)],axis=0) for j in ids])
        if len(pts)<15:continue
        opt=least_squares(lambda center:np.linalg.norm(pts-center,axis=1)-radius,c,
                          bounds=(c-5,c+5),loss='soft_l1',f_scale=1.)
        res=np.linalg.norm(pts-opt.x,axis=1)-radius;support=int((abs(res)<=2).sum())
        rms=float(np.sqrt(np.mean(res**2)))
        if support<15 or rms>2.5 or np.any(opt.active_mask):continue
        candidates.append(dict(cx=float(opt.x[0]),cy=float(opt.x[1]),hough_radius_px=float(cc[2]),
            ring_rms_px=rms,angular_bins=support,fit_bins=len(pts),support_bins=[j for j,r in zip(ids,res) if abs(r)<=2]))
    return max(candidates,key=lambda c:(c['angular_bins'],-c['ring_rms_px'])) if candidates else None

def protected():
    paths=list(OLD.rglob('*'))+list((ROOT/'g1_spin/src').rglob('*.py'))
    paths += list((ARCHIVE/'raw').glob('*.raw'))
    for run in NAMES:
        folder=ARCHIVE/'decoded'/NAMES[run];m=json.loads((folder/'manifest.json').read_text())
        paths += [folder/'manifest.json']+[folder/p for p in m['parts']]
    paths += [ARCHIVE/'recording_labels.json',ARCHIVE/'verification.json']
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths)) if p.is_file() and '__pycache__' not in str(p)}

def geometry():
    if (OUT/'geometry.json').exists():raise FileExistsError('Preserve geometry')
    write_json(OUT/'source_hashes.json',protected())
    write_json(OUT/'config.json',{str(k):v for k,v in CONFIG.items()})
    allrows=[];anchorrows=[]
    fig,axes=plt.subplots(3,2,figsize=(12,15))
    for rowidx,run in enumerate([2,3,4]):
        e,t0=load(run);cfg=CONFIG[run];bytime={}
        for starts in BLOCKS[run]:
            for block in starts:
                seed=np.array(cfg['seed'])
                for t in range(round((block-.025)*1e6),round((block+.235)*1e6)+1,10000):
                    q,_=select(e,t0,t-2500,5000);result=fit(q,seed,cfg['radius'],cfg['sector'])
                    g=dict(run=run,time_s=t/1e6,start_s=(t-2500)/1e6,end_s=(t+2500)/1e6,
                        event_median_s=float(np.median(q['t']-t0)*1e-6) if len(q) else None,
                        geometry_valid=result is not None,cx=None,cy=None,radius_px=cfg['radius'],
                        hough_radius_px=None,ring_rms_px=None,angular_bins=0,fit_bins=0,support_bins=[])
                    if result:g.update(result);seed=np.array([result['cx'],result['cy']])
                    bytime[t]=g
        # Anchor observations are independently fitted at the exact specified times.
        for col,s in enumerate(ANCHORS[run]):
            t=round(s*1e6);q,_=select(e,t0,t-2500,5000)
            result=fit(q,cfg['seed'],cfg['radius'],cfg['sector']);center=np.array(cfg['seed']) if result is None else np.array([result['cx'],result['cy']])
            box=np.r_[np.floor(center-100).astype(int),np.floor(center+100).astype(int)]
            ax=axes[rowidx,col];ax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]])
            if result:ax.add_patch(plt.Circle(center,cfg['radius'],color='lime',fill=False,lw=.8))
            ax.set_title(f'run{run} {s:.3f}s | r={cfg["radius"]:g} | '+(f'{result["angular_bins"]} bins RMS {result["ring_rms_px"]:.2f}' if result else 'NO GEOMETRY'))
            anchorrows.append(dict(run=run,time_s=s,result=result,raw_events=len(q),radius_px=cfg['radius'],sector=cfg['sector']))
        allrows.extend(bytime[t] for t in sorted(bytime))
        print('geometry',run,len(bytime),sum(g['geometry_valid'] for g in bytime.values()),flush=True)
    write_json(OUT/'geometry.json',allrows);write_json(OUT/'geometry_anchors.json',anchorrows)
    table(OUT/'geometry.csv',[{k:v for k,v in g.items() if k!='support_bins'} for g in allrows])
    fig.suptitle('Pre-spin raw 5ms geometry anchors | candidates, not truth');fig.tight_layout();fig.savefig(OUT/'geometry_anchors.png',dpi=150);plt.close(fig)

def freeze():
    if (OUT/'freeze.json').exists():raise FileExistsError('Preserve freeze')
    assert (OUT/'geometry_review.json').exists()
    files=['protocol.md','run.py','config.json','geometry.json','geometry.csv','geometry_anchors.json','geometry_review.json','source_hashes.json']
    write_json(OUT/'freeze.json',dict(files={p:digest(OUT/p) for p in files},python=platform.python_version(),
        estimator_sha256=digest(OLD/'frozen/spherical_cmax.py'),max_rpm=3500,initial=None,
        planned_extension_windows=36,planned_ablation_windows=8,planned_recovered_windows=2,planned_result_rows=110))

def check_freeze():
    frozen=json.loads((OUT/'freeze.json').read_text())
    for p,h in frozen['files'].items():assert digest(OUT/p)==h,p

def prepare():
    check_freeze();(OUT/'windows').mkdir(exist_ok=True);audit=[]
    if (OUT/'window_audit.json').exists():raise FileExistsError('Preserve prepared windows')
    for run in NAMES:
        jobs=[]
        if run in BLOCKS:
            for part,starts in zip(['development','evaluation'],BLOCKS[run]):
                for b in starts:
                    for offset in [.025,.095,.165]:jobs.append(('extension',round((b+offset)*1e6),part,round((b+.1)*1e6)))
        if run==1:jobs += [('recovered',725000,'supplement',800000),('recovered',5765000,'supplement',5700000)]
        jobs += [('ablation',round(s*1e6),'diagnostic',None) for s in SPEED.get(run,[])]
        if not jobs:continue
        e,t0=load(run);gt,gxy=observations(run)
        for kind,begin,part,anchor in jobs:
            end=begin+10000;ref=(begin+end)//2;name=f'run{run}_{begin:08d}_{kind}'
            raw,a=select(e,t0,begin,10000);xy=np.column_stack([raw['x'],raw['y']]).astype(float);times=raw['t'].astype(np.int64)-t0
            radius=CONFIG[run]['radius'] if kind=='extension' else 60.
            variants=['B1','B2'] if kind=='extension' else (['B0','B1','B2'] if kind=='recovered' else list(ALPHAS))
            row=dict(window=name,run=run,kind=kind,partition=part,start_us=begin,end_us=end,reference_us=ref,
                anchor_us=anchor,radius_px=radius,raw_count=len(raw),geometry_status='candidate',failure=None,
                optimizer_count=0,common_count=0,supports={},common_fractions={},variants=variants,
                input_sha256=None,hypothesis_only=True,physics_validated=False)
            try:
                if kind=='ablation':
                    oldname=('rotation' if run==1 else 'translation')+f'_{begin:08d}'
                    d=np.load(OLD/'windows'/f'{oldname}.npz');ids=d['ids'];common_ids=d['common_ids']
                    assert np.array_equal(e[ids]['t'].astype(np.int64)-t0,d['time_us'])
                    cs={v:d['B1']+alpha*(d['B2']-d['B1']) for v,alpha in ALPHAS.items()}
                    xy=d['xy'];times=d['time_us'];polarity=d['polarity']
                    row['supports']={v:len(ids) for v in variants};row['common_fractions']={v:1. for v in variants}
                    # Preserve the original common set; no new mask or resampling for alpha.
                    assert all(np.all(np.linalg.norm(xy-c,axis=1)<radius) for c in cs.values())
                else:
                    # Full-window endpoints and every actual event query are checked.
                    centers(gt,gxy,[begin,end,ref]);c1=centers(gt,gxy,[ref])[0]
                    cs={'B1':np.tile(c1,(len(raw),1)),'B2':centers(gt,gxy,times)}
                    if kind=='recovered':cs['B0']=np.tile(centers(gt,gxy,[anchor])[0],(len(raw),1))
                    dist={v:np.linalg.norm(xy-c,axis=1) for v,c in cs.items()}
                    row['supports']={v:int(np.count_nonzero(d<=.8*radius)) for v,d in dist.items()}
                    limit=.8*radius if kind=='extension' else 44.4
                    common=np.logical_and.reduce([d<=limit for d in dist.values()]);sel=np.flatnonzero(common);common_ids=sel+a
                    row['common_fractions']={v:len(sel)/n if n else None for v,n in row['supports'].items()}
                    if len(sel)>5000:sel=sel[np.linspace(0,len(sel)-1,5000).astype(int)]
                    ids=sel+a;xy=xy[sel];times=times[sel];polarity=raw['p'][sel];cs={v:c[sel] for v,c in cs.items()}
                row.update(common_count=len(common_ids),optimizer_count=len(ids))
                path=OUT/'windows'/f'{name}.npz'
                np.savez_compressed(path,ids=ids,common_ids=common_ids,xy=xy,time_us=times,polarity=polarity,origin_us=t0,**cs)
                row['input_sha256']=digest(path)
            except ValueError as exc:row.update(geometry_status='unavailable',failure=str(exc))
            audit.append(row)
        print('prepared',run,len(jobs),flush=True)
    assert len(audit)==46
    write_json(OUT/'window_audit.json',audit);table(OUT/'window_audit.csv',audit)

def run_estimate(points,seconds):
    calls=[];original=cmax.minimize
    def trace(*args,**kwargs):
        result=original(*args,**kwargs)
        calls.append(dict(seed=np.asarray(args[1]).tolist(),omega=result.x.tolist(),fun=float(result.fun),success=bool(result.success),iterations=int(result.nit)))
        return result
    cmax.minimize=trace;start=perf_counter()
    try:result=cmax.estimate(points,seconds,max_rpm=3500,initial=None)
    finally:cmax.minimize=original
    elapsed=perf_counter()-start
    if len(seconds):
        dt=seconds-(seconds[0]+seconds[-1])/2
        for c in calls:c['contrast']=float(np.var(cmax.iwe(cmax.warp(points,dt,c['omega']))/len(points)))
    return result,calls,elapsed

def optimize():
    check_freeze();path=OUT/'results.jsonl'
    if path.exists():raise FileExistsError('Preserve results')
    rows=[]
    with path.open('w') as out:
        for w in json.loads((OUT/'window_audit.json').read_text()):
            d=None if w['geometry_status']=='unavailable' else np.load(OUT/'windows'/f'{w["window"]}.npz')
            for v in w['variants']:
                result=None;calls=[];elapsed=0.
                if d is not None:
                    rel=d['xy']-d[v];points=cmax.lift(rel[:,0],rel[:,1],0,0,w['radius_px'])
                    result,calls,elapsed=run_estimate(points,(d['time_us']-w['start_us'])*1e-6)
                row={k:w[k] for k in ['window','run','kind','partition','start_us','end_us','radius_px','geometry_status','optimizer_count']}
                omega=None if result is None else (result['omega_rad_s'] if result['valid'] else result['quality'].get('rejected_candidate_omega_rad_s'))
                row.update(variant=v,valid=bool(result and result['valid']),failure_reason='geometry_unavailable' if result is None else result['failure_reason'],
                    omega=omega,axis=None if omega is None or np.linalg.norm(omega)<1e-12 else (np.asarray(omega)/np.linalg.norm(omega)).tolist(),
                    candidate_rpm=None if omega is None else float(np.linalg.norm(omega)*60/(2*np.pi)),
                    rpm=None if result is None else result['spin_rpm'],gain=None if result is None else result['quality'].get('relative_improvement'),
                    runtime_s=elapsed,result=result,optimizer_candidates=calls,hypothesis_only=True,physics_validated=False)
                out.write(json.dumps(row,allow_nan=False)+'\n');out.flush();rows.append(row)
            print(w['window'],w['geometry_status'],w['optimizer_count'],flush=True)
    assert len(rows)==110
    table(OUT/'results.csv',[{k:v for k,v in r.items() if k not in ['result','optimizer_candidates']} for r in rows])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['geometry','freeze','prepare','optimize']);args=parser.parse_args()
    globals()[args.stage]()
