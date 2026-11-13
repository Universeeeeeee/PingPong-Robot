"""Frozen, paired B0/B1/B2 experiment. Geometry hypotheses are not ground truth."""
import argparse
import csv
import hashlib
import importlib.util
import json
import platform
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from explore import OUT, ROOT, DATA, RUNS, load

spec=importlib.util.spec_from_file_location('frozen_cmax',OUT/'frozen/spherical_cmax.py')
cmax=importlib.util.module_from_spec(spec);spec.loader.exec_module(cmax)
BASE=['B0','B1','B2']
SENS={f'{dim}{offset:+d}':(dim,offset) for dim in ['cx','cy','r'] for offset in [-2,-1,1,2]}
BLOCKS={'rotation':([.7,1.4,2.1],[3.5,4.2,5.6]),'static':([.5,1.5],[2.5,3.2]),
        'translation':([4.,4.4,5.],[5.6,6.,6.4])}
SENS_TIMES={'rotation':[1495000,3595000,5695000], 'translation':[4095000,5095000,6095000]}

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def write_json(p,value):p.write_text(json.dumps(value,indent=2,allow_nan=False))

def table(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def observations(key):
    rows=list(csv.DictReader((OUT/f'{key}_geometry.csv').open()))
    good=[r for r in rows if r['geometry_valid']=='True']
    t=np.array([round(float(r['time_s'])*1e6) for r in good],np.int64)
    xy=np.array([[float(r['cx']),float(r['cy'])] for r in good])
    return t,xy

def centers(t,xy,query):
    query=np.asarray(query,dtype=np.int64)
    if np.any(query<t[0]) or np.any(query>t[-1]):raise ValueError('geometry_outside_coverage')
    lo=max(0,int(np.searchsorted(t,query.min(),side='right'))-1)
    hi=min(len(t)-1,int(np.searchsorted(t,query.max(),side='left')))
    if np.any(np.diff(t[lo:hi+1])>30000):raise ValueError('geometry_gap_over_30ms')
    return np.column_stack([np.interp(query-t[0],t-t[0],xy[:,dim]) for dim in range(2)])

def protected_hashes():
    paths=[OUT/'protocol.md',OUT/'geometry.py',OUT/'run.py',OUT/'rotation_geometry.csv',OUT/'control_geometry.csv']
    paths+=sorted((OUT/'frozen').glob('*.py'))
    for key in RUNS:
        folder=DATA/RUNS[key];paths.append(folder/'manifest.json')
        m=json.loads((folder/'manifest.json').read_text());paths.extend(folder/p for p in m['parts'])
        paths.append(ROOT/'g1_spin/data/archive_20261105/raw'/f'{RUNS[key]}.raw')
    return {str(p.relative_to(ROOT)):digest(p) for p in paths}

def prepare():
    (OUT/'windows').mkdir(exist_ok=True)
    records=[]
    for key in RUNS:
        e,t0=load(key);gt,gxy=observations(key)
        labels=['rotation'] if key=='rotation' else ['static','translation']
        for state in labels:
            for partition,starts in zip(['development','evaluation'],BLOCKS[state]):
                for block in starts:
                    anchor=round((block+.1)*1e6)
                    for offset in [.025,.095,.165]:
                        begin=round((block+offset)*1e6);end=begin+10000;ref=(begin+end)//2
                        name=f'{state}_{begin:08d}'
                        a,b=np.searchsorted(e['t'],[t0+begin,t0+end]);raw=e[a:b]
                        xy=np.column_stack([raw['x'],raw['y']]).astype(float)
                        times=raw['t'].astype(np.int64)-t0
                        row=dict(window=name,state=state,partition=partition,run=key,start_us=begin,end_us=end,
                            geometry_reference_us=ref,block_anchor_us=anchor,raw_count=len(raw),
                            geometry_status='candidate',geometry_failure=None,common_count=0,optimizer_count=0,
                            B0_support=None,B1_support=None,B2_support=None,common_fraction_B2=None,
                            center_motion_px=None,sensitivity=begin in SENS_TIMES.get(state,[]),
                            hypothesis_only=True,physics_validated=False)
                        try:
                            if state=='static':
                                cs={v:np.tile([714.,355.],(len(raw),1)) for v in BASE}
                                row['geometry_status']='manual_stationary_position_prior_not_detected'
                            else:
                                # Require the full nominal window, not only the timestamps of selected events.
                                centers(gt,gxy,np.array([begin,end,anchor,ref]))
                                c0=centers(gt,gxy,[anchor])[0];c1=centers(gt,gxy,[ref])[0]
                                cs={'B0':np.tile(c0,(len(raw),1)),'B1':np.tile(c1,(len(raw),1)),
                                    'B2':centers(gt,gxy,times)}
                            dist={v:np.linalg.norm(xy-cs[v],axis=1) for v in BASE}
                            for v in BASE:row[v+'_support']=int(np.count_nonzero(dist[v]<=48))
                            common=np.logical_and.reduce([dist[v]<=44.4 for v in BASE])
                            selected=np.flatnonzero(common);common_ids=selected+a
                            row['common_count']=len(selected)
                            row['common_fraction_B2']=len(selected)/row['B2_support'] if row['B2_support'] else None
                            if len(selected)>5000:selected=selected[np.linspace(0,len(selected)-1,5000).astype(int)]
                            row['optimizer_count']=len(selected)
                            row['center_motion_px']=float(np.linalg.norm(cs['B2'][-1]-cs['B2'][0])) if len(raw) else 0.
                            np.savez_compressed(OUT/'windows'/f'{name}.npz',ids=selected+a,common_ids=common_ids,
                                xy=xy[selected],time_us=times[selected],polarity=raw['p'][selected],
                                **{v:cs[v][selected] for v in BASE},origin_us=t0)
                        except ValueError as exc:
                            row['geometry_status']='unavailable';row['geometry_failure']=str(exc)
                        records.append(row)
    table(OUT/'window_audit.csv',records);write_json(OUT/'window_audit.json',records)
    environment=dict(python=platform.python_version(),platform=platform.platform(),numpy=np.__version__,
                     scipy=scipy.__version__,max_rpm=3500,window_us=10000,protected_hashes=protected_hashes())
    write_json(OUT/'environment.json',environment)
    for state in BLOCKS:
        rr=[r for r in records if r['state']==state]
        print(state,'planned',len(rr),'geometry_available',sum(r['geometry_status']!='unavailable' for r in rr),
              '>=100 common events',sum(r['optimizer_count']>=100 for r in rr),flush=True)

def run_estimate(points,seconds):
    calls=[];original=cmax.minimize
    def trace(*args,**kwargs):
        result=original(*args,**kwargs)
        calls.append(dict(seed=np.asarray(args[1]).tolist(),omega=result.x.tolist(),fun=float(result.fun),
                          success=bool(result.success),iterations=int(result.nit)))
        return result
    cmax.minimize=trace
    start=perf_counter()
    try:result=cmax.estimate(points,seconds,max_rpm=3500,initial=None)
    finally:cmax.minimize=original
    elapsed=perf_counter()-start
    if len(seconds):
        dt=seconds-(seconds[0]+seconds[-1])/2
        zero=float(np.var(cmax.iwe(points)/len(points)))
        for call in calls:
            call['contrast']=float(np.var(cmax.iwe(cmax.warp(points,dt,call['omega']))/len(points)))
    else:zero=None
    return result,calls,elapsed,zero

def optimize():
    env=json.loads((OUT/'environment.json').read_text());assert env['protected_hashes']==protected_hashes()
    audit=json.loads((OUT/'window_audit.json').read_text());rows=[]
    detail_file=OUT/'results.jsonl'
    if detail_file.exists():raise FileExistsError('Preserve prior run; do not overwrite results')
    with detail_file.open('w') as out:
        for w in audit:
            variants=BASE+list(SENS) if w['sensitivity'] else BASE
            data=None if w['geometry_status']=='unavailable' else np.load(OUT/'windows'/f'{w["window"]}.npz')
            for variant in variants:
                row={k:w[k] for k in ['window','state','partition','start_us','end_us','geometry_status','optimizer_count']}
                row.update(variant=variant,valid=False,failure_reason=None,wx=None,wy=None,wz=None,rpm=None,
                           candidate_rpm=None,zero_contrast=None,contrast=None,gain=None,runtime_s=0.,
                           competitor_count=0,optimizer_reference_us=None)
                if data is None:
                    row['failure_reason']='geometry_unavailable';result=None;calls=[]
                else:
                    center=data[variant if variant in BASE else 'B2'].copy();radius=60.
                    if variant in SENS:
                        dim,offset=SENS[variant]
                        if dim=='r':radius+=offset
                        else:center[:,0 if dim=='cx' else 1]+=offset
                    relative=data['xy']-center
                    points=cmax.lift(relative[:,0],relative[:,1],0.,0.,radius)
                    seconds=(data['time_us']-w['start_us'])*1e-6
                    result,calls,elapsed,zero=run_estimate(points,seconds)
                    row.update(valid=result['valid'],failure_reason=result['failure_reason'],runtime_s=elapsed,
                               zero_contrast=zero,competitor_count=len(calls))
                    if len(seconds):row['optimizer_reference_us']=float((int(data['time_us'][0])+int(data['time_us'][-1]))/2)
                    candidate=result['omega_rad_s'] if result['valid'] else result['quality'].get('rejected_candidate_omega_rad_s')
                    if candidate is not None:row['candidate_rpm']=float(np.linalg.norm(candidate)*60/(2*np.pi))
                    if result['valid']:row.update(wx=candidate[0],wy=candidate[1],wz=candidate[2],rpm=result['spin_rpm'])
                    row['contrast']=result['quality'].get('contrast_after')
                    row['gain']=result['quality'].get('relative_improvement')
                detail=dict(**row,result=result,optimizer_candidates=calls,hypothesis_only=True,physics_validated=False)
                out.write(json.dumps(detail,allow_nan=False)+'\n');out.flush();rows.append(row)
            print(w['window'],w['geometry_status'],w['optimizer_count'],'events',flush=True)
    table(OUT/'results.csv',rows)
    assert env['protected_hashes']==protected_hashes()
    print('DONE',len(rows),'rows; protected hashes unchanged',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','optimize']);args=parser.parse_args()
    prepare() if args.stage=='prepare' else optimize()
