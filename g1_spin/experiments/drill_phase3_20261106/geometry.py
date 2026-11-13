"""New block geometry, all old inputs copied read-only, and pre-spin freeze."""
import argparse
import csv
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import OUT,P1,P2,base,load,select,picture,fit,centers,digest,protected,write_json,table,array_hash

def geometry():
    if (OUT/'geometry.json').exists():raise FileExistsError('Preserve new geometry')
    write_json(OUT/'source_hashes.json',protected());cfg=json.loads((OUT/'config.json').read_text());allrows=[]
    fig,axes=plt.subplots(2,4,figsize=(17,10))
    for ax,c in zip(axes.flat,cfg):
        e,t0=load(c['run']);b=c['block_s'];seed=np.array(c['seed'],float);rows=[]
        for t in range(round((b-.025)*1e6),round((b+.235)*1e6)+1,10000):
            q,_=select(e,t0,t-2500,5000);result=fit(q,seed,c['radius_px'],c['sector_deg'])
            row=dict(run=c['run'],block_s=b,time_s=t/1e6,start_us=t-2500,end_us=t+2500,
                geometry_valid=result is not None,cx=None,cy=None,radius_px=c['radius_px'],hough_radius_px=None,
                ring_rms_px=None,angular_bins=0,fit_bins=0,support_bins=[])
            if result:row.update(result);seed=np.array([result['cx'],result['cy']])
            rows.append(row)
        allrows+=rows;s=round((b+.1)*1e6);q,_=select(e,t0,s-2500,5000)
        good=[r for r in rows if r['geometry_valid']]
        center=np.array(c['seed'],float) if not good else np.array([[r['cx'],r['cy']] for r in good])[np.argmin(abs(np.array([r['time_s'] for r in good])-(b+.1)))]
        box=np.r_[np.floor(center-105).astype(int),np.ceil(center+105).astype(int)]
        ax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]])
        if good:ax.add_patch(plt.Circle(center,c['radius_px'],fill=False,color='lime',lw=.8))
        ax.set_title(f'run{c["run"]} block{b:g} | r{c["radius_px"]} | {len(good)}/{len(rows)} candidates',fontsize=10)
        print('geometry',c['run'],b,len(good),'/',len(rows),flush=True)
    write_json(OUT/'geometry.json',allrows);table(OUT/'geometry.csv',[{k:v for k,v in r.items() if k!='support_bins'} for r in allrows])
    fig.suptitle('Before spin: raw block anchors with nearest accepted geometry candidate; missing remains missing')
    fig.tight_layout();fig.savefig(OUT/'new_geometry_review.png',dpi=150);plt.close(fig)

def prepare():
    if (OUT/'window_audit.json').exists():raise FileExistsError('Preserve prepared inputs')
    (OUT/'windows').mkdir(exist_ok=True);audit=[]
    old1=json.loads((P1/'window_audit.json').read_text());old2=json.loads((P2/'window_audit.json').read_text())
    jobs=[]
    for w in old1:
        run=1 if w['run']=='rotation' else 5;begin=w['start_us'];source=P1/'windows'/f'{w["window"]}.npz';geo='phase1_candidate' if w['state']!='static' else 'phase1_manual_static_prior'
        status=w['geometry_status'];failure=w['geometry_failure'];replacement=False
        if run==1 and begin in [725000,5765000]:
            v=next(v for v in old2 if v['kind']=='recovered' and v['start_us']==begin)
            source=P2/'windows'/f'{v["window"]}.npz';status=v['geometry_status'];failure=v['failure'];geo='phase2_recovered_phase1_candidate';replacement=True
        jobs.append(dict(run=run,start_us=begin,end_us=w['end_us'],condition={'rotation':'fixed_vertical','static':'static_prior','translation':'translation_early'}[w['state']],
            partition=w['partition'],geometry_source=geo,radius_px=60.,geometry_status='unavailable' if status=='unavailable' else 'available',
            geometry_failure=failure,source_file=str(source),recovered_replacement=replacement,source_mask='phase1_B0B1B2_intersection44.4px'))
    for w in old2:
        if w['kind']!='extension':continue
        jobs.append(dict(run=w['run'],start_us=w['start_us'],end_us=w['end_us'],condition={2:'camera_facing',3:'oblique_fixed',4:'reverse_fixed'}[w['run']],
            partition=w['partition'],geometry_source='phase2_extension_candidate',radius_px=w['radius_px'],geometry_status='unavailable' if w['geometry_status']=='unavailable' else 'available',
            geometry_failure=w['failure'],source_file=str(P2/'windows'/f'{w["window"]}.npz'),recovered_replacement=False,source_mask='phase2_B1B2_intersection0.8r'))
    for w in jobs:
        name=f'run{w["run"]}_{w["start_us"]:08d}';row=dict(window=name,**w,event_count=0,common_count=0,input_hash=None,input_file=None,input_file_sha256=None,center=None)
        if w['geometry_status']=='available':
            d=np.load(w['source_file']);path=OUT/'windows'/f'{name}.npz'
            np.savez_compressed(path,ids=d['ids'],common_ids=d['common_ids'],xy=d['xy'],time_us=d['time_us'],polarity=d['polarity'],center=d['B1'][0] if len(d['ids']) else np.zeros(2),origin_us=d['origin_us'],radius_px=w['radius_px'])
            row.update(event_count=len(d['ids']),common_count=len(d['common_ids']),input_hash=array_hash(d['ids'],d['xy'],d['time_us'],d['polarity'],d['B1'],np.array(w['radius_px'])),
                input_file=str(path.relative_to(OUT)),input_file_sha256=digest(path),center=d['B1'][0].tolist() if len(d['ids']) else None)
        audit.append(row)
    cfg=json.loads((OUT/'config.json').read_text());geom=json.loads((OUT/'geometry.json').read_text())
    for run in [4,5]:
        e,t0=load(run)
        for c in [c for c in cfg if c['run']==run]:
            g=[g for g in geom if g['run']==run and g['block_s']==c['block_s'] and g['geometry_valid']]
            ts=np.array([round(r['time_s']*1e6) for r in g]);coords=np.array([[r['cx'],r['cy']] for r in g])
            for offset in [.025,.095,.165]:
                begin=round((c['block_s']+offset)*1e6);end=begin+10000;ref=(begin+end)//2;name=f'run{run}_{begin:08d}'
                raw,a=select(e,t0,begin,10000);pixels=np.column_stack([raw['x'],raw['y']]).astype(float);times=raw['t'].astype(np.int64)-t0
                row=dict(window=name,run=run,start_us=begin,end_us=end,condition='dynamic_total_motion' if run==4 else 'translation_late',partition='fixed_extension',geometry_source='phase3_new_block_candidate',
                    radius_px=float(c['radius_px']),geometry_status='available',geometry_failure=None,source_file=None,recovered_replacement=False,source_mask='phase3_B1_0.8r',
                    block_s=c['block_s'],event_count=0,common_count=0,input_hash=None,input_file=None,input_file_sha256=None,center=None,raw_count=len(raw))
                try:
                    centers(ts,coords,[begin,end,ref]);center=centers(ts,coords,[ref])[0]
                    sel=np.flatnonzero(np.linalg.norm(pixels-center,axis=1)<=.8*c['radius_px']);common_ids=sel+a
                    if len(sel)>5000:sel=sel[np.linspace(0,len(sel)-1,5000).astype(int)]
                    path=OUT/'windows'/f'{name}.npz';ids=sel+a;xy=pixels[sel];tt=times[sel];pol=raw['p'][sel]
                    np.savez_compressed(path,ids=ids,common_ids=common_ids,xy=xy,time_us=tt,polarity=pol,center=center,origin_us=t0,radius_px=float(c['radius_px']))
                    row.update(event_count=len(ids),common_count=len(common_ids),input_hash=array_hash(ids,xy,tt,pol,np.tile(center,(len(ids),1)),np.array(float(c['radius_px']))),
                        input_file=str(path.relative_to(OUT)),input_file_sha256=digest(path),center=center.tolist())
                except ValueError as exc:row.update(geometry_status='unavailable',geometry_failure=str(exc))
                audit.append(row)
    assert len(audit)==108 and len({w['window'] for w in audit})==108
    write_json(OUT/'window_audit.json',audit)
    keys=list(dict.fromkeys(k for w in audit for k in w));table(OUT/'window_audit.csv',[{k:w.get(k) for k in keys} for w in audit])
    print('prepared',len(audit),'available geometry',sum(w['geometry_status']=='available' for w in audit),flush=True)

def freeze():
    if (OUT/'freeze.json').exists():raise FileExistsError('Preserve freeze')
    (OUT/'frozen').mkdir(exist_ok=True)
    for name in ['spherical_cmax.py','optical_flow.py']:shutil.copy2(P1/'frozen'/name,OUT/'frozen'/name)
    names=['protocol.md','config.json','geometry_prior_review.json','geometry_final_review.json','geometry.json','window_audit.json','source_hashes.json',
        'geometry.py','common.py','flow.py','methods.py','run.py','test_math.py','periodicity.py','math_verification.json','periodicity.json']
    names += ['frozen/spherical_cmax.py','frozen/optical_flow.py']
    write_json(OUT/'freeze.json',dict(files={p:digest(OUT/p) for p in names},windows=108,routes=['G','A','S'],result_rows=324,period_attempts=4,
        caveat='No independent physical truth; candidate/default distinctions retained.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['geometry','prepare','freeze']);args=parser.parse_args();globals()[args.stage]()
