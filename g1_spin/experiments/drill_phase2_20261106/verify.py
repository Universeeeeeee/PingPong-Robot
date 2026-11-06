"""Independent complete-file event/geometry/xyz checks for all 46 planned windows."""
import hashlib
import json
import numpy as np
from run import OUT,OLD,ROOT,NAMES,ALPHAS,load,centers,cmax,digest,write_json,check_freeze

def interpolation(ts,xy,query):
    values=[]
    for t in query:
        right=int(np.searchsorted(ts,t,side='left'))
        if right<len(ts) and ts[right]==t:values.append(xy[right]);continue
        if right==0 or right==len(ts):raise ValueError('geometry_outside_coverage')
        left=right-1
        if ts[right]-ts[left]>30000:raise ValueError('geometry_gap_over_30ms')
        a=(int(t)-int(ts[left]))/(int(ts[right])-int(ts[left]))
        values.append(xy[left]*(1-a)+xy[right]*a)
    return np.asarray(values).reshape(-1,2)

def xyz(pixels,center,radius):
    xx=(pixels[:,0]-center[:,0])/radius;yy=(pixels[:,1]-center[:,1])/radius
    assert np.all(xx*xx+yy*yy<=1+1e-12)
    return np.column_stack([xx,yy,-np.sqrt(np.maximum(0,1-xx*xx-yy*yy))])

def sha_array(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def main():
    import csv
    check_freeze()
    ts=np.array([0,10000,100000,110000]);xy=np.column_stack([ts/1000,ts/2000])
    # The unused middle 90ms gap must not reject two locally supported queries.
    assert np.allclose(centers(ts,xy,[5000,105000]),[[5,2.5],[105,52.5]])
    assert np.array_equal(centers(ts,xy,[10000,100000]),xy[[1,2]])
    try:centers(ts,xy,[50000])
    except ValueError as e:assert str(e)=='geometry_gap_over_30ms'
    else:raise AssertionError('Actually crossed gap was accepted')
    audit=json.loads((OUT/'window_audit.json').read_text());details=[json.loads(l) for l in (OUT/'results.jsonl').read_text().splitlines()]
    assert len(audit)==46 and len(details)==110
    assert len({(r['window'],r['variant']) for r in details})==110
    checks=[];max_center=0.;max_xyz=0.;replay=None;contexts=0
    for run in NAMES:
        e,t0=load(run)
        if run in [1,5]:
            gs=list(csv.DictReader((OLD/('rotation_geometry.csv' if run==1 else 'control_geometry.csv')).open()))
        else:gs=[g for g in json.loads((OUT/'geometry.json').read_text()) if g['run']==run]
        good=[g for g in gs if str(g['geometry_valid'])=='True']
        ts=np.array([round(float(g['time_s'])*1e6) for g in good]);coords=np.array([[float(g['cx']),float(g['cy'])] for g in good]);ix=np.argsort(ts);ts=ts[ix];coords=coords[ix]
        for w in [w for w in audit if w['run']==run]:
            rawids=np.flatnonzero((e['t']>=t0+w['start_us'])&(e['t']<t0+w['end_us']))
            raw=e[rawids];times=raw['t'].astype(np.int64)-t0;pixels=np.column_stack([raw['x'],raw['y']]).astype(float)
            failure=None
            try:
                if w['kind']=='ablation':
                    oldname=('rotation' if run==1 else 'translation')+f'_{w["start_us"]:08d}'
                    old=np.load(OLD/'windows'/f'{oldname}.npz')
                    ids=old['ids'];common_ids=old['common_ids'];expected_pixels=np.column_stack([e[ids]['x'],e[ids]['y']]).astype(float)
                    tt=e[ids]['t'].astype(np.int64)-t0
                    c1=interpolation(ts,coords,[w['reference_us']])[0]
                    cb2=interpolation(ts,coords,tt)
                    cs={v:c1[None,:]+a*(cb2-c1[None,:]) for v,a in ALPHAS.items()}
                else:
                    interpolation(ts,coords,[w['start_us'],w['end_us'],w['reference_us']])
                    c1=interpolation(ts,coords,[w['reference_us']])[0]
                    cs={'B1':np.tile(c1,(len(raw),1)),'B2':interpolation(ts,coords,times)}
                    if w['kind']=='recovered':cs['B0']=np.tile(interpolation(ts,coords,[w['anchor_us']])[0],(len(raw),1))
                    limit=.8*w['radius_px'] if w['kind']=='extension' else 44.4
                    mask=np.ones(len(raw),bool)
                    for c in cs.values():mask&=np.sum((pixels-c)**2,axis=1)<=limit**2
                    common_ids=rawids[mask];sel=np.flatnonzero(mask)
                    if len(sel)>5000:sel=sel[np.linspace(0,len(sel)-1,5000).astype(int)]
                    ids=rawids[sel];expected_pixels=pixels[sel];tt=times[sel];cs={v:c[sel] for v,c in cs.items()}
            except ValueError as exc:failure=str(exc)
            if failure:
                assert w['geometry_status']=='unavailable' and w['failure']==failure
                assert all(r['failure_reason']=='geometry_unavailable' for r in details if r['window']==w['window'])
                checks.append(dict(window=w['window'],status='verified_failure',reason=failure));continue
            assert w['geometry_status']=='candidate'
            path=OUT/'windows'/f'{w["window"]}.npz';assert digest(path)==w['input_sha256'];d=np.load(path)
            assert np.array_equal(d['ids'],ids);assert np.array_equal(d['common_ids'],common_ids)
            assert np.array_equal(d['xy'],expected_pixels);assert np.array_equal(d['time_us'],tt)
            assert np.array_equal(d['polarity'],e[ids]['p']);assert len(ids)==w['optimizer_count']
            hashes={}
            for v,c in cs.items():
                err=float(np.max(abs(d[v]-c))) if len(ids) else 0.;assert err<1e-10;max_center=max(max_center,err)
                independent=xyz(expected_pixels,c,w['radius_px']);rel=d['xy']-d[v]
                actual=cmax.lift(rel[:,0],rel[:,1],0,0,w['radius_px'])
                err3=float(np.max(abs(actual-independent))) if len(ids) else 0.;assert err3<1e-10;max_xyz=max(max_xyz,err3)
                assert np.all(independent[:,2]<=0);assert np.allclose(np.linalg.norm(independent,axis=1),1,atol=1e-12)
                hashes[v]=sha_array(actual)
                if replay is None and len(ids)>=100 and v=='B1':
                    row=next(r for r in details if r['window']==w['window'] and r['variant']==v)
                    result=cmax.estimate(actual,(tt-w['start_us'])*1e-6,max_rpm=3500,initial=None)
                    assert result==row['result'];replay=dict(window=w['window'],variant=v,exact_match=True)
            checks.append(dict(window=w['window'],status='verified_input',events=len(ids),xyz_sha256=hashes))
        for context in json.loads((OUT/'raw_diagnostics.json').read_text()):
            if not context['window'].startswith(f'run{run}_'):continue
            ids=np.flatnonzero((e['t']>=t0+context['raw_context_start_us'])&(e['t']<t0+context['raw_context_end_us']))
            path=OUT/context['raw_context_file'];assert digest(path)==context['raw_context_sha256']
            saved=np.load(path);assert np.array_equal(saved['ids'],ids);assert np.array_equal(saved['events'],e[ids])
            assert int(saved['origin_us'])==t0;contexts+=1
        del e
    protected=json.loads((OUT/'source_hashes.json').read_text())
    for rel,h in protected.items():assert digest(ROOT/rel)==h,rel
    report=dict(status='passed',planned_windows=46,result_rows=110,verified_input_windows=sum(c['status']=='verified_input' for c in checks),
        verified_failure_windows=sum(c['status']=='verified_failure' for c in checks),max_center_difference_px=max_center,
        max_xyz_difference=max_xyz,negative_visible_z_and_unit_norm_checked=True,
        tests=['unused middle gap accepted','actual >30ms interpolation gap rejected','exact knots do not cross gap'],
        uninstrumented_replay=replay,raw_60ms_contexts_verified=contexts,protected_files_unchanged=len(protected),checks=checks,
        caveat='Independent computation checks; no independent center, axis or RPM truth.')
    write_json(OUT/'verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
