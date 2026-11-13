"""Independent input/geometry reconstruction plus an uninstrumented estimator replay."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
DATA=ROOT/'g1_spin/data/archive_20261105/decoded'
NAMES={'rotation':'DSP03F100025_20261105_201058','control':'DSP03F100025_20261105_215030'}

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def interpolation(ts,coords,query):
    result=[]
    for t in query:
        right=int(np.searchsorted(ts,t,side='left'))
        if right<len(ts) and ts[right]==t:result.append(coords[right]);continue
        assert right>0 and right<len(ts)
        left=right-1;assert ts[right]-ts[left]<=30000
        alpha=(int(t)-int(ts[left]))/(int(ts[right])-int(ts[left]))
        result.append(coords[left]*(1-alpha)+coords[right]*alpha)
    return np.asarray(result)

def main():
    env=json.loads((OUT/'environment.json').read_text())
    for rel,h in env['protected_hashes'].items():assert digest(ROOT/rel)==h,rel
    for rel,h in json.loads((OUT/'baseline_hashes.json').read_text()).items():assert digest(ROOT/rel)==h,rel
    audit=json.loads((OUT/'window_audit.json').read_text());checks=[];max_center_error=0.
    for key in NAMES:
        folder=DATA/NAMES[key];m=json.loads((folder/'manifest.json').read_text())
        e=np.concatenate([np.load(folder/p,allow_pickle=False) for p in m['parts']]);origin=int(e['t'][0])
        geom=list(csv.DictReader((OUT/f'{key}_geometry.csv').open()));geom=[r for r in geom if r['geometry_valid']=='True']
        ts=np.array([round(float(r['time_s'])*1e6) for r in geom],np.int64)
        coords=np.array([[float(r['cx']),float(r['cy'])] for r in geom])
        for w in [w for w in audit if w['run']==key and w['geometry_status']!='unavailable']:
            # Whole-file boolean time selection independently checks searchsorted-based preparation.
            idx=np.flatnonzero((e['t']>=origin+w['start_us'])&(e['t']<origin+w['end_us']))
            raw=e[idx];t=raw['t'].astype(np.int64)-origin
            xy=np.column_stack([raw['x'],raw['y']]).astype(float)
            if w['state']=='static':centers={v:np.tile([714.,355.],(len(raw),1)) for v in ['B0','B1','B2']}
            else:
                c0=interpolation(ts,coords,[w['block_anchor_us']])[0]
                c1=interpolation(ts,coords,[w['geometry_reference_us']])[0]
                centers={'B0':np.tile(c0,(len(raw),1)),'B1':np.tile(c1,(len(raw),1)),
                         'B2':interpolation(ts,coords,t)}
            common=np.ones(len(raw),bool)
            for v in centers:common&=np.sum((xy-centers[v])**2,axis=1)<=44.4**2
            common_ids=idx[common];sel=np.flatnonzero(common)
            if len(sel)>5000:sel=sel[np.linspace(0,len(sel)-1,5000).astype(int)]
            d=np.load(OUT/'windows'/f'{w["window"]}.npz',allow_pickle=False)
            assert np.array_equal(d['common_ids'],common_ids)
            assert np.array_equal(d['ids'],idx[sel]);assert np.array_equal(d['xy'],xy[sel])
            assert np.array_equal(d['time_us'],t[sel]);assert np.array_equal(d['polarity'],raw['p'][sel])
            for v in centers:
                err=float(np.max(np.abs(d[v]-centers[v][sel]))) if len(sel) else 0.
                assert err<1e-10;max_center_error=max(max_center_error,err)
                ref=centers[v][sel][len(sel)//2] if len(sel) else np.zeros(2)
                aligned=xy[sel]-centers[v][sel]+ref
                assert np.allclose((aligned-ref)/60,(xy[sel]-centers[v][sel])/60,atol=1e-12)
            checks.append(dict(window=w['window'],matched_events=len(sel)))
    details=[json.loads(l) for l in (OUT/'results.jsonl').read_text().splitlines()]
    assert len(details)==216 and len({(r['window'],r['variant']) for r in details})==216
    # Replay exactly one completed main call without the minimize recorder.
    chosen=next(r for r in details if r['state']=='rotation' and r['variant']=='B2' and r['optimizer_count']>=100)
    d=np.load(OUT/'windows'/f'{chosen["window"]}.npz');relative=d['xy']-d['B2']
    spec=importlib.util.spec_from_file_location('audit_cmax',OUT/'frozen/spherical_cmax.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    p=mod.lift(relative[:,0],relative[:,1],0,0,60)
    result=mod.estimate(p,(d['time_us']-chosen['start_us'])*1e-6,max_rpm=3500,initial=None)
    assert result==chosen['result'],'Recorder changed estimator output'
    velocity=np.array([120.,-70.]);shift=(d['time_us']-d['time_us'][0])[:,None]*1e-6*velocity
    recovery=(d['xy']+shift)-shift
    equivariance_error=float(np.max(np.abs(recovery-d['xy'])))
    assert equivariance_error<1e-10
    report=dict(status='passed',paired_windows_verified=len(checks),max_center_reconstruction_error_px=max_center_error,
        recorder_replay_window=chosen['window'],recorder_replay_exact_match=True,
        synthetic_coordinate_shift_recovery_px=equivariance_error,protected_files=len(env['protected_hashes']),
        result_rows=len(details),checks=checks,
        caveat='Verifies computation, event identity and immutability, not physical center or spin accuracy.')
    (OUT/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
