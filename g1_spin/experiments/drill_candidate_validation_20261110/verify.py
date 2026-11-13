"""Independent raw indexing, common-mask reconstruction and baseline eligibility."""
import importlib.util
import json
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('candidate_trial',OUT/'run.py');trial=importlib.util.module_from_spec(s);s.loader.exec_module(trial)
rows=[json.loads(s) for s in (OUT/'results.jsonl').read_text().splitlines()]
jobs=json.loads((OUT/'windows.json').read_text());assert len(rows)==88
assert len({(r['window'],r['method']) for r in rows})==88
checks=[];eligibility=[]
for run in trial.BLOCKS:
    folder=trial.p2.ARCHIVE/'decoded'/trial.p2.NAMES[run]
    m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']]
    offsets=np.cumsum([0]+[len(p) for p in parts]);origin=int(parts[0]['t'][0]);ts,xy=trial.p2.observations(run)
    sm=xy.copy()
    for i in range(1,len(ts)-1):
        if np.all(np.diff(ts[i-1:i+2])<=30000):sm[i]=np.sort(xy[i-1:i+2],axis=0)[1]
    for j in [j for j in jobs if j['run']==run]:
        begin,end=j['begin_us'],j['end_us'];name=f'run{run}_{end:08d}';rr=[r for r in rows if r['window']==name]
        raw=[];ids=[]
        for i,p in enumerate(parts):
            mask=(p['t']>=origin+begin)&(p['t']<origin+end)
            raw.append(p[mask]);ids.append(np.flatnonzero(mask)+offsets[i])
        raw=np.concatenate(raw);ids=np.concatenate(ids);t=raw['t'].astype(np.int64)-origin
        pix=np.column_stack([raw['x'],raw['y']]);radius=trial.RADIUS[run]
        if j['condition']=='static':fixed=np.array([714.,355.]);centers=np.tile(fixed,(len(t),1))
        else:
            fixed=np.array([np.interp(end-5000,ts,xy[:,d]) for d in range(2)])
            centers=np.column_stack([np.interp(t,ts,sm[:,d]) for d in range(2)])
        common=(np.sum((pix-fixed)**2,axis=1)<=(.8*radius)**2)&(np.sum((pix-centers)**2,axis=1)<=(.8*radius)**2)
        for r in rr:
            path=OUT/'inputs'/f"{name}_{r['method']}.npz"
            if not path.exists():continue
            mask=common&(t>=end-(10000 if r['method']=='baseline' else 20000))
            selected=np.flatnonzero(mask)
            assert len(selected)==r['common_count']
            if len(selected)>5000:selected=selected[np.linspace(0,len(selected)-1,5000).astype(int)]
            data=np.load(path)
            assert np.array_equal(ids[selected],data['ids'])
            assert np.array_equal(raw[selected],data['events'])
            expected=np.tile(fixed,(len(selected),1)) if r['method']=='baseline' else centers[selected]
            assert np.allclose(expected,data['centers'],rtol=0,atol=1e-9)
            if r['result']['valid']:
                omega=r['result']['omega_rad_s'];assert np.isclose(np.linalg.norm(omega)*60/(2*np.pi),r['result']['spin_rpm'])
            checks.append(dict(window=name,method=r['method'],event_count=len(selected),raw_and_mask_and_center_verified=True))
        # Check standalone 10ms fixed center could be prepared when joint coverage failed.
        baseline=next(r for r in rr if r['method']=='baseline')
        if baseline['result']['failure_reason'] in ['geometry_gap_over_30ms','geometry_outside_coverage']:
            queries=[end-10000,end-5000,end]
            can=True
            for q in queries:
                k=np.searchsorted(ts,q)
                if k==len(ts) or q<ts[0] or (ts[k]!=q and ts[k]-ts[k-1]>30000):can=False
            eligibility.append(dict(window=name,baseline_10ms_geometry_available=can))
for p,h in json.loads((OUT/'freeze.json').read_text()).items():assert trial.old.digest(trial.ROOT/p)==h
result=dict(result_rows=88,prepared_inputs=len(checks),protected_hashes_unchanged=True,checks=checks)
trial.old.save(OUT/'verification.json',result);trial.old.save(OUT/'baseline_eligibility.json',eligibility)
print(json.dumps(dict(result_rows=88,prepared_inputs=len(checks),baseline_eligibility=eligibility)))
# Supplement uses standalone 10ms geometry, so verify it outside the shared-mask lane.
supp_path=OUT/'baseline_supplement.json'
if supp_path.exists():
    supplemental=[]
    for row in json.loads(supp_path.read_text()):
        data=np.load(OUT/'inputs'/f"{row['window']}_standalone_baseline.npz")
        run=row['run'];folder=trial.p2.ARCHIVE/'decoded'/trial.p2.NAMES[run]
        m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']]
        off=0;actual=[];allids=[];allraw=[];origin=int(parts[0]['t'][0])
        for part in parts:
            ids=data['ids'];local=ids[(ids>=off)&(ids<off+len(part))]-off;actual.append(part[local])
            mask=(part['t']>=origin+row['begin_us'])&(part['t']<origin+row['end_us'])
            allraw.append(part[mask]);allids.append(np.flatnonzero(mask)+off);off+=len(part)
        assert np.array_equal(np.concatenate(actual),data['events'])
        raw=np.concatenate(allraw);ids=np.concatenate(allids);ts,xy=trial.p2.observations(run)
        center=np.array([np.interp(row['end_us']-5000,ts,xy[:,d]) for d in range(2)])
        pix=np.column_stack([raw['x'],raw['y']]);sel=np.flatnonzero(np.sum((pix-center)**2,axis=1)<=(.8*trial.RADIUS[run])**2)
        if len(sel)>5000:sel=sel[np.linspace(0,len(sel)-1,5000).astype(int)]
        assert np.array_equal(ids[sel],data['ids']) and np.allclose(data['centers'],center,rtol=0,atol=1e-9)
        supplemental.append(dict(window=row['window'],raw_mask_center_verified=True))
    trial.old.save(OUT/'supplement_verification.json',supplemental)
    print('Supplement inputs independently verified',len(supplemental))
