"""Reconstruct physical event IDs independently of experiment window extraction."""
import hashlib
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
NAMES = {1:'201058',4:'203829',5:'215030'}
rows = [json.loads(s) for s in (OUT/'results.jsonl').read_text().splitlines()]
assert len(rows) == 126
assert len({(r['window'],r['method']) for r in rows}) == 126
checks=[]
for run, suffix in NAMES.items():
    folder=ROOT/f'g1_spin/data/archive_20261105/decoded/DSP03F100025_20261105_{suffix}'
    manifest=json.loads((folder/'manifest.json').read_text())
    parts=[np.load(folder/f,mmap_mode='r') for f in manifest['parts']]
    offsets=np.cumsum([0]+[len(p) for p in parts])
    origin=int(parts[0]['t'][0])
    for path in sorted((OUT/'inputs').glob(f'run{run}_*.npz')):
        d=np.load(path); ids=d['ids']; common=d['common_ids']
        group=[r for r in rows if r['window']==path.stem]
        assert len(group)==3
        r=group[0]; begin,end=r['begin_us'],r['end_us']
        actual=[]; raw_pixels=[];raw_ids=[]
        for i,p in enumerate(parts):
            loc=ids[(ids>=offsets[i])&(ids<offsets[i+1])]-offsets[i]
            actual.append(p[loc])
            mask=(p['t']>=origin+begin)&(p['t']<origin+end)
            raw=p[mask];raw_pixels.append(np.column_stack([raw['x'],raw['y']]))
            raw_ids.append(np.flatnonzero(mask)+offsets[i])
        actual=np.concatenate(actual)
        assert np.array_equal(actual,d['events'])
        allids=np.concatenate(raw_ids)
        assert len(allids)==r['raw_count']
        assert np.isin(common,allids).all()
        expected=common if len(common)<=5000 else common[np.linspace(0,len(common)-1,5000).astype(int)]
        assert np.array_equal(ids,expected)
        assert len(ids)==r['optimizer_count'] and len(common)==r['common_count']
        for m in ['fixed','linear','median3']:
            rel=np.column_stack([actual['x'],actual['y']])-d[m]
            q=np.sum((rel/60)**2,axis=1)
            assert np.all(q<=.64+1e-12)
        for item in group:
            if item['valid']:
                assert np.isclose(np.linalg.norm(item['omega'])*60/(2*np.pi),item['rpm'])
            else:
                assert item['rpm'] is None and item['omega'] is None
        checks.append(dict(window=path.stem,events=len(ids),raw_ids_verified=True))
for filename,h in json.loads((OUT/'freeze.json').read_text()).items():
    assert hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==h
result=dict(rows=len(rows),unique_windows=len({r['window'] for r in rows}),
            prepared_windows=len(checks),checks=checks,protected_hashes_unchanged=True)
(OUT/'independent_verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
