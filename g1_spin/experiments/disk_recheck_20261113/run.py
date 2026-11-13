"""Fixed geometry/time/search-bound ablation for the original eccentric disk."""
import importlib.util,json,csv,hashlib
from pathlib import Path
from time import perf_counter
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];DATA=ROOT/'g1_spin/data/decoded/local_analysis_222403'
s=importlib.util.spec_from_file_location('cmax',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py');cmax=importlib.util.module_from_spec(s);s.loader.exec_module(cmax)
MODELS=ROOT/'g1_spin/experiments/real_center_sensitivity_20261102/models.json'
SEGMENTS={'early_turn':(1.6,1.9),'middle':(6.2,6.5),'late':(11.2,11.5)}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def center(model,t):return np.array(model['center_anchor_px'])+((np.asarray(t)/1e6)-model['anchor_s'])[...,None]*np.array(model['velocity_px_s'])

def main():
    if (OUT/'results.jsonl').exists():raise FileExistsError('Preserve output')
    manifest=json.loads((DATA/'manifest.json').read_text());paths=[DATA/'manifest.json']+[DATA/f for f in manifest['parts']]+[MODELS,OUT/'run.py',OUT/'protocol.md',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths};save('freeze.json',hashes)
    e=np.concatenate([np.load(DATA/f,allow_pickle=False) for f in manifest['parts']]);origin=int(e['t'][0]);t=e['t'].astype(np.int64)-origin
    assert not np.any(np.diff(t)<0);models=json.loads(MODELS.read_text());(OUT/'inputs').mkdir(exist_ok=True);rows=[]
    with (OUT/'results.jsonl').open('x') as stream:
        for segment,(start,end) in SEGMENTS.items():
            for offset in [.05,.15,.25]:
                mid=round((start+offset)*1e6)
                for ms in [20,50,100]:
                    lo,hi=mid-ms*500,mid+ms*500;a,b=np.searchsorted(t,[lo,hi]);raw=e[a:b];times=t[a:b];pix=np.column_stack([raw['x'],raw['y']]).astype(float)
                    cs={'fixed50':np.tile(center(models[segment]['50'],mid),(len(raw),1)),
                        'linear50':center(models[segment]['50'],times),'linear100':center(models[segment]['100'],times)}
                    rr=[48.]
                    # All main comparisons share support including the radius sensitivity on middle 50ms windows.
                    if offset==.15 and ms==50:rr=[46.,48.,50.]
                    mask=(pix[:,0]>=475)&(pix[:,0]<638)&(pix[:,1]>=143)&(pix[:,1]<240)
                    for c in cs.values():mask&=np.linalg.norm(pix-c,axis=1)<=.8*min(rr)
                    ix=np.flatnonzero(mask);common=ix+a
                    if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
                    name=f'{segment}_{mid:08d}_{ms}ms';np.savez_compressed(OUT/'inputs'/f'{name}.npz',events=raw[ix],ids=ix+a,common_ids=common,origin_us=origin,**{k:v[ix] for k,v in cs.items()})
                    jobs=[(m,48.,600) for m in cs]
                    if offset==.15 and ms==50:jobs += [(m,48.,3500) for m in cs]+[('linear50',r,600) for r in [46.,50.]]
                    for method,radius,bound in jobs:
                        rel=pix[ix]-cs[method][ix];begin=perf_counter();result=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,radius),(times[ix]-lo)*1e-6,max_rpm=bound)
                        row=dict(window=name,segment=segment,mid_us=mid,window_ms=ms,method=method,radius=radius,max_rpm=bound,raw_count=len(raw),common_count=len(common),input_count=len(ix),elapsed_s=perf_counter()-begin,result=result)
                        rows.append(row);stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                    print(name,len(ix),flush=True)
    assert len(rows)==96 and all(digest(ROOT/p)==h for p,h in hashes.items());save('run_check.json',dict(rows=96,windows=27,source_hashes_unchanged=True))

if __name__=='__main__':main()
