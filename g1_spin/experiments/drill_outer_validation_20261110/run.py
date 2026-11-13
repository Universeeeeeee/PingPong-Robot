"""End-to-end offline candidate/reference frontend comparison, fixed parameters."""
import importlib.util
import json
import sys
from pathlib import Path
from time import perf_counter
import numpy as np
from policy import choose

OUT=Path(__file__).resolve().parent
SRC=OUT.parent/'drill_outer_contour_20261110'
sys.path.insert(0,str(SRC))
s=importlib.util.spec_from_file_location('outer_front',SRC/'run.py');front=importlib.util.module_from_spec(s);s.loader.exec_module(front)
base=front.base;old=front.old;p2=base.p2;cmax=base.prev.cmax
BLOCKS={1:[.7,2.1,3.5,5.6],3:[1.4,2.6,4.2,6.2],4:[.8,1.6,2.4,3.2],5:[4.,4.4,5.6,6.4]}
RADIUS={1:60.,3:59.,4:60.,5:60.}


def main():
    if (OUT/'geometry.json').exists():raise FileExistsError('Preserve output')
    jobs=[dict(run=run,end_us=round((b+.205)*1e6),condition='translation' if run==5 else 'rotation') for run,blocks in BLOCKS.items() for b in blocks]
    jobs += [dict(run=5,end_us=905000,condition='static'),dict(run=5,end_us=2905000,condition='static')]
    previous=json.loads((front.SOURCE/'measurements.json').read_text())
    for j in jobs:
        assert not any(x['run']==j['run'] and max(j['end_us']-20000,x['mid_us']-x['duration_us']//2)<min(j['end_us'],x['mid_us']+x['duration_us']//2) for x in previous)
    old.save(OUT/'windows.json',jobs)
    paths=[OUT/'run.py',OUT/'protocol.md',OUT/'policy.py',OUT/'windows.json',SRC/'contour.py',front.SOURCE/'run.py',
           p2.OUT/'run.py',p2.OUT/'geometry.json',p2.OLD/'rotation_geometry.csv',p2.OLD/'control_geometry.csv',old.OUT/'run.py',
           base.ROOT/'g1_spin/src/g1_spin/spherical_cmax.py',old.OUT/'reference.json']
    hashes={str(p.relative_to(base.ROOT)):old.digest(p) for p in paths};old.save(OUT/'freeze.json',hashes)
    (OUT/'geometry_inputs').mkdir(exist_ok=True);(OUT/'spin_inputs').mkdir(exist_ok=True)
    geometries=[];routes=[]
    for run in BLOCKS:
        folder=p2.ARCHIVE/'decoded'/p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']];origin=int(parts[0]['t'][0]);ots,oxy=p2.observations(run)
        for j in [x for x in jobs if x['run']==run]:
            end=j['end_us'];name=f'run{run}_{end:08d}';static=j['condition']=='static';radius=RADIUS[run];theta=np.deg2rad(35 if run==3 else 90);direction=[np.cos(theta),np.sin(theta)]
            try:anchor=np.array([714.,355.]) if static else p2.centers(ots,oxy,[end-10000])[0]
            except ValueError:anchor=None
            local=[]
            if anchor is not None:
                for mid in range(end-50000,end+30001,10000):
                    raw,ids=old.read_window(parts,origin,mid-2500,mid+2500);pix=np.column_stack([raw['x'],raw['y']]);mask=np.linalg.norm(pix-anchor,axis=1)<=radius+190
                    raw=raw[mask];ids=ids[mask];candidate=front.extract(raw,anchor,radius,direction);points=np.array(candidate['points']).reshape(-1,2)
                    fit=base.fit_circle(points,anchor,radius,12);record=dict(case=name,run=run,mid_us=mid,anchor=anchor.tolist(),radius=radius,candidate=candidate,fit=fit)
                    path=OUT/'geometry_inputs'/f'{name}_{mid:08d}.npz';np.savez_compressed(path,events=raw,ids=ids,origin_us=origin)
                    record['input_file']=str(path.relative_to(OUT));record['input_sha256']=old.digest(path);local.append(record);geometries.append(record)
            good=[r for r in local if r['fit']['valid']];nts=[r['mid_us'] for r in good]
            a,b=choose(ots,nts,end,static)
            shaft=dict(valid=False,reason='no_anchor')
            if anchor is not None and run in [1,3,4]:
                raw,_=old.read_window(parts,origin,end-20000,end)
                shaft=base.prev.prev.shaft(raw,anchor,radius,direction)
            routes.append(dict(**j,case=name,anchor=None if anchor is None else anchor.tolist(),outer_valid=len(good),old_route=a,new_route=b,shaft=shaft))
        del parts
    old.save(OUT/'geometry.json',geometries);old.save(OUT/'routes.json',routes)
    old.save(OUT/'routing_freeze.json',{f:old.digest(OUT/f) for f in ['geometry.json','routes.json']})
    results=[]
    for run in BLOCKS:
        folder=p2.ARCHIVE/'decoded'/p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']];origin=int(parts[0]['t'][0]);ots,oxy=p2.observations(run)
        for j in [x for x in routes if x['run']==run]:
            end=j['end_us'];name=j['case'];radius=RADIUS[run];cache={}
            for lane in ['old','new']:
                route=j[lane+'_route'];cached=route in cache
                if not cached:
                    if route=='unavailable':value=dict(result=dict(valid=False,spin_rpm=None,axis_camera=None,omega_rad_s=None,failure_reason='geometry_unavailable'),elapsed_s=None,input_file=None)
                    else:
                        begin=end-(10000 if route=='old10' else 20000)
                        raw,ids=old.read_window(parts,origin,begin,end);t=raw['t'].astype(np.int64)-origin;pixels=np.column_stack([raw['x'],raw['y']]).astype(float)
                        if route=='old10':cs=np.tile(p2.centers(ots,oxy,[end-5000])[0],(len(raw),1))
                        elif route=='old20' and j['condition']=='static':cs=np.tile([714.,355.],(len(raw),1))
                        elif route=='old20':cs=p2.centers(ots,old.smooth(ots,oxy),t)
                        else:
                            gg=[g for g in geometries if g['case']==name and g['fit']['valid']];ts=np.array([g['mid_us'] for g in gg]);xy=np.array([g['fit']['center'] for g in gg]);cs=p2.centers(ts,old.smooth(ts,xy),t)
                        ix=np.flatnonzero(np.linalg.norm(pixels-cs,axis=1)<=.8*radius);common=ids[ix]
                        if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
                        path=OUT/'spin_inputs'/f'{name}_{route}.npz';np.savez_compressed(path,events=raw[ix],ids=ids[ix],common_ids=common,centers=cs[ix],origin_us=origin,radius=radius)
                        start=perf_counter();rel=pixels[ix]-cs[ix];result=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,radius),(t[ix]-begin)*1e-6,max_rpm=3500)
                        value=dict(result=result,elapsed_s=perf_counter()-start,input_file=str(path.relative_to(OUT)),input_sha256=old.digest(path),event_count=len(ix))
                    cache[route]=value
                results.append(dict(**j,lane=lane,route=route,cached_same_route=cached,**cache[route]))
            print(name,j['old_route'],j['new_route'],flush=True)
        del parts
    old.save(OUT/'results.json',results)
    assert len(results)==36 and all(old.digest(base.ROOT/p)==h for p,h in hashes.items())
    print('18 windows / 36 states complete; sources unchanged')

if __name__=='__main__':main()
