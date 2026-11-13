"""One offline entry for three frozen spin estimators. No hardware control."""
import argparse,json,platform
from pathlib import Path
from time import perf_counter
import numpy as np
from common import ROOT,Events,read_geometry,smooth,save,digest,verify_lock
from frozen import cmax,contour,geometry,policy


def evaluate(manifest_path,out):
    verify_lock();manifest_path=Path(manifest_path).resolve();m=json.loads(manifest_path.read_text())
    if m['schema_version']!=1 or m['geometry_frame']!='event_camera_pixels' or m['geometry_clock']!='event_relative_us':raise ValueError('Unsupported geometry frame/clock/schema')
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    event_path=(manifest_path.parent/m['events']).resolve();geom_path=(manifest_path.parent/m['geometry_csv']).resolve()
    events=Events(event_path,m['event_origin_us']);ts,xy=read_geometry(geom_path)
    if len(set(w['id'] for w in m['windows']))!=len(m['windows']):raise ValueError('Duplicate window ids')
    files=[manifest_path,geom_path]+events.paths
    input_hashes={str(p):digest(p) for p in set(files)}
    save(out/'provenance.json',dict(inputs=input_hashes,code_lock=json.loads((ROOT/'freeze.json').read_text()),python=platform.python_version(),numpy=np.__version__))
    results=[]
    for job in m['windows']:
        end=job['end_us'];radius=float(job['radius_px']);anchor=np.array(job['anchor_px'],float);direction=np.array(job['shaft_direction_xy'],float)
        if not isinstance(end,int) or end<50000 or radius<=0 or not np.isfinite(radius) or anchor.shape!=(2,) or direction.shape!=(2,) or not np.isfinite(anchor).all() or not np.isfinite(direction).all() or np.linalg.norm(direction)==0:raise ValueError('Invalid window geometry')
        if any(ch not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for ch in job['id']):raise ValueError('Invalid window id')
        static=job.get('fixed_center_prior',False)
        outer=[];frontend_start=perf_counter()
        for mid in range(end-50000,end+30001,10000):
            raw,_=events.window(mid-2500,mid+2500)
            candidate=contour.extract(raw,anchor,radius,direction)
            fit=geometry.fit_circle(np.array(candidate['points']).reshape(-1,2),anchor,radius,12)
            outer.append(dict(t_us=mid,fit=fit,selected_pixels=candidate['points']))
        frontend_s=perf_counter()-frontend_start
        good=[x for x in outer if x['fit']['valid']];nts=np.array([x['t_us'] for x in good],np.int64);nxy=np.array([x['fit']['center'] for x in good]).reshape(-1,2)
        old_route,new_route=policy.choose(ts,nts,end,static)
        baseline='old10' if static or policy.covers(ts,end-10000,end) else 'unavailable'
        routes={'baseline10':baseline,'smooth20_fallback':old_route,'outer20_fallback':new_route}
        save(out/(job['id']+'_geometry.json'),dict(routes=routes,outer=outer,frontend_s=frontend_s))
        cache={}
        for variant,route in routes.items():
            cached=route in cache
            if not cached:
                if route=='unavailable':value=dict(valid=False,spin_rpm=None,axis_camera=None,omega_rad_s=None,failure_reason='geometry_unavailable',compute_s=None,input_count=0)
                else:
                    begin=end-(10000 if route=='old10' else 20000);raw,ids=events.window(begin,end);t=raw['t'].astype(np.int64)-events.origin;pix=np.column_stack([raw['x'],raw['y']]).astype(float)
                    if static and route.startswith('old'):centers=np.tile(anchor,(len(raw),1))
                    elif route=='old10':centers=np.tile([np.interp(end-5000,ts,xy[:,d]) for d in range(2)],(len(raw),1))
                    else:
                        gt,gxy=(nts,nxy) if route=='outer20' else (ts,xy);gxy=smooth(gt,gxy)
                        centers=np.column_stack([np.interp(t,gt,gxy[:,d]) for d in range(2)])
                    ix=np.flatnonzero(np.linalg.norm(pix-centers,axis=1)<=.8*radius);common_count=len(ix)
                    if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
                    start=perf_counter();rel=pix[ix]-centers[ix];value=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,radius),(t[ix]-begin)*1e-6,max_rpm=3500)
                    value.update(compute_s=perf_counter()-start,input_count=len(ix),common_count=common_count)
                    np.savez_compressed(out/f"{job['id']}_{route}_input.npz",ids=ids[ix],events=raw[ix],centers=centers[ix])
                cache[route]=value
            reference=job.get('reference_rpm');metrics={}
            if reference is not None:
                if reference['kind'] not in ['same_setting_approximate','synchronized_reference','synthetic_truth']:raise ValueError('Unknown reference kind')
                metrics=dict(reference_kind=reference['kind'],reference_rpm=reference['value'],abs_reference_deviation_rpm=abs(value['spin_rpm']-reference['value']) if value['valid'] else None)
            results.append(dict(window=job['id'],condition=job['condition'],variant=variant,route=route,cached_same_route=cached,frontend_s=frontend_s if variant=='outer20_fallback' else 0.,**value,**metrics))
        print(job['id'],routes,flush=True)
    save(out/'results.json',results)
    summary=[]
    for condition in sorted(set(r['condition'] for r in results)):
        for variant in routes:
            rr=[r for r in results if r['condition']==condition and r['variant']==variant];valid=[r for r in rr if r['valid']]
            summary.append(dict(condition=condition,variant=variant,planned=len(rr),valid=len(valid),valid_fraction=len(valid)/len(rr),rpm_median=float(np.median([r['spin_rpm'] for r in valid])) if valid else None,failures={reason:sum(r['failure_reason']==reason for r in rr) for reason in sorted(set(r['failure_reason'] for r in rr if r['failure_reason']))}))
    save(out/'summary.json',summary)
    if any(digest(Path(p))!=sha for p,sha in input_hashes.items()):raise ValueError('Input changed during evaluation')
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('manifest');p.add_argument('--out',required=True);a=p.parse_args();evaluate(a.manifest,a.out)
