"""Truth-isolated runner: events+geometry only, frozen stage3 algorithms."""
import json
from time import perf_counter
import numpy as np
from common import OUT,flow,module,digest,write_json,array_hash

methods=module('stage4_frozen_methods',OUT/'frozen/stage3_methods.py')

def prepared(scene):
    data=np.load(OUT/'inputs'/f'{scene}.npz');geometry=json.loads((OUT/'inputs'/f'{scene}_geometry.json').read_text())
    ids=np.arange(len(data['x']));hold=(data['x'].astype(int)+2*data['y'].astype(int))%5==0
    selected=ids[~hold]
    if len(selected)>5000:selected=selected[np.linspace(0,len(selected)-1,5000).astype(int)]
    xy=np.column_stack([data['x'],data['y']]).astype(float)
    times=(data['t_us']-geometry['window_start_us'])*1e-6;dt=(data['t_us']-geometry['reference_us'])*1e-6
    velocity=np.array(geometry['center_velocity_px_s']);compensated=xy-dt[:,None]*velocity
    return data,geometry,selected,ids[hold],xy,times,compensated

def search_counted(points,times,initial):
    actual=0;original=methods.minimize;original_g=methods.cmax.estimate
    def count(*args,**kwargs):
        nonlocal actual
        actual+=1;return original(*args,**kwargs)
    def no_g(*args,**kwargs):raise AssertionError('No G fallback permitted')
    methods.minimize=count;methods.cmax.estimate=no_g
    try:search=methods.axis_search(points,times,initial)
    finally:methods.minimize=original;methods.cmax.estimate=original_g
    assert actual==1;search['actual_minimize_calls']=actual;search['G_fallback_called']=False
    return search

def raw_rhs_init(planes,times,geometry):
    gradients=planes['gradients'];velocity=np.array(geometry['center_velocity_px_s']);radius=geometry['radius_px']
    center=np.array(geometry['center_px']);reference=(geometry['reference_us']-geometry['window_start_us'])*1e-6
    t=times[planes['seed_indices']]+planes['coefficients'][:,2]
    xy=planes['pixel_centers']-center-(t-reference)[:,None]*velocity
    s=methods.cmax.lift(xy[:,0],xy[:,1],0,0,radius)
    matrix=radius*np.cross(s,np.column_stack([gradients,np.zeros(len(s))]));rhs=1-gradients@velocity
    if len(s):
        omega,_,rank,singular=np.linalg.lstsq(matrix,rhs,rcond=None)
        condition=float(singular[0]/singular[-1]) if len(singular)==3 and singular[-1]>1e-15 else None
        residual=float(np.sqrt(np.mean((matrix@omega-rhs)**2)))
    else:omega=np.zeros(3);rank=0;singular=[];condition=None;residual=None
    status='available'
    if len(s)<12:status='insufficient_flow_constraints'
    elif rank<3 or condition is None or condition>100:status='degenerate_flow_system'
    elif not np.isfinite(omega).all() or np.linalg.norm(omega)<1e-6:status='nonfinite_or_zero_initial_direction'
    return dict(available=status=='available',failure_reason=None if status=='available' else status,
        planes=len(s),rank=int(rank),condition=condition,raw_omega_rad_s=omega.tolist(),residual_rms=residual,
        rhs_min=float(rhs.min()) if len(rhs) else None,rhs_max=float(rhs.max()) if len(rhs) else None,
        point_reference='Each raw plane coordinate centroid and its fitted centroid timestamp, with supplied center trajectory',
        formula='r*(s cross [gx,gy,0]) dot omega = 1 - gradient dot center_velocity')

def main():
    if (OUT/'synthetic_results.json').exists():raise FileExistsError('Preserve synthetic results')
    (OUT/'runner').mkdir(exist_ok=True);records=[]
    for scene in ['base','reverse','second','zero','translation','mixed']:
        data,geometry,ids,hold,xy,times,comp=prepared(scene);center=np.array(geometry['center_px']);radius=geometry['radius_px']
        input_hash=array_hash(ids,comp[ids],times[ids],data['polarity'][ids],center,np.array(radius))
        relative=comp[ids]-center;points=methods.cmax.lift(relative[:,0],relative[:,1],0,0,radius)
        start=perf_counter();planes,meta=flow.extract(comp[ids],times[ids],data['polarity'][ids],center,radius);extraction_s=perf_counter()-start
        np.savez_compressed(OUT/'runner'/f'{scene}_planes.npz',**planes,training_event_ids=ids,holdout_event_ids=hold)
        # The complete event extraction is a fixed support diagnostic, not replacement inputs.
        complete,complete_meta=flow.extract(comp,times,data['polarity'],center,radius)
        full_initializers={m:flow.initialize(complete['points'],complete['gradients'],radius,m) for m in ['A','S']}
        rawplanes,rawmeta=flow.extract(xy[ids],times[ids],data['polarity'][ids],center,radius)
        rawinit=raw_rhs_init(rawplanes,times[ids],geometry)
        np.savez_compressed(OUT/'runner'/f'{scene}_raw_planes.npz',**rawplanes,training_event_ids=ids)
        routes=[]
        for method in ['A','S','G']:
            init=None;search=None;outcome=None
            if method=='G':
                search=methods.g_search(points,times[ids]);outcome=methods.decision(points,times[ids],search,'G')
            else:
                init=flow.initialize(planes['points'],planes['gradients'],radius,method)
                if init['available'] and methods.precheck(points,times[ids]) is None:
                    search=search_counted(points,times[ids],init['raw_omega_rad_s']);outcome=methods.decision(points,times[ids],search,method)
            routes.append(dict(method='A_literal' if method=='A' else method,input_hash=input_hash,initializer=init,search=search,decision=outcome,
                terminal_reason=(None if outcome['measurement_valid'] else outcome['common_quality']['failure_reason']) if outcome is not None else init['failure_reason']))
        record=dict(scene=scene,events=len(data['x']),training_events=len(ids),holdout_events=len(hold),input_hash=input_hash,
            input_event_file_sha256=digest(OUT/'inputs'/f'{scene}.npz'),geometry_file_sha256=digest(OUT/'inputs'/f'{scene}_geometry.json'),
            extraction_metadata=meta,extraction_s=extraction_s,complete_event_extraction_metadata=complete_meta,
            complete_event_initializers=full_initializers,raw_extraction_metadata=rawmeta,raw_rhs_S_initializer=rawinit,
            planes_sha256=digest(OUT/'runner'/f'{scene}_planes.npz'),routes=routes,
            truth_access='Runner never opens truth/ or generation configuration; geometry center velocity is explicit supplied geometric input.')
        records.append(record);print(scene,len(data['x']),'train',len(ids),'planes',meta['accepted_planes'],
            [(r['method'],r['terminal_reason']) for r in routes],flush=True)
    # Old planar fixture remains a standalone rank-degeneracy control.
    i,j=np.meshgrid(np.arange(-6,7),np.arange(-6,7));pixel=np.column_stack([640+i.ravel(),360+j.ravel()]).astype(float)
    ticks=10000+100*(i.ravel()+6)+150*(j.ravel()+6);order=np.argsort(ticks,kind='stable');pixel=pixel[order];ticks=ticks[order]
    planes,meta=flow.extract(pixel,(ticks-ticks[0])*1e-6,np.ones(len(ticks),np.uint8),[640,360],60)
    deg=dict(events=len(ticks),metadata=meta,initializers={m:flow.initialize(planes['points'],planes['gradients'],60,m) for m in ['A','S']})
    write_json(OUT/'synthetic_results.json',records);write_json(OUT/'rank_degeneracy_control.json',deg)

if __name__=='__main__':main()
