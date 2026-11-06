"""T-only fitting boundary and separately evaluated H/V predictions."""
import json
from time import perf_counter
import numpy as np
from helpers import flow,methods,hash_arrays
from frontends import extract
from prediction import regions,in_box,roi_from_source,normal_support,local_holdout,future,compare_future

def Tfit(Traw,Ttime,Tpol,Tids,center,velocity,radius,frontend):
    start=perf_counter();comp=Traw-(Ttime-.005)[:,None]*velocity
    planes,metadata=extract(comp,Ttime,Tpol,Traw.astype(int),Tids,center,radius,frontend);flow_s=perf_counter()-start
    relative=comp-center;points=methods.cmax.lift(relative[:,0],relative[:,1],0,0,radius)
    return planes,metadata,points,flow_s

def S_raw_diagnostic(Traw,Ttime,Tpol,Tids,center,velocity,radius,frontend):
    raw,metadata=extract(Traw,Ttime,Tpol,Traw.astype(int),Tids,center,radius,frontend)
    t=Ttime[raw['seed_indices']]+raw['coefficients'][:,2]
    local=raw['pixel_centers']-center-(t-.005)[:,None]*velocity
    points=methods.cmax.lift(local[:,0],local[:,1],0,0,radius)
    gradient=raw['gradients'];matrix=radius*np.cross(points,np.column_stack([gradient,np.zeros(len(gradient))]));rhs=1-gradient@velocity
    if len(points):
        omega,_,rank,singular=np.linalg.lstsq(matrix,rhs,rcond=None);condition=float(singular[0]/singular[-1]) if len(singular)==3 and singular[-1]>1e-15 else None
        residual=float(np.sqrt(np.mean((matrix@omega-rhs)**2)))
    else:omega=np.zeros(3);rank=0;condition=None;residual=None
    available=len(points)>=12 and rank==3 and condition is not None and condition<=100 and np.linalg.norm(omega)>=1e-6
    return dict(planes=len(points),available=bool(available),rank=int(rank),condition=condition,raw_omega_rad_s=omega.tolist(),residual_rms=residual,
        rhs_min=float(rhs.min()) if len(rhs) else None,rhs_max=float(rhs.max()) if len(rhs) else None,
        formula='r*(s cross [gx,gy,0]) dot omega = 1-gradient dot velocity',native_scalar_residual_pass=bool(residual is not None and residual<=.35),
        status='semantic_diagnostic_only_not_frontend_selection',metadata=metadata),raw

def route_search(points,times,planes,radius,method):
    start=perf_counter();init=flow.initialize(planes['points'],planes['gradients'],radius,method);ls_s=perf_counter()-start
    failure=methods.precheck(points,times)
    result=dict(initializer=init,search=None,decision=None,ls_s=ls_s,search_s=0.,status=failure or init['failure_reason'])
    if failure or not init['available']:return result
    original=methods.minimize;original_g=methods.cmax.estimate;count=0
    def counted(*args,**kwargs):
        nonlocal count
        count+=1;return original(*args,**kwargs)
    def forbidden(*args,**kwargs):raise AssertionError('Frontend A/S cannot fall back to G')
    methods.minimize=counted;methods.cmax.estimate=forbidden
    try:search=methods.axis_search(points,times,init['raw_omega_rad_s'])
    finally:methods.minimize=original;methods.cmax.estimate=original_g
    assert count==1;search['actual_NM_calls']=count;search['G_fallback_called']=False
    decision=methods.decision(points,times,search,method)
    result.update(search=search,decision=decision,search_s=search['search_s'],status=decision['common_quality']['failure_reason'])
    return result

def split(raw_xy,times):
    hold=(np.floor(raw_xy[:,0]/2).astype(int)+np.floor(raw_xy[:,1]/2).astype(int))%2==0
    return dict(T=(times<.006)&~hold,H=(times<.006)&hold,V=times>=.006)

def evaluate_route(omega,bundle,Tmask,source_mask,center,velocity,radius,roi):
    raw=bundle['xy'];times=bundle['times'];pol=bundle['polarity'];ids=bundle['ids'];v=times>=.006
    src=Tmask&source_mask
    return future(raw[src],times[src],pol[src],ids[src],raw[v],times[v],pol[v],ids[v],center,velocity,radius,omega,roi)

def accepted_prediction(fitted,comparison,translation_comparison):
    return bool(fitted['decision'] and fitted['decision']['measurement_valid'] and comparison and comparison['engineering_pass'] and translation_comparison and translation_comparison['engineering_pass'])

def process_case(bundle,lane,frontend,include_LOO=True):
    raw=bundle['xy'];times=bundle['times'];pol=bundle['polarity'];ids=bundle['ids'];center=bundle['center'];radius=bundle['radius'];velocity=np.array(lane['velocity'],float)
    masks=split(raw,times);T=masks['T'];H=masks['H'];V=masks['V'];selected=regions(raw[T],center)
    Traw=raw[T];Ttime=times[T];Tpol=pol[T];Tids=ids[T]
    planes,meta,points,flow_s=Tfit(Traw,Ttime,Tpol,Tids,center,velocity,radius,frontend)
    support=normal_support(planes,Traw,selected);compH=raw[H]-(times[H]-.005)[:,None]*velocity;compT=Traw-(Ttime-.005)[:,None]*velocity
    holdout=local_holdout(planes,compT,Ttime,Tpol,compH,times[H],pol[H],ids[H])
    H_region_counts=[]
    for region in selected:
        count=int(np.count_nonzero(in_box(raw[H],region['bbox'])&np.isin(ids[H],holdout['predicted_H_ids'])))
        H_region_counts.append(dict(region=region['region'],unique_predicted_H=count))
    holdout['region_support']=H_region_counts;holdout['regions_with_at_least5_H']=sum(r['unique_predicted_H']>=5 for r in H_region_counts)
    holdout['engineering_pass']=holdout['engineering_pass'] and holdout['regions_with_at_least5_H']>=3
    source=(times>=.005)&(times<.006);roi=roi_from_source(raw[T]);zero=evaluate_route([0,0,0],bundle,T,source,center,np.zeros(2),radius,roi)
    translation=evaluate_route([0,0,0],bundle,T,source,center,velocity,radius,roi)
    routes={};allLOO=[]
    for name,method in [('A_literal','A'),('S_comp','S')]:
        fitted=route_search(points,Ttime,planes,radius,method);future_result=None;comparison=None;translation_comparison=None;loo=[]
        if fitted['search'] is not None:
            omega=fitted['search']['optimizer_candidate_omega']
            future_result=evaluate_route(omega,bundle,T,source,center,velocity,radius,roi);comparison=compare_future(future_result,zero)
            translation_comparison=compare_future(future_result,translation)
            if include_LOO:
                for region in selected:
                    remove=in_box(Traw,region['bbox'],4);keep=~remove
                    p,m,s,cost=Tfit(Traw[keep],Ttime[keep],Tpol[keep],Tids[keep],center,velocity,radius,frontend)
                    localfit=route_search(s,Ttime[keep],p,radius,method)
                    held=T&in_box(raw,region['bbox']);held_source=held&source;localroi=roi_from_source(raw[held])
                    base=evaluate_route([0,0,0],bundle,held,source,center,np.zeros(2),radius,localroi)
                    trans=evaluate_route([0,0,0],bundle,held,source,center,velocity,radius,localroi)
                    prediction=None;cmp=None;transcmp=None
                    if localfit['search'] is not None:
                        prediction=evaluate_route(localfit['search']['optimizer_candidate_omega'],bundle,held,source,center,velocity,radius,localroi)
                        cmp=compare_future(prediction,base)
                        transcmp=compare_future(prediction,trans)
                    entry=dict(region=region,removed_T_ids=Tids[remove].tolist(),remaining_T_ids=Tids[keep].tolist(),
                        extraction_metadata={k:v for k,v in m.items() if k!='seed_diagnostics'},fitted=localfit,flow_s=cost,
                        evaluator_only_source_ids=ids[held_source].tolist(),fixed_roi=localroi,zero_future=base,translation_future=trans,model_future=prediction,comparison=cmp,model_vs_translation_only=transcmp,
                        accepted_and_predictively_supported=accepted_prediction(localfit,cmp,transcmp))
                    loo.append(entry)
        fitted.update(future_V=future_result,comparison_V=comparison,model_vs_translation_only=translation_comparison,leave_region_out=loo,
            pipeline_s=flow_s+fitted['ls_s']+fitted['search_s'],flow_s=flow_s)
        decision=fitted['decision'];usable=bool(decision and decision['measurement_valid'])
        evaluable=sum(bool(e['comparison'] and e['comparison']['evaluable']) for e in loo)
        passed=sum(e['accepted_and_predictively_supported'] for e in loo)
        fitted['qualification']=dict(global_search_accepted=usable,region_support=support['sufficient'],local_H=holdout['engineering_pass'],
            future_V=bool(comparison and comparison['engineering_pass'] and translation_comparison and translation_comparison['engineering_pass']),LOO_evaluable=evaluable,LOO_passed=passed,
            pass_rotation=bool(usable and support['sufficient'] and holdout['engineering_pass'] and comparison and comparison['engineering_pass'] and translation_comparison and translation_comparison['engineering_pass'] and passed>=3))
        routes[name]=fitted
    rawdiagnostic,rawplanes=S_raw_diagnostic(Traw,Ttime,Tpol,Tids,center,velocity,radius,frontend)
    # Physical ID maps and each surface/cluster support are retained for independent leak audit.
    result=dict(case=bundle['case'],source=bundle['source'],frontend=frontend,lane=lane,counts={k:int(v.sum()) for k,v in masks.items()},
        T_ids=ids[T].tolist(),H_ids=ids[H].tolist(),V_ids=ids[V].tolist(),coordinate_input_hash=hash_arrays(Tids,compT,Ttime,Tpol,center,np.array(radius)),
        regions=selected,normal_support=support,extraction_metadata=meta,holdout_H=holdout,
        zero_future=zero,translation_only_future=translation,translation_vs_zero=compare_future(translation,zero),routes=routes,
        S_raw_RHS=rawdiagnostic,geometry_caveat='Velocity is a frozen geometry hypothesis, not independent physical truth.')
    arrays=dict(**planes,T_physical_ids=Tids,H_physical_ids=ids[H],V_physical_ids=ids[V])
    return result,arrays,rawplanes
