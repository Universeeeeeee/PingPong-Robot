"""Frozen G and explicit directed-axis A/S searches, with separate quality semantics."""
import importlib.util
from time import perf_counter
import numpy as np
from scipy.optimize import minimize
from common import P1,OUT

# Before freeze the exact phase1 source supports analytic tests; after freeze use local frozen copy.
source=OUT/'frozen/spherical_cmax.py'
if not source.exists():source=P1/'frozen/spherical_cmax.py'
spec=importlib.util.spec_from_file_location('stage3_frozen_cmax',source)
cmax=importlib.util.module_from_spec(spec);spec.loader.exec_module(cmax)
BOUND=3500*2*np.pi/60

def score(points,dt,omega):return float(np.var(cmax.iwe(cmax.warp(points,dt,omega))/len(points)))

def precheck(points,times):
    if len(points)<100:return 'insufficient_events'
    if not np.isfinite(points).all() or not np.isfinite(times).all():return 'nonfinite_geometry_or_times'
    if np.any(np.diff(times)<0):return 'nonmonotonic_timestamps'
    if np.ptp(times)<=1e-6:return 'insufficient_time_span'
    if np.max(abs(np.linalg.norm(points,axis=1)-1))>1e-5:return 'nonunit_sphere'
    return None

def common_quality(points,times,omega,contrast,zero,converged):
    pre=precheck(points,times)
    if pre:return dict(passed=False,failure_reason=pre,gain=None)
    gain=contrast/max(zero,1e-12)-1;speed=float(np.linalg.norm(omega))
    reason=None
    if speed*np.ptp(times)<.03 or gain<.08:reason='no_rotation_evidence'
    elif not converged:reason='optimization_not_converged'
    elif speed>=BOUND*.985:reason='search_boundary'
    elif np.linalg.eigvalsh(np.cov(points.T))[-2]<1e-5:reason='spatially_degenerate_features'
    return dict(passed=reason is None,failure_reason=reason,gain=float(gain),displacement_rad=float(speed*np.ptp(times)))

def g_search(points,times):
    calls=[];objective_calls=0;contrast_calls=0;original=cmax.minimize;original_iwe=cmax.iwe
    def counted_iwe(*args,**kwargs):
        nonlocal contrast_calls
        contrast_calls+=1;return original_iwe(*args,**kwargs)
    def trace(fun,seed,**kwargs):
        local_calls=0
        def counted(omega):
            nonlocal objective_calls,local_calls
            objective_calls+=1;local_calls+=1;return fun(omega)
        result=original(counted,seed,**kwargs)
        calls.append(dict(seed=np.asarray(seed).tolist(),omega=result.x.tolist(),fun=float(result.fun),
            converged=bool(result.success),iterations=int(result.nit),objective_calls=local_calls))
        return result
    cmax.minimize=trace;cmax.iwe=counted_iwe;start=perf_counter()
    try:native=cmax.estimate(points,times,max_rpm=3500,initial=None)
    finally:cmax.minimize=original;cmax.iwe=original_iwe
    elapsed=perf_counter()-start
    candidate=native['omega_rad_s'] if native['valid'] else native['quality'].get('rejected_candidate_omega_rad_s')
    zero=native['quality'].get('contrast_before');contrast=native['quality'].get('contrast_after')
    return dict(native=native,optimizer_candidate_omega=candidate,optimizer_candidate_contrast=contrast,zero_contrast=zero,
        selected_seed=native['quality'].get('initial_omega_rad_s'),magnitude_scan=native['quality'].get('candidate_scan',[]),
        optimizer_converged=native['quality'].get('converged'),nm_calls=len(calls),objective_calls=objective_calls,contrast_calls=contrast_calls,
        local_solutions=calls,search_s=elapsed,multistart_check='performed_by_frozen_G')

def axis_search(points,times,initial):
    axis=np.asarray(initial,float);axis/=np.linalg.norm(axis);dt=times-(times[0]+times[-1])/2
    n_speed=int(np.clip(np.ceil(BOUND*np.ptp(times)/.4),8,30));scans=[];contrast_calls=0;objective_calls=0;history=[]
    def contrast(omega):
        nonlocal contrast_calls
        contrast_calls+=1;return score(points,dt,omega)
    start=perf_counter()
    for magnitude in np.r_[0.,np.linspace(BOUND/n_speed,BOUND,n_speed)]:
        scans.append(dict(rpm=float(magnitude*60/(2*np.pi)),omega=(axis*magnitude).tolist(),contrast=contrast(axis*magnitude)))
    zero=scans[0]['contrast'];best=max(scans,key=lambda s:s['contrast']);seed=np.array(best['omega'])
    def objective(omega):
        nonlocal objective_calls
        objective_calls+=1;speed=np.linalg.norm(omega)
        if speed>BOUND:return 1+(speed/BOUND-1)**2
        value=contrast(omega);history.append(value);return -value/max(zero,1e-12)
    simplex=np.vstack([seed,seed+np.eye(3)*max(2.,BOUND/40)])
    local=minimize(objective,seed,method='Nelder-Mead',options=dict(initial_simplex=simplex,maxiter=350,xatol=.05,fatol=1e-6))
    after=contrast(local.x);elapsed=perf_counter()-start
    return dict(native=None,optimizer_candidate_omega=local.x.tolist(),optimizer_candidate_contrast=after,zero_contrast=zero,
        selected_seed=seed.tolist(),magnitude_scan=scans,optimizer_converged=bool(local.success),nm_calls=1,
        objective_calls=objective_calls,contrast_calls=contrast_calls,local_solutions=[dict(seed=seed.tolist(),omega=local.x.tolist(),fun=float(local.fun),
            converged=bool(local.success),iterations=int(local.nit),objective_calls=objective_calls)],objective_history=history,
        search_s=elapsed,multistart_check='not_performed_single_start_not_a_pass')

def decision(points,times,search,method):
    candidate=search['optimizer_candidate_omega'];zero=search['zero_contrast'];contrast=search['optimizer_candidate_contrast']
    if candidate is None:
        pre=precheck(points,times)
        return dict(common_quality=dict(passed=False,failure_reason=pre or 'no_optimizer_candidate',gain=None),zero_won=None,
            candidate_omega=None,candidate_axis=None,candidate_rpm=None,candidate_contrast=None,measurement_valid=False,
            measurement_omega=None,measurement_axis=None,measurement_rpm=None)
    zero_won=zero>=contrast
    selected=np.zeros(3) if zero_won else np.asarray(candidate,float);selected_score=zero if zero_won else contrast
    quality=common_quality(points,times,selected,selected_score,zero,search['optimizer_converged'])
    # G keeps its extra native multistart rejection; A/S cannot count an absent rivalry check as passed.
    valid=quality['passed'] and (method!='G' or search['native']['valid'])
    speed=float(np.linalg.norm(selected));axis=None if speed<1e-12 else (selected/speed).tolist()
    return dict(common_quality=quality,zero_won=bool(zero_won),candidate_omega=selected.tolist(),candidate_axis=axis,
        candidate_rpm=speed*60/(2*np.pi),candidate_contrast=float(selected_score),measurement_valid=bool(valid),
        measurement_omega=selected.tolist() if valid else None,measurement_axis=axis if valid else None,measurement_rpm=speed*60/(2*np.pi) if valid else None)
