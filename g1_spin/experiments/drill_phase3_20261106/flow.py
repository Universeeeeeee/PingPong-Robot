"""Shared frozen-neighborhood planes and two explicitly different initializers."""
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree

def skew(points):
    points=np.asarray(points,float);x,y,z=points.T;zero=np.zeros(len(points))
    return np.stack([np.stack([zero,-z,y],axis=1),np.stack([z,zero,-x],axis=1),np.stack([-y,x,zero],axis=1)],axis=1)

def lift_velocity(points,flow_pixels_s,radius):
    x,y,z=np.asarray(points,float).T;fx,fy=np.asarray(flow_pixels_s,float).T
    return np.column_stack([fx/radius,fy/radius,-(x*fx+y*fy)/(radius*z)])

def extract(xy,times,polarity,center,radius):
    """Same 4px/2ms/300seed filters as frozen optical_flow; preserve each accepted fit."""
    xy=np.asarray(xy,float);times=np.asarray(times,float);polarity=np.asarray(polarity)
    empty=dict(points=np.empty((0,3)),gradients=np.empty((0,2)),pixel_centers=np.empty((0,2)),seed_indices=np.empty(0,np.int64),
        coefficients=np.empty((0,3)),rms_s=np.empty(0),time_std_s=np.empty(0),spatial_ratio=np.empty(0),
        neighbor_offsets=np.zeros(1,np.int64),neighbor_indices=np.empty(0,np.int64))
    if len(xy)<100:return empty,dict(accepted_planes=0,seeds_attempted=0,failures={'insufficient_events':1})
    tree=cKDTree(xy);seeds=np.linspace(0,len(xy)-1,min(300,len(xy)),dtype=int)
    pts=[];grads=[];centers=[];accepted=[];coefficients=[];rmses=[];stds=[];ratios=[];neighbors=[];offsets=[0];fail=Counter()
    for i in seeds:
        ix=np.asarray(tree.query_ball_point(xy[i],4),dtype=int)
        ix=ix[(np.abs(times[ix]-times[i])<=.002)&(polarity[ix]==polarity[i])]
        if len(ix)<10:fail['too_few_neighbors']+=1;continue
        coords=xy[ix];c=coords.mean(axis=0);centered=coords-c
        singular=np.linalg.svd(centered,compute_uv=False)
        if singular[0]<1e-6 or singular[1]/singular[0]<.08:fail['spatial_plane_degenerate']+=1;continue
        dt=times[ix]-times[i];std=float(np.std(dt))
        if std<1e-6:fail['insufficient_local_time_span']+=1;continue
        design=np.column_stack([centered,np.ones(len(ix))]);coeff,_,rank,_=np.linalg.lstsq(design,dt,rcond=None)
        rms=float(np.sqrt(np.mean((design@coeff-dt)**2)))
        if rank<3 or rms>max(20e-6,.15*std):fail['plane_residual_or_rank']+=1;continue
        sx,sy=(c-center)/radius
        if sx*sx+sy*sy>=.9**2:fail['outside_flow_inner_sphere']+=1;continue
        if np.linalg.norm(coeff[:2])<1e-8:fail['vanishing_gradient']+=1;continue
        pts.append([sx,sy,-np.sqrt(1-sx*sx-sy*sy)]);grads.append(coeff[:2]);centers.append(c);accepted.append(i)
        coefficients.append(coeff);rmses.append(rms);stds.append(std);ratios.append(float(singular[1]/singular[0]));neighbors.extend(ix.tolist());offsets.append(len(neighbors))
    if not pts:return empty,dict(accepted_planes=0,seeds_attempted=len(seeds),failures=dict(fail))
    data=dict(points=np.array(pts),gradients=np.array(grads),pixel_centers=np.array(centers),seed_indices=np.array(accepted,np.int64),
        coefficients=np.array(coefficients),rms_s=np.array(rmses),time_std_s=np.array(stds),spatial_ratio=np.array(ratios),
        neighbor_offsets=np.array(offsets,np.int64),neighbor_indices=np.array(neighbors,np.int64))
    return data,dict(accepted_planes=len(pts),seeds_attempted=len(seeds),failures=dict(fail))

def system(points,gradients,radius,method):
    points=np.asarray(points,float);gradients=np.asarray(gradients,float)
    if method=='S':return radius*np.cross(points,np.column_stack([gradients,np.zeros(len(points))])),np.ones(len(points)),None
    f=gradients/np.sum(gradients*gradients,axis=1)[:,None]
    velocity=lift_velocity(points,f,radius)
    projection=velocity[:,:,None]*velocity[:,None,:]/np.sum(velocity*velocity,axis=1)[:,None,None]
    return (projection@(-skew(points))).reshape(-1,3),(projection@velocity[:,:,None]).reshape(-1),projection

def initialize(points,gradients,radius,method):
    n=len(points)
    if n:
        matrix,rhs,_=system(points,gradients,radius,method)
        omega,_,rank,singular=np.linalg.lstsq(matrix,rhs,rcond=None)
        condition=float(singular[0]/singular[-1]) if len(singular)==3 and singular[-1]>1e-15 else None
        residual=float(np.sqrt(np.mean((matrix@omega-rhs)**2)))
    else:omega=np.zeros(3);rank=0;singular=np.empty(0);condition=None;residual=None
    failure=None
    if n<12:failure='insufficient_flow_constraints'
    elif rank<3 or condition is None or condition>100:failure='degenerate_flow_system'
    elif not np.isfinite(omega).all() or np.linalg.norm(omega)<1e-6:failure='nonfinite_or_zero_initial_direction'
    native=failure
    if method=='S' and not native and residual>.35:native='inconsistent_flow_constraints'
    return dict(method=method,available=failure is None,failure_reason=failure,planes=n,rank=int(rank),condition=condition,
        singular_values=singular.tolist(),raw_omega_rad_s=omega.tolist(),residual_rms=residual,
        residual_units='dimensionless' if method=='S' else 'unit_sphere_per_second',
        native_scalar_status=('accepted' if native is None else native) if method=='S' else 'not_applicable',
        caveat='A/S share rank/count/condition/direction gates; no residual acceptance gate except separately reported scalar native status.')
