"""Same-support time regression and pixel-normal regression; no 3D mixed-unit TLS."""
import numpy as np

def representation(xy,times,ids,kind):
    if kind=='R_all':return np.arange(len(times))
    if kind!='R_latest':raise ValueError(kind)
    order=np.lexsort((ids,times,xy[:,1],xy[:,0]));pixels=xy[order]
    last=np.r_[np.any(pixels[1:]!=pixels[:-1],axis=1),True]
    return np.sort(order[last]) if len(order) else np.empty(0,int)

def canonical(normal,speed,intercept):
    sign=1.
    if normal[0]<-1e-12 or (abs(normal[0])<=1e-12 and normal[1]<0):sign=-1.
    return normal*sign,float(speed*sign),float(intercept*sign)

def fit(xy,times,kind):
    xy=np.asarray(xy,float);times=np.asarray(times,float);center=xy.mean(axis=0) if len(xy) else np.zeros(2);t0=float(times.mean()) if len(times) else 0.
    result=dict(model=kind,events=len(times),unique_pixels=len(np.unique(xy,axis=0)),available=False,failure_reason=None,
        center=center.tolist(),reference_s=t0,normal=None,speed_px_s=None,intercept_px=None,old_time_RMS_us=None,old_time_gate_us=None,
        normal_residual_RMS_px=None,Ctt_s2=None,eigenvalues_px2=None,eigengap_px2=None,condition=None)
    if len(times)<3:result['failure_reason']='insufficient_points';return result
    X=xy-center;dt=times-t0;Cxx=X.T@X/len(X);Cxt=X.T@dt/len(X);Ctt=float(dt@dt/len(X));result['Ctt_s2']=Ctt
    if Ctt<=1e-12:result['failure_reason']='insufficient_time_variance';return result
    if kind=='M_time':
        if np.linalg.matrix_rank(np.column_stack([X,np.ones(len(X))]))<3:result['failure_reason']='time_model_spatial_design_rank_deficient';return result
        matrix=np.column_stack([X,np.ones(len(X))]);coeff,_,rank,singular=np.linalg.lstsq(matrix,dt,rcond=None);gradient=coeff[:2];norm=float(np.linalg.norm(gradient))
        rms=float(np.sqrt(np.mean((matrix@coeff-dt)**2)));result['old_time_RMS_us']=rms*1e6;result['old_time_gate_us']=max(20e-6,.15*np.sqrt(Ctt))*1e6
        result['condition']=float(singular[0]/singular[-1]);result['gradient_s_px']=gradient.tolist()
        if norm<1e-10:result['failure_reason']='time_gradient_zero_unidentifiable_normal';return result
        normal=gradient/norm;speed=1/norm;intercept=-coeff[2]/norm
    elif kind=='M_spatial':
        matrix=Cxx-np.outer(Cxt,Cxt)/Ctt;eigenvalues,eigenvectors=np.linalg.eigh(matrix);gap=float(eigenvalues[1]-eigenvalues[0]);result['eigenvalues_px2']=eigenvalues.tolist();result['eigengap_px2']=gap
        if gap<=1e-8*max(1.,float(np.trace(Cxx))):result['failure_reason']='normal_direction_not_unique';return result
        normal=eigenvectors[:,0];speed=float(normal@Cxt/Ctt);intercept=0.
    else:raise ValueError(kind)
    normal,speed,intercept=canonical(normal,speed,intercept)
    residual=X@normal-speed*dt-intercept
    result.update(available=True,normal=normal.tolist(),speed_px_s=speed,intercept_px=intercept,normal_residual_RMS_px=float(np.sqrt(np.mean(residual*residual))),
        local_support_at_least12_pixels=result['unique_pixels']>=12)
    return result

def residual(model,xy,times):
    return abs((np.asarray(xy)-model['center'])@np.array(model['normal'])-model['speed_px_s']*(np.asarray(times)-model['reference_s'])-model['intercept_px'])

def equal_unit_distances(xy,times,polarity,distance):
    bins=np.floor(np.asarray(times)/.0005+1e-10).astype(int);keys=np.column_stack([xy.astype(int),polarity,bins]);unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    count=np.bincount(inverse);sums=np.bincount(inverse,weights=distance);means=sums/count
    unique_pixels=len(np.unique(xy,axis=0))
    return dict(units=len(unique),unique_pixels=unique_pixels,time_bins=len(np.unique(bins)),mean=float(np.mean(means)) if len(means) else None,
        median=float(np.median(means)) if len(means) else None,p95=float(np.percentile(means,95)) if len(means) else None,
        max=float(np.max(means)) if len(means) else None,unit_keys=unique.tolist(),unit_mean_distance_px=means.tolist(),unit_events=count.tolist(),
        evaluable=bool(unique_pixels>=4 and len(np.unique(bins))>=2))
