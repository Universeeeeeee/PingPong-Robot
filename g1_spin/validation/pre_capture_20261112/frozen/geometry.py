import numpy as np
from scipy.optimize import least_squares

def fit_circle(points,anchor,radius,minimum):
    if len(points)<minimum:return dict(valid=False,reason='insufficient_bins',bins=len(points))
    fit=least_squares(lambda c:np.linalg.norm(points-c,axis=1)-radius,anchor,bounds=(anchor-10,anchor+10),loss='soft_l1',f_scale=1.)
    delta=fit.x-points;dist=np.linalg.norm(delta,axis=1);jac=delta/dist[:,None]
    singular=np.linalg.svd(jac,compute_uv=False);condition=float(singular[0]/singular[-1])
    rms=float(np.sqrt(np.mean((dist-radius)**2)))
    reason='bound' if np.any(fit.active_mask) else 'residual' if rms>2.5 else 'condition' if condition>10 else None
    return dict(valid=reason is None,reason=reason,bins=len(points),center=fit.x.tolist(),rms_px=rms,condition=condition)
