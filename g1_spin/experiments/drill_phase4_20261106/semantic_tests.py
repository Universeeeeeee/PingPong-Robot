"""A_literal retained; covector dual and raw translation RHS are analytic comparisons."""
import numpy as np
from common import OUT,flow,write_json

def dual_system(points,gradients,radius,rhs):
    lifted_covector=radius*np.column_stack([gradients,np.zeros(len(points))])
    h=lifted_covector-points*np.sum(points*lifted_covector,axis=1)[:,None]
    h2=np.sum(h*h,axis=1);P=h[:,:,None]*h[:,None,:]/h2[:,None,None]
    velocity=rhs[:,None]*h/h2[:,None]
    return (P@(-flow.skew(points))).reshape(-1,3),velocity.reshape(-1),h,P

def main():
    xy=np.array([[.6,0],[0,.6],[-.5,.2],[.2,-.5],[.4,.4],[-.3,-.4]])
    points=np.column_stack([xy,-np.sqrt(1-np.sum(xy*xy,axis=1))]);truth=np.array([1.,2.,3.]);radius=120.
    normals=np.array([[1.,1.],[1.,-.3],[.2,1.],[.7,-1.],[1.,.1],[.4,1.]])
    cdot=np.array([400.,-300.]);vel=radius*np.cross(truth,points)[:,:2]+cdot
    gradients=normals/np.sum(normals*vel,axis=1)[:,None];rhs=1-gradients@cdot
    scalar=radius*np.cross(points,np.column_stack([gradients,np.zeros(len(points))]))
    matrix,b,h,P=dual_system(points,gradients,radius,rhs)
    assert np.allclose(scalar@truth,rhs,atol=1e-12)
    assert np.allclose(np.sum(h*np.cross(truth,points),axis=1),rhs,atol=1e-12)
    assert np.allclose(matrix@truth,b,atol=1e-12)
    assert np.allclose(P,P.transpose(0,2,1)) and np.allclose(P@P,P)
    # Perturbed RHS demonstrates that stacking P induces the explicit 1/||h|| weights.
    noisy=rhs+np.array([.02,-.03,.01,.04,-.02,.05]);dual,d_rhs,_,_=dual_system(points,gradients,radius,noisy)
    dual_w=np.linalg.lstsq(dual,d_rhs,rcond=None)[0];weight=1/np.linalg.norm(h,axis=1)
    weighted=np.linalg.lstsq(scalar*weight[:,None],noisy*weight,rcond=None)[0]
    unweighted=np.linalg.lstsq(scalar,noisy,rcond=None)[0]
    assert np.allclose(dual_w,weighted,atol=1e-12) and not np.allclose(weighted,unweighted,atol=1e-5)
    # Exact translation-only plane, b=0: keep its nonzero covector/projector, zero tangent velocity.
    g=np.array([[1/500.,-1/1500.]]);s=points[[0]];zero_rhs=1-g@cdot
    assert np.allclose(zero_rhs,0,atol=1e-12)
    zero_m,zero_b,zero_h,zero_P=dual_system(s,g,radius,zero_rhs)
    assert np.linalg.norm(zero_h)>0 and np.isfinite(zero_P).all() and np.allclose(zero_b,0)
    # f_normal - translation is generally not a valid residual normal-flow vector.
    f=g/np.sum(g*g,axis=1)[:,None];bad=f-cdot
    assert np.linalg.norm(bad)>1 and np.allclose(np.sum(g*bad,axis=1),0,atol=1e-12)
    # Whole-coordinate transformation changes the time-plane gradient by dividing by b.
    physical=np.array([.005,.004]);translation=np.array([30.,-20.]);bcoef=1-physical@translation
    corrected=physical/bcoef;sample_xy=np.array([[1.,2.],[3.,-1.],[-2.,4.],[5.,3.]])
    t=sample_xy@physical+.002;coords=sample_xy-t[:,None]*translation
    assert np.allclose(coords@corrected+.002/bcoef,t)
    # Original stage3 non-equivalence counterexample remains unchanged.
    s=np.array([[.6,0,-.8]]);original_g=np.array([[-5.,-5.]]);omega=np.array([0.,1.,1.])
    a,ab,projection=flow.system(s,original_g,1.,'A');sm,sb,_=flow.system(s,original_g,1.,'S')
    mismatch=(projection@np.cross(omega,s)[:,:,None]).reshape(-1)-ab
    assert np.allclose(sm@omega,1) and np.isclose(np.linalg.norm(mismatch),.24597367497455452)
    result=dict(status='passed',tests=['raw_gradient_translation_RHS','covector_tangent_scalar_per_constraint',
        'projection_stack_equals_explicit_weighted_not_unweighted_LS','zero_RHS_keeps_projector_without_dividing_by_b',
        'normal_flow_minus_translation_is_not_full_flow','coordinate_then_refit_exact_plane_identity','stage3_A_literal_counterexample_retained'],
        raw_RHS_true_omega=truth.tolist(),raw_scalar_recovery=np.linalg.lstsq(scalar,rhs,rcond=None)[0].tolist(),
        dual_LS=dual_w.tolist(),weighted_scalar_LS=weighted.tolist(),unweighted_scalar_LS=unweighted.tolist(),
        zero_RHS=zero_rhs.tolist(),zero_RHS_h=zero_h.tolist(),zero_RHS_P=zero_P.tolist(),bad_normal_subtraction=bad.tolist(),
        preserved_A_literal_mismatch_norm=float(np.linalg.norm(mismatch)),
        scope='Analytic semantics only; no fourth real-data algorithm or changes to stage3 A_literal.')
    write_json(OUT/'semantic_verification.json',result);print(result)

if __name__=='__main__':main()
