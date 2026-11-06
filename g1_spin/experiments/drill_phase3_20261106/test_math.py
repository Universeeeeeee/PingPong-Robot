"""Fixed analytic counterexamples: do not alter A to force agreement with truth."""
import json
import numpy as np
from common import OUT,write_json
from flow import skew,lift_velocity,system,initialize
from methods import common_quality

def main():
    tests=[]
    s=np.array([[.6,0.,-.8]]);omega=np.array([0.,1.,1.]);truev=np.cross(omega,s);g=np.array([[-5.,-5.]])
    f=g/np.sum(g*g,axis=1)[:,None];v=lift_velocity(s,f,1.);matrix,rhs,P=system(s,g,1.,'A');sm,sb,_=system(s,g,1.,'S')
    mismatch=(P@truev[:,:,None]).reshape(-1,3)-v
    assert np.allclose(truev,[[-.8,.6,-.6]]) and np.allclose(v,[[-.1,-.1,-.075]])
    assert np.allclose(sm@omega,[1.]);assert np.allclose(mismatch,[[-.15365853658536585,-.15365853658536585,-.1152439024390244]])
    assert np.isclose(np.linalg.norm(mismatch),.245973675,atol=1e-9)
    tests.append('fixed_single_point_counterexample_preserved')
    xy=np.array([[.6,0],[0,.6],[-.5,.2],[.2,-.5],[.4,.4],[-.3,-.4]])
    points=np.column_stack([xy,-np.sqrt(1-np.sum(xy*xy,axis=1))]);trueomega=np.array([1.,2.,3.])
    normals=np.array([[1.,1.],[1.,-.3],[.2,1.],[.7,-1.],[1.,.1],[.4,1.]])
    fullv=np.cross(trueomega,points);gradients=normals/np.sum(normals*fullv[:,:2],axis=1)[:,None]
    recovered={};ranks={}
    for method in ['A','S']:
        matrix,rhs,projection=system(points,gradients,1.,method);w,_,rank,_=np.linalg.lstsq(matrix,rhs,rcond=None)
        recovered[method]=w.tolist();ranks[method]=int(rank);assert rank==3
        if method=='S':assert np.allclose(w,trueomega,atol=1e-12)
        else:assert np.allclose(w,[1.10926831,1.81278968,3.75706699],atol=1e-8)
    tests.append('six_point_full_rank_S_exact_A_differs_as_specified')
    projection=fullv[:,:,None]*fullv[:,None,:]/np.sum(fullv*fullv,axis=1)[:,None,None]
    matrix=(projection@(-skew(points))).reshape(-1,3);rhs=(projection@fullv[:,:,None]).reshape(-1)
    full_rec=np.linalg.lstsq(matrix,rhs,rcond=None)[0];assert np.allclose(full_rec,trueomega,atol=1e-12)
    tests.append('full_tangential_velocity_projection_recovers_true_omega')
    for pro in [P,projection,system(points,gradients,1.,'A')[2]]:
        assert np.allclose(pro,pro.transpose(0,2,1));assert np.allclose(pro@pro,pro)
    tests.append('P_symmetric_and_idempotent')
    def sphere(pixel,r):
        xy=pixel/r;return np.r_[xy,-np.sqrt(1-np.dot(xy,xy))]
    radius=60.;pixel=np.array([36.,12.]);point=sphere(pixel,radius);pixelvelocity=np.array([120.,-70.]);epsilon=1e-7
    fd=(sphere(pixel+epsilon*pixelvelocity,radius)-sphere(pixel-epsilon*pixelvelocity,radius))/(2*epsilon)
    analytic=lift_velocity(point[None,:],pixelvelocity[None,:],radius)[0]
    assert np.allclose(fd,analytic,atol=1e-8);assert abs(np.dot(point,analytic))<1e-12
    J=np.vstack([np.eye(2),-point[:2]/point[2]])
    assert not np.allclose(J.T@J,np.eye(2));tests.append('sphere_lift_finite_difference_tangent_and_JtJ_not_identity')
    for method in ['A','S']:
        a,b,_=system(points,gradients,1.,method);aa,bb,_=system(points,gradients/10,10.,method)
        assert np.allclose(a,aa) and np.allclose(b,bb)
        neg_a,neg_b,_=system(points,-gradients,1.,method)
        assert np.allclose(np.linalg.lstsq(neg_a,neg_b,rcond=None)[0],-np.array(recovered[method]))
        deg=initialize(np.tile(points[0],(12,1)),np.tile(gradients[0],(12,1)),1.,method)
        assert not deg['available'] and deg['failure_reason']=='degenerate_flow_system'
        assert initialize(points,gradients,1.,method)['failure_reason']=='insufficient_flow_constraints'
        assert initialize(np.empty((0,3)),np.empty((0,2)),1.,method)['failure_reason']=='insufficient_flow_constraints'
    tests.append('units_sign_count_rank_and_empty_degeneracy')
    assert common_quality(np.tile(points,(20,1)),np.linspace(0,.01,120),np.zeros(3),1.,1.,True)['failure_reason']=='no_rotation_evidence'
    tests.append('zero_candidate_is_no_evidence_not_valid_measurement')
    report=dict(status='passed',tests=tests,single_point=dict(sphere=s.tolist(),trueomega=omega.tolist(),truev=truev.tolist(),g=g.tolist(),
        lifted_normal_velocity=v.tolist(),scalar_lhs=(sm@omega).tolist(),projection_mismatch=mismatch.tolist(),mismatch_norm=float(np.linalg.norm(mismatch))),
        six_point=dict(trueomega=trueomega.tolist(),A_omega=recovered['A'],S_omega=recovered['S'],ranks=ranks,full_velocity_P_recovery=full_rec.tolist()),
        finite_difference_max_error=float(np.max(abs(fd-analytic))),JtJ=(J.T@J).tolist(),
        interpretation='The literal normal-flow Euclidean 3D P interpretation differs from the scalar image-plane aperture system because J transpose J is not identity. This does not establish an error in undisclosed Sony implementation.')
    write_json(OUT/'math_verification.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
