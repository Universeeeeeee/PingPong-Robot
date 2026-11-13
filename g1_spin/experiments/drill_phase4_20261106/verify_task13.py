"""Independent same-event chain checks; do not generate new scenes or change estimates."""
import ast
import json
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial.transform import Rotation
from common import OUT,ROOT,P3,flow,module,digest,write_json

def objective(points,times,omega):
    dt=times-(times[0]+times[-1])/2;warped=Rotation.from_rotvec(-dt[:,None]*np.asarray(omega)).apply(points)
    lat=np.arcsin(np.clip(warped[:,2],-1,1));lon=np.arctan2(warped[:,1],warped[:,0]);y=(lat/np.pi+.5)*63;x=((lon/(2*np.pi)+.5)*128)%128
    iy=np.floor(y).astype(int);ix=np.floor(x).astype(int);fy=y-iy;fx=x-ix;hist=np.zeros((64,128))
    for dy in [0,1]:
        for dx in [0,1]:np.add.at(hist,(np.minimum(iy+dy,63),(ix+dx)%128),(fy if dy else 1-fy)*(fx if dx else 1-fx))
    return float(np.var(gaussian_filter(hist,1,truncate=2,mode=('reflect','wrap'))/len(points)))

def independent_ls(points,gradients,radius,method,rhs=None):
    matrix=[];target=[]
    for i,(point,gradient) in enumerate(zip(points,gradients)):
        x,y,z=point;a,b=gradient
        if method=='S':matrix.append(radius*np.cross(point,[a,b,0.]));target.append(1. if rhs is None else rhs[i])
        else:
            u,v=gradient/np.dot(gradient,gradient);velocity=np.array([u,v,-(x*u+y*v)/z])/radius
            P=np.outer(velocity,velocity)/np.dot(velocity,velocity);cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
            matrix.extend(P@(-cross));target.extend(P@velocity)
    matrix=np.array(matrix).reshape(-1,3);target=np.array(target)
    return np.linalg.lstsq(matrix,target,rcond=None)[0]

def main():
    frozen=json.loads((OUT/'synthetic_freeze.json').read_text())
    for rel,h in frozen['files'].items():assert digest(OUT/rel)==h,rel
    for rel,h in json.loads((OUT/'source_hashes.json').read_text()).items():assert digest(ROOT/rel)==h,rel
    # All manifest-covered stage3 outputs and their prior protected sources remain unchanged.
    manifest=json.loads((P3/'deliverable_manifest.json').read_text());prior=json.loads((P3/'source_hashes.json').read_text())
    for rel,value in manifest['files'].items():assert digest(P3/rel)==value['sha256'],rel
    for rel,h in prior.items():assert digest(ROOT/rel)==h,rel
    result=json.loads((OUT/'synthetic_results.json').read_text());checks=[];max_plane_error=0.;max_score_error=0.;score_count=0
    for scene in ['base','reverse','second','zero','translation','mixed']:
        record=next(r for r in result if r['scene']==scene);data=np.load(OUT/'inputs'/f'{scene}.npz');geometry=json.loads((OUT/'inputs'/f'{scene}_geometry.json').read_text())
        assert data['x'].dtype.kind in 'iu' and data['y'].dtype.kind in 'iu' and data['t_us'].dtype==np.int64
        assert np.all(np.diff(data['t_us'])>=0) and np.isin(data['polarity'],[0,1]).all()
        plane=np.load(OUT/'runner'/f'{scene}_planes.npz');train=plane['training_event_ids'];hold=plane['holdout_event_ids'];times=(data['t_us']-10000)*1e-6
        xy=np.column_stack([data['x'],data['y']]).astype(float);center=np.array(geometry['center_px']);velocity=np.array(geometry['center_velocity_px_s']);radius=geometry['radius_px']
        xy-=((data['t_us']-15000)*1e-6)[:,None]*velocity
        assert np.all((data['x'][train].astype(int)+2*data['y'][train].astype(int))%5!=0)
        assert np.all((data['x'][hold].astype(int)+2*data['y'][hold].astype(int))%5==0)
        expected=(xy[train]-center)/radius;points=np.column_stack([expected,-np.sqrt(1-np.sum(expected*expected,axis=1))]);assert np.all(points[:,2]<0) and np.allclose(np.linalg.norm(points,axis=1),1)
        for k,seed in enumerate(plane['seed_indices']):
            ids=plane['neighbor_indices'][plane['neighbor_offsets'][k]:plane['neighbor_offsets'][k+1]]
            local=xy[train[ids]];centroid=local.mean(axis=0);X=np.column_stack([local-centroid,np.ones(len(ids))]);dt=times[train[ids]]-times[train[seed]]
            # Positive-definite normal equations independently check the accepted full-rank local fit.
            coeff=np.linalg.solve(X.T@X,X.T@dt);err=float(np.max(abs(coeff-plane['coefficients'][k])));max_plane_error=max(max_plane_error,err);assert err<1e-12
            assert np.all(data['polarity'][train[ids]]==data['polarity'][train[seed]])
            assert np.all(np.sum((local-xy[train[seed]])**2,axis=1)<=16+1e-10)
            assert np.all(abs(times[train[ids]]-times[train[seed]])<=.002+1e-15)
            assert not np.intersect1d(train[ids],hold).size
        for route in record['routes']:
            init=route['initializer'];search=route['search'];method='A' if route['method']=='A_literal' else route['method']
            if init is not None and len(plane['points']):
                reconstructed=independent_ls(plane['points'],plane['gradients'],radius,method)
                assert np.allclose(reconstructed,init['raw_omega_rad_s'],atol=1e-10,rtol=1e-11)
            if search is not None and search['optimizer_candidate_omega'] is not None:
                candidates=[dict(omega=search['optimizer_candidate_omega'],contrast=search['optimizer_candidate_contrast'])]
                if method!='G':
                    assert search['actual_minimize_calls']==1 and not search['G_fallback_called'];candidates+=search['magnitude_scan']
                    direction=np.array(init['raw_omega_rad_s']);direction/=np.linalg.norm(direction)
                    assert candidates[1]['rpm']==0 and np.allclose(candidates[1]['omega'],0)
                    for scan in search['magnitude_scan']:
                        assert np.linalg.norm(np.cross(scan['omega'],direction))<1e-8 and np.dot(scan['omega'],direction)>=-1e-10
                for candidate in candidates:
                    calculated=objective(points,times[train],candidate['omega']);error=abs(calculated-candidate['contrast']);max_score_error=max(max_score_error,error);score_count+=1
                    assert np.isclose(calculated,candidate['contrast'],atol=1e-15,rtol=1e-10)
            elif method!='G':assert not init['available']
        checks.append(dict(scene=scene,integer_events=len(xy),training_events=len(train),holdout_events=len(hold),planes=len(plane['points']),truth_not_read_by_verifier_runner_input_reconstruction=True))
    # AST and file-access code boundary: runner has no evaluator or generator imports and no truth paths.
    tree=ast.parse((OUT/'run_synthetic.py').read_text());imports=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.ImportFrom):imports.append(node.module)
        elif isinstance(node,ast.Import):imports.extend(n.name for n in node.names)
    assert not set(imports)&{'generate','evaluate'}
    opened_truth=[n for n in ast.walk(tree) if isinstance(n,ast.Constant) and isinstance(n.value,str) and n.value.startswith('truth/')]
    assert not opened_truth
    eval_result=json.loads((OUT/'synthetic_evaluation.json').read_text());assert eval_result['status']=='passed' and len(eval_result['S_three_full_chain_gates'])==3
    assert json.loads((OUT/'semantic_verification.json').read_text())['status']=='passed'
    verification=dict(status='passed',scenes=6,existing_frozen_algorithms_unchanged=True,stage3_manifest_files_unchanged=len(manifest['files']),
        prior_protected_sources_unchanged=len(prior),same_integer_events_to_plane_LS_and_search_verified=True,
        runner_truth_access='AST imports and paths checked; evaluator truth never opened by runner',
        accepted_local_planes_independently_rebuilt=True,max_local_coefficient_error=max_plane_error,
        independent_IWE_score_checks=score_count,max_IWE_score_difference=max_score_error,three_S_full_chain_gates_passed=True,
        checks=checks,scope='Task1+3 synthetic checkpoint; no real front-end experiments or new real-data accuracy claims.')
    write_json(OUT/'verification_task13.json',verification);print(json.dumps(verification,indent=2))

if __name__=='__main__':main()
