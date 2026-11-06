"""Evaluator-only truth access: no truth values enter extraction, initialization or search."""
import json
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import OUT,flow,write_json,table

def angle(a,b):
    if min(np.linalg.norm(a),np.linalg.norm(b))<1e-12:return None
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b),-1,1))))

def error(omega,truth):
    if omega is None:return dict(axis_error_deg=None,speed_relative_error=None,omega_error_norm=None)
    omega=np.array(omega);speed=np.linalg.norm(truth)
    return dict(axis_error_deg=angle(omega,truth),speed_relative_error=None if speed<1e-12 else float(np.linalg.norm(omega)/speed-1),
        omega_error_norm=float(np.linalg.norm(omega-truth)))

def lift(pixels,radius):
    xy=pixels/radius;return np.column_stack([xy,-np.sqrt(1-np.sum(xy*xy,axis=1))])

def summary(values):
    return dict(count=len(values),median=float(np.median(values)) if len(values) else None,p95=float(np.percentile(values,95)) if len(values) else None,max=float(np.max(values)) if len(values) else None)

def heldout_motion(points,times,labels,train_ids,hold_ids,omega,radius):
    q=Rotation.from_rotvec(-(times-.005)[:,None]*np.array(omega)).apply(points)
    rows=[]
    for edge in np.unique(labels):
        training=train_ids[labels[train_ids]==edge];hold=hold_ids[labels[hold_ids]==edge]
        if len(training)<3 or not len(hold):continue
        _,_,vh=np.linalg.svd(q[training],full_matrices=False);normal=vh[-1]
        residual=radius*np.abs(q[hold]@normal)
        rows.append(dict(edge_id=int(edge),training_events=len(training),holdout_events=len(hold),normal_learned_from_training_only=normal.tolist(),
            residual_px=summary(residual)))
    return rows

def main():
    results=json.loads((OUT/'synthetic_results.json').read_text());evaluation=[];rowtable=[];loo=[];local_rows=[]
    fig,axes=plt.subplots(2,3,figsize=(15,10));gate_checks=[];max_root=0.;max_rounded=0.;polarity_mismatches=0
    for result,ax in zip(results,axes.flat):
        scene=result['scene'];truth=np.load(OUT/'truth'/f'{scene}.npz');data=np.load(OUT/'inputs'/f'{scene}.npz')
        geo=json.loads((OUT/'inputs'/f'{scene}_geometry.json').read_text());radius=geo['radius_px'];center=np.array(geo['center_px'])
        omega=truth['omega'];velocity=truth['velocity'];labels=truth['edge_ids'];exact=truth['exact_time_s'];times=(data['t_us']-10000)*1e-6
        xy=np.column_stack([data['x'],data['y']]).astype(float);comp=xy-(times-.005)[:,None]*velocity;points=lift(comp-center,radius)
        runner=np.load(OUT/'runner'/f'{scene}_planes.npz');train=runner['training_event_ids'];hold=runner['holdout_event_ids']
        assert not np.intersect1d(train,hold).size and np.array_equal(np.r_[np.sort(np.r_[train,hold])],np.arange(len(xy)))
        # Independently recompute exact crossing, visibility, material arc and derivative polarity.
        if len(xy):
            s=lift(xy-center-(exact-.005)[:,None]*velocity,radius)
            q=Rotation.from_rotvec(-(exact-.005)[:,None]*omega).apply(s);normals=truth['normals'][labels]
            residual=abs(np.sum(q*normals,axis=1));max_root=max(max_root,float(residual.max()));assert residual.max()<1e-10
            patch=np.arccos(np.clip(np.sum(q*truth['anchors'][labels],axis=1),-1,1));assert patch.max()<=.14+1e-10
            assert np.all(s[:,2]<0)
            ds=np.column_stack([np.full(len(s),-velocity[0]/radius),np.full(len(s),-velocity[1]/radius),
                (s[:,0]*velocity[0]+s[:,1]*velocity[1])/(radius*s[:,2])])
            dq=Rotation.from_rotvec(-(exact-.005)[:,None]*omega).apply(ds-np.cross(omega,s))
            polarity=(np.sum(normals*dq,axis=1)>0).astype(np.uint8);polarity_mismatches+=int(np.count_nonzero(polarity!=data['polarity']))
            max_rounded=max(max_rounded,float(np.max(abs(times-exact))));assert np.max(abs(times-exact))<=.5e-6+1e-15
        # Plane physical-edge mapping is evaluation metadata only, after all runner results exist.
        plane_edges=[];plane_purity=[]
        for k in range(len(runner['points'])):
            ix=runner['neighbor_indices'][runner['neighbor_offsets'][k]:runner['neighbor_offsets'][k+1]]
            counts=Counter(labels[train[ix]].tolist());edge,count=counts.most_common(1)[0];plane_edges.append(edge);plane_purity.append(count/len(ix))
        plane_edges=np.array(plane_edges);per_edge=dict(Counter(map(int,plane_edges)));assert not plane_purity or min(plane_purity)==1.
        # Raw heldout events were not fit or optimized; validate local time prediction without labels.
        hold_xy=comp[hold];hold_time=times[hold];local_rms=[];local_ratios=[]
        tree=cKDTree(hold_xy)
        for k,pixel in enumerate(runner['pixel_centers']):
            seed=runner['seed_indices'][k];event=train[seed];ii=np.array(tree.query_ball_point(pixel,4),int)
            ii=ii[(abs(hold_time[ii]-times[event])<=.002)&(data['polarity'][hold[ii]]==data['polarity'][event])]
            if not len(ii):continue
            predicted=(hold_xy[ii]-pixel)@runner['gradients'][k]+runner['coefficients'][k,2]
            delta=hold_time[ii]-times[event];rms=float(np.sqrt(np.mean((predicted-delta)**2)));local_rms.append(rms*1e6)
            scale=max(20e-6,.15*runner['time_std_s'][k]);local_ratios.append(rms/scale)
            local_rows.append(dict(scene=scene,plane=k,physical_edge_eval_only=int(plane_edges[k]),holdout_events=len(ii),rms_us=rms*1e6,
                training_residual_gate_us=scale*1e6,holdout_rms_over_training_gate=rms/scale))
        routes=[]
        for route in result['routes']:
            init=route['initializer'];final=route['decision'];initial=None if init is None else init['raw_omega_rad_s']
            measurement=None if final is None else final['measurement_omega'];initial_error=error(initial,omega);final_error=error(measurement,omega)
            hold_metrics=None if measurement is None else heldout_motion(points,times,labels,train,hold,measurement,radius)
            init_metrics=None if initial is None else heldout_motion(points,times,labels,train,hold,initial,radius)
            info=dict(method=route['method'],initial_error=initial_error,final_error=final_error,
                final_heldout_motion=hold_metrics,initial_heldout_motion=init_metrics,terminal_reason=route['terminal_reason'])
            routes.append(info);rowtable.append(dict(scene=scene,method=route['method'],events=len(xy),planes=len(runner['points']),
                physical_edges=len(per_edge),rank=None if init is None else init['rank'],condition=None if init is None else init['condition'],
                init_axis_error_deg=initial_error['axis_error_deg'],init_speed_relative_error=initial_error['speed_relative_error'],
                final_axis_error_deg=final_error['axis_error_deg'],final_speed_relative_error=final_error['speed_relative_error'],
                measurement_valid=bool(final and final['measurement_valid']),terminal_reason=route['terminal_reason'],
                final_heldout_edge_residual_median_px=summary([p['residual_px']['median'] for p in hold_metrics])['median'] if hold_metrics else None))
            if scene in ['base','reverse','second'] and route['method']=='S':
                passed=initial_error['axis_error_deg']<=10 and abs(initial_error['speed_relative_error'])<=.15 and final_error['axis_error_deg']<=5 and abs(final_error['speed_relative_error'])<=.1
                assert passed and final['measurement_valid'];gate_checks.append(dict(scene=scene,S_full_chain_passed=True,initial_error=initial_error,final_error=final_error))
        # Leave one physical edge out of the already extracted planes; evaluation-only group split.
        if len(plane_edges):
            matrix=radius*np.cross(runner['points'],np.column_stack([runner['gradients'],np.zeros(len(plane_edges))]))
            for edge in sorted(per_edge):
                keep=plane_edges!=edge;om,_,rank,sing=np.linalg.lstsq(matrix[keep],np.ones(keep.sum()),rcond=None)
                residual=matrix[~keep]@om-1
                loo.append(dict(scene=scene,heldout_physical_edge=edge,fit_planes=int(keep.sum()),heldout_planes=int((~keep).sum()),
                    rank=int(rank),condition=float(sing[0]/sing[-1]) if len(sing)==3 and sing[-1]>1e-15 else None,
                    estimated_omega=om.tolist(),axis_error_deg=angle(om,omega),heldout_aperture_RMS=float(np.sqrt(np.mean(residual**2))),
                    truth_grouping_used_only_by_evaluator=True))
        raw_init=result['raw_rhs_S_initializer'];comp_init=next(r for r in result['routes'] if r['method']=='S')['initializer']
        comparison=dict(raw_rhs_error=error(raw_init['raw_omega_rad_s'],omega),compensated_refit_error=error(comp_init['raw_omega_rad_s'],omega),
            raw_rhs_omega=raw_init['raw_omega_rad_s'],compensated_refit_omega=comp_init['raw_omega_rad_s'],
            raw_planes=raw_init['planes'],compensated_planes=comp_init['planes'],
            difference_norm=float(np.linalg.norm(np.array(raw_init['raw_omega_rad_s'])-comp_init['raw_omega_rad_s'])),
            caveat='Same raw events and fixed split; full coordinate correction changes neighbors and fits. Initializer values are diagnostic, not accepted zero measurements.')
        evaluation.append(dict(scene=scene,truth_omega=omega.tolist(),truth_velocity=velocity.tolist(),physical_edge_plane_counts=per_edge,
            minimum_plane_edge_purity=min(plane_purity,default=None),local_holdout_time_RMS_us=summary(local_rms),
            local_holdout_RMS_over_training_gate=summary(local_ratios),zero_motion_heldout_baseline=heldout_motion(points,times,labels,train,hold,[0,0,0],radius) if len(xy) else [],
            routes=routes,translation_path_comparison=comparison))
        if len(xy):ax.scatter(xy[:,0],xy[:,1],c=times*1000,s=1,cmap='viridis',vmin=0,vmax=10)
        ax.set(aspect='equal',xlim=(525,755),ylim=(475,245),title=f'{scene} | {len(xy)} integer crossings, {len(runner["points"])} planes',xlabel='x (px)',ylabel='y (px)')
    assert polarity_mismatches==0 and len(gate_checks)==3
    report=dict(status='passed',S_three_full_chain_gates=gate_checks,independent_max_root_residual=max_root,
        independent_integer_time_rounding_max_s=max_rounded,polarity_mismatches=polarity_mismatches,scenes=evaluation,
        caveat='Geometric edge-crossing implementation check; no DVS noise/refractory/contrast model, not a real-data accuracy claim.')
    write_json(OUT/'synthetic_evaluation.json',report);table(OUT/'synthetic_metrics.csv',rowtable);table(OUT/'heldout_local_motion.csv',local_rows);table(OUT/'leave_one_edge_out.csv',loo)
    fig.suptitle('Fixed ideal material-edge scenes | evaluator only; holdout IDs never fit or optimized')
    fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(OUT/'synthetic_events.png',dpi=150);plt.close(fig)
    print('S full chain passed on all3 rotations; roots/polarity/holdout/evaluator separation checks passed')

if __name__=='__main__':main()
