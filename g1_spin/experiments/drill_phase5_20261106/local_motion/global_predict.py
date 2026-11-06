"""LOPO evaluator compares already frozen omega and2D translation on identical heldpatch templates."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

OUT=Path(__file__).resolve().parent;P4=OUT.parents[1]/'drill_phase4_20261106'
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False))
def summary(a):
    return dict(n=len(a),mean=float(np.mean(a)) if len(a) else None,median=float(np.median(a)) if len(a) else None,p95=float(np.percentile(a,95)) if len(a) else None)
def by_unit(xy,pol,dist,times):
    groups={}
    for p,s,d,t in zip(xy,pol,dist,times):groups.setdefault((int(p[0]),int(p[1]),int(s),int(np.floor(t/.0005+1e-10))),[]).append(float(d))
    return [float(np.mean(v)) for v in groups.values()]

def main():
    real=OUT/'real';path=real/'development_global_LOPO_predictions.json'
    if path.exists():raise FileExistsError('Preserve fixed heldpatch comparison')
    systems=json.loads((real/'development_T_global_systems.json').read_text());manifest=json.loads((real/'development_inputs.json').read_text());old=json.loads((P4/'real_frontend/input_manifest.json').read_text());out=[]
    for system in [s for s in systems if s['global_training_sufficient']]:
        case=next(c for c in manifest if c['case']==system['case']);source=next(c for c in old if c['case']==case['case']);context=np.load(P4/'real_frontend'/source['input_file']);center=np.array(case['center']);radius=case['radius_px'];patch_records=[]
        for leave in system['leave_patch_out']:
            patch=next(p for p in case['patches'] if p['patch']==leave['heldout_patch']);data=np.load(OUT/patch['T_file']);chosen=(data['time_s']>=.005)&(data['time_s']<.006)
            xy=data['xy'][chosen];times=data['time_s'][chosen];pol=data['polarity'][chosen];ids=data['ids'][chosen]
            local=(xy-center)/radius;invalid=np.sum(local*local,axis=1)>=1
            box=patch['box'];roi=[box[0]-15,box[1]-15,box[2]+15,box[3]+15]
            models={'zero':np.zeros(3),'rotation':np.array(leave['rotation']['omega_or_velocity'])};velocity=np.array(leave['translation']['omega_or_velocity']) if leave['translation']['omega_or_velocity'] is not None else None
            segments=[];aggregates={name:{'forward':[],'reverse':[]} for name in ['zero','rotation','translation']}
            for segment in range(8):
                low=.006+segment*.0005;high=low+.0005;t=(low+high)/2
                keep=(context['times']>=low-1e-12)&(context['times']<high-1e-12)&(context['xy'][:,0]>=roi[0])&(context['xy'][:,0]<roi[2])&(context['xy'][:,1]>=roi[1])&(context['xy'][:,1]<roi[3])
                target=context['xy'][keep];tp=context['polarity'][keep];tid=context['ids'][keep]
                entry=dict(segment=segment,target_time_s=t,source_ids=ids.tolist(),target_ids=tid.tolist(),fixed_roi=roi,status='insufficient_template_or_target',scores={})
                if invalid.any():entry['status']='geometry_invalid_source'
                elif len(xy)>=4 and len(target)>=4:
                    point=np.column_stack([local,-np.sqrt(1-np.sum(local*local,axis=1))]);entry['status']='evaluable'
                    for name in ['zero','rotation','translation']:
                        if name=='translation':
                            if velocity is None:continue
                            prediction=xy+(t-times)[:,None]*velocity;visible=np.ones(len(xy),bool)
                        else:
                            warped=Rotation.from_rotvec((t-times)[:,None]*models[name]).apply(point);prediction=warped[:,:2]*radius+center;visible=warped[:,2]<0
                        forward=np.full(len(xy),3.);reverse=np.full(len(target),3.)
                        for polarity in [0,1]:
                            si=np.flatnonzero((pol==polarity)&visible);ti=np.flatnonzero(tp==polarity)
                            if len(si) and len(ti):forward[si]=np.minimum(cKDTree(target[ti]).query(prediction[si])[0],3.);reverse[ti]=np.minimum(cKDTree(prediction[si]).query(target[ti])[0],3.)
                        f=by_unit(xy,pol,forward,times);r=by_unit(target,tp,reverse,context['times'][keep]);aggregates[name]['forward']+=f;aggregates[name]['reverse']+=r
                        entry['scores'][name]=dict(forward=summary(f),reverse=summary(r),equal_direction_mean=float((np.mean(f)+np.mean(r))/2),backface_sources=int((~visible).sum()))
                segments.append(entry)
            scores={}
            for name,data2 in aggregates.items():
                f=data2['forward'];r=data2['reverse'];scores[name]=dict(forward=summary(f),reverse=summary(r),equal_direction_mean=float((np.mean(f)+np.mean(r))/2) if f and r else None)
            rot=scores['rotation']['equal_direction_mean'];trans=scores['translation']['equal_direction_mean']
            patch_records.append(dict(heldout_patch=patch['patch'],training_observation_ids=leave['training_observation_ids'],source_T_ids=ids.tolist(),
                no_heldout_patch_used_in_system=True,segments=segments,aggregate=scores,
                rotation_minus_translation_mean=None if rot is None or trans is None else rot-trans,
                comparison_scope='Same heldpatch T-last1ms source and all fixed ROI/context V targets; no correspondence or target selection by model. Unit-weighted pixel/polarity per0.5ms segment, distance cap3px.',
                accepted_rotation_explanation=False,reason='Local validation does not support>=3patch nonzero motion; these system predictions remain diagnostic, metric not recalibrated here.'))
        out.append(dict(case=case['case'],representation=system['representation'],model_kind=system['model_kind'],patches=patch_records,
            source='Frozen T-only all-qualified global systems; no H/future filtering or new fits during evaluator'))
    write(path,out);print('Saved',len(out),'frozen-system comparisons;',sum(len(r['patches']) for r in out),'heldpatch templates, no new fit')

if __name__=='__main__':main()
