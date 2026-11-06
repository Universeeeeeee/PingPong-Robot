"""T-only fixed bootstrap qualification and global systems, before heldout evaluation."""
import json
import hashlib
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False))

def classification(row):
    m=row['model'];b=row['bootstrap'];valid=bool(m['available'] and m['unique_pixels']>=12 and b['valid_draws']>=24 and b['normal_spread95_deg'] is not None and b['normal_spread95_deg']<=15 and b['speed_deviation95_px_s'] is not None)
    result=dict(training_qualified=valid,training_motion_class='unidentifiable',speed_threshold_px_s=None,bootstrap_valid_minimum=24,
        reason='insufficient_T_support_or_time_normal_bootstrap' if not valid else None)
    if not valid:return result
    v=m['speed_px_s'];delta=b['speed_deviation95_px_s'];threshold=max(50.,2*delta);CI=b['speed_CI95'];result['speed_threshold_px_s']=threshold
    if abs(v)>threshold:result['training_motion_class']='candidate_nonzero'
    elif abs(v)<=50 and CI[0]<=0<=CI[1]:result['training_motion_class']='static_compatible'
    else:result['reason']='speed_not_nonzero_or_static_identifiable'
    return result

def least_squares(matrix,rhs,rank_required):
    if not len(rhs):return dict(rank=0,condition=None,omega_or_velocity=None,RMS=None,available=False)
    solution,_,rank,singular=np.linalg.lstsq(matrix,rhs,rcond=None);cond=float(singular[0]/singular[-1]) if len(singular)>=rank_required and singular[-1]>1e-12 else None
    return dict(rank=int(rank),condition=cond,omega_or_velocity=solution.tolist(),RMS=float(np.sqrt(np.mean((matrix@solution-rhs)**2))),
        available=bool(rank==rank_required and cond is not None and cond<=100))

def build_system(rows,center,radius):
    qualified=[r for r in rows if r['T_qualification']['training_qualified']];patches=sorted(set(r['patch'] for r in qualified));matrix=[];normals=[];rhs=[];missing=[]
    for r in qualified:
        m=r['model'];local=(np.array(m['center'])-center)/radius
        if local@local>=1:missing.append(r['observation_id']);continue
        s=np.r_[local,-np.sqrt(1-local@local)];n=np.array(m['normal']);matrix.append(radius*np.cross(s,[n[0],n[1],0.]));normals.append(n);rhs.append(m['speed_px_s'])
    matrix=np.array(matrix).reshape(-1,3);normals=np.array(normals).reshape(-1,2);rhs=np.array(rhs)
    if missing:
        invalid=dict(rank=0,condition=None,omega_or_velocity=None,RMS=None,available=False,failure_reason='geometry_invalid')
        return dict(T_qualified_count=len(qualified),independent_spatial_patches=len(patches),observation_ids=[r['observation_id'] for r in qualified],
            spherical_geometry_invalid=missing,rotation=invalid.copy(),translation=invalid.copy(),leave_patch_out=[],
            global_training_sufficient=False,rotation_chain_allowed=False,translation_LOPO_rank2=False,
            failure_reason='geometry_invalid',fit_boundary='Whole system rejected; no silent observation deletion or mismatched LOPO row masks.')
    rotation=least_squares(matrix,rhs,3);translation=least_squares(normals,rhs,2);lopo=[]
    for patch in patches:
        keep=np.array([r['patch']!=patch for r in qualified]);rot=least_squares(matrix[keep],rhs[keep],3);trans=least_squares(normals[keep],rhs[keep],2)
        lopo.append(dict(heldout_patch=patch,training_observation_ids=[r['observation_id'] for r,k in zip(qualified,keep) if k],rotation=rot,translation=trans))
    sufficient=bool(len(qualified)>=12 and len(patches)>=4 and not missing and rotation['available'] and all(r['rotation']['rank']==3 for r in lopo))
    # Same raw covector constraints retain actual estimated v for static-compatible rows; never forced0.
    return dict(T_qualified_count=len(qualified),independent_spatial_patches=len(patches),observation_ids=[r['observation_id'] for r in qualified],
        spherical_geometry_invalid=missing,rotation=rotation,translation=translation,leave_patch_out=lopo,
        global_training_sufficient=sufficient,rotation_chain_allowed=sufficient,
        translation_LOPO_rank2=bool(lopo and all(r['translation']['rank']==2 for r in lopo)),
        equations='r*(s cross [nx,ny,0]) dot omega = vn; n dot u = vn; fixed B1 cdot0',
        fit_boundary='All and only T-qualified rows; no H/future/permutation filtering; diagnostic LS may be reported below gate, no rotation claim or CMax.')

def main():
    real=OUT/'real';path=real/'development_training_qualification.json'
    if path.exists():raise FileExistsError('Preserve T-only qualification/global freeze')
    rows=json.loads((real/'development_training.json').read_text());manifest=json.loads((real/'development_inputs.json').read_text())
    for index,r in enumerate(rows):r['observation_id']=index;r['T_qualification']=classification(r)
    systems=[]
    for case in manifest:
        for rep in ['R_all','R_latest']:
            for kind in ['M_time','M_spatial']:
                local=[r for r in rows if r['case']==case['case'] and r['representation']==rep and r['model_kind']==kind]
                sys=build_system(local,np.array(case['center']),case['radius_px']);systems.append(dict(case=case['case'],representation=rep,model_kind=kind,**sys))
    write(path,rows);write(real/'development_T_global_systems.json',systems)
    write(real/'development_qualification_freeze.json',dict(stage='T-only fit/qualification/global systems before H/future evaluation',
        files={f:digest(OUT/f) for f in ['real/development_training.json','real/development_training_freeze.json','qualify_training.py',
        'real/development_training_qualification.json','real/development_T_global_systems.json']},local_rows=len(rows),systems=len(systems),
        all_constraints_selected_from_T_only=True,H_future_never_read=True))
    print('FROZEN T qualifications',len(rows),'systems',len(systems),'global sufficient',sum(s['global_training_sufficient'] for s in systems))

if __name__=='__main__':main()
