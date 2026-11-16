"""Diagnose coarse normal constraints without rerunning or selecting CMAX fits."""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

OUT = Path(__file__).resolve().parent
FRONT = OUT.parent/'drill_phase4_20261106/real_frontend'


def unit(vector):
    magnitude = np.linalg.norm(vector)
    return vector/magnitude if magnitude>1e-12 else None


def angle(first, second):
    if first is None or second is None:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(first),unit(second)),-1,1))))


def solve(matrix, rhs):
    omega, _, rank, singular = np.linalg.lstsq(matrix,rhs,rcond=None)
    return dict(omega=omega.tolist(), rpm=float(np.linalg.norm(omega)*60/(2*np.pi)),
                axis=unit(omega).tolist() if unit(omega) is not None else None,
                rank=int(rank), condition=float(singular[0]/singular[-1]) if rank==3 else None,
                residual_rms_px_s=float(np.sqrt(np.mean((matrix@omega-rhs)**2))))


def main():
    pilots=json.loads((OUT/'coarse_seed/results.json').read_text())
    manifest={r['case']:r for r in json.loads((FRONT/'input_manifest.json').read_text())}
    rows=[]
    for pilot in pilots:
        if pilot['variant']!='clean' or not pilot['initializer']['available']:
            continue
        obs=[r for r in pilot['initializer']['observations'] if r['eligible']]
        matrix=np.array([r['constraint'] for r in obs])
        rhs=np.array([r['rhs_px_s'] for r in obs])
        patch=np.array([r['patch'] for r in obs])
        polarity=np.array([r['polarity'] for r in obs])
        base=solve(matrix,rhs)
        record=dict(case=pilot['case'],baseline=base,variants=[],leave_one_patch_out=[])
        residual=matrix@np.array(base['omega'])-rhs
        # Robust scale from the current residual only, independent of spin labels.
        scale=max(1.,1.4826*np.median(abs(residual-np.median(residual))))
        robust=least_squares(lambda w:matrix@w-rhs,base['omega'],loss='soft_l1',f_scale=scale)
        record['variants'].append(dict(name='soft_l1',scale_px_s=float(scale),omega=robust.x.tolist(),
                                       rpm=float(np.linalg.norm(robust.x)*60/(2*np.pi)),
                                       axis_change_deg=angle(robust.x,base['omega'])))
        for pol in [0,1]:
            selected=polarity==pol
            if selected.sum()>=3:
                result=solve(matrix[selected],rhs[selected])
                result.update(name=f'polarity_{pol}',constraints=int(selected.sum()),axis_change_deg=angle(result['omega'],base['omega']))
                record['variants'].append(result)
        for excluded in np.unique(patch):
            result=solve(matrix[patch!=excluded],rhs[patch!=excluded])
            result.update(excluded_patch=int(excluded),axis_change_deg=angle(result['omega'],base['omega']))
            record['leave_one_patch_out'].append(result)
        entry=manifest[pilot['case']]
        lane=next((r for r in entry['lanes'] if r['name']=='secant' and r['available']),None)
        if lane:
            velocity=np.array(lane['velocity'])
            center=np.array(entry['center']);radius=entry['radius_px']
            corrected_matrix=[];corrected_rhs=[]
            for obs_row in obs:
                fit=obs_row['fit'];normal=np.r_[fit['normal'],0]
                at_center=center+velocity*(fit['reference_s']-lane['reference_time_s'])
                xy=(np.array(fit['center'])-at_center)/radius
                point=np.r_[xy,-np.sqrt(1-xy@xy)]
                corrected_matrix.append(radius*np.cross(point,normal))
                corrected_rhs.append(fit['speed_px_s']-np.dot(velocity,fit['normal']))
            result=solve(np.array(corrected_matrix),np.array(corrected_rhs))
            result.update(name='supplied_secant_center_motion',center_velocity_px_s=velocity.tolist(),
                          axis_change_deg=angle(result['omega'],base['omega']))
            record['variants'].append(result)
        rows.append(record)
    (OUT/'constraint_sensitivity.json').write_text(json.dumps(dict(
        scope='Input-constraint diagnostics only. G and spin references not used. Robustness checks do not select a production initializer.',
        rows=rows),indent=2,allow_nan=False)+'\n')
    for row in rows:
        print(row['case'],'base',round(row['baseline']['rpm'],1),'LOPO max angle',
              round(max(v['axis_change_deg'] for v in row['leave_one_patch_out']),1),
              'variants',[(v['name'],round(v['rpm'],1),round(v['axis_change_deg'],1)) for v in row['variants']])

if __name__=='__main__':
    main()
