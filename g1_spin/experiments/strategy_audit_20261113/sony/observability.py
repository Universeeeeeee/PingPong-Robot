"""Geometry-only upper-bound observability with hypothetical full 2D velocities.

No real flow measurements, confidence intervals or true RPM errors are inferred.
"""
import json
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent
REPO=OUT.parents[3]

def geometry(points):
    x,y,z=points.T
    n=len(points);A=np.zeros((2*n,5))
    A[0::2,0]=1;A[1::2,1]=1
    A[0::2,3]=z;A[0::2,4]=-y
    A[1::2,2]=-z;A[1::2,4]=x
    singular=np.linalg.svd(A,compute_uv=False)
    col_norm=np.linalg.norm(A,axis=0)
    normalized=np.linalg.svd(A/col_norm,compute_uv=False)
    # Portion of transverse-spin column energy removable by free translation.
    explained=float(np.mean(z)**2/np.mean(z*z))
    return dict(support_points=n,z_mean=float(z.mean()),z_std=float(z.std()),
        translation_explains_transverse_spin_column_energy=explained,
        unabsorbed_fraction=1-explained,svd_singular=singular.tolist(),
        condition=float(singular[0]/singular[-1]),column_normalized_condition=float(normalized[0]/normalized[-1]))

files=[('disk_early','g1_spin/experiments/disk_recheck_20261113/inputs/early_turn_01750000_50ms.npz','linear50',48),
 ('disk_middle','g1_spin/experiments/disk_recheck_20261113/inputs/middle_06350000_50ms.npz','linear50',48),
 ('disk_late','g1_spin/experiments/disk_recheck_20261113/inputs/late_11350000_50ms.npz','linear50',48),
 ('disk_high_translation','g1_spin/experiments/disk_recheck_20261113/inputs/high_motion_00450000_50ms.npz','linear50',48),
 ('drill_fixed','g1_spin/experiments/drill_candidate_validation_20261110/inputs/run1_00835000_candidate.npz','centers',60)]
results=[]
for label,path,key,radius in files:
    path=REPO/path
    with np.load(path) as d:
        e=d['events'];c=d[key];uv=np.column_stack([e['x'],e['y']]).astype(float)-c
    # One support point per 4px body-relative image cell, reducing repeat-event weighting.
    cells=np.floor(uv/4).astype(int);ids=np.unique(cells,axis=0,return_inverse=True)[1]
    xy=np.array([np.median(uv[ids==i],axis=0) for i in range(ids.max()+1)])/radius
    z=-np.sqrt(1-np.sum(xy**2,axis=1));p=np.column_stack([xy,z])
    results.append(dict(label=label,source=str(path.relative_to(REPO)),**geometry(p)))
# Analytic centered cap control makes the degeneration explicit.
for rho in [.1,.3,.5,.8]:
    t=np.linspace(0,2*np.pi,360,endpoint=False)
    radii=np.linspace(.001,rho,30)
    p=np.vstack([np.column_stack([rr*np.cos(t),rr*np.sin(t),np.full(len(t),-np.sqrt(1-rr*rr))]) for rr in radii])
    results.append(dict(label=f'synthetic_cap_rho_{rho}',**geometry(p)))
(OUT/'observability.json').write_text(json.dumps(results,indent=2))
print(json.dumps([{k:r[k] for k in ['label','support_points','translation_explains_transverse_spin_column_energy','condition','column_normalized_condition']} for r in results],indent=2))
