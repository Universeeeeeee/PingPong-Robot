"""Controlled geometry-only point trajectories. Not a physical event simulation."""
import argparse
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.spatial.transform import Rotation
from common import save,verify_lock
from frozen import cmax

AXES=[[0,1,0],[0,0,1],[.6,.7,.4]]
CASES=[('ideal',0.),('center_px',.5),('center_px',1.),('center_px',2.),
       ('velocity_px_s',50.),('velocity_px_s',100.),('lag_ms',-5.),('lag_ms',-1.),('lag_ms',1.),('lag_ms',5.),
       ('radius_px',-1.),('radius_px',1.),('radius_px',3.),('drop_fraction',.5)]

def center(t):return np.column_stack([128+300*t+2*np.sin(2*np.pi*8*t),128+100*t+np.cos(2*np.pi*8*t)])

def generate(axis,seed):
    rng=np.random.default_rng(seed);landmarks=rng.normal(size=(24,3));landmarks/=np.linalg.norm(landmarks,axis=1,keepdims=True)
    ts=np.arange(168)/8400;t=ts+.1;omega=np.array(axis,float);omega=omega/np.linalg.norm(omega)*(247*2*np.pi/60)
    p=np.einsum('tij,kj->tki',Rotation.from_rotvec(ts[:,None]*omega).as_matrix(),landmarks)
    visible=p[:,:,2]<-.12;t=np.broadcast_to(t[:,None],visible.shape)[visible];uv=p[:,:,:2][visible]*60+center(t)
    return uv,t,omega

def hypothesis(t,kind,value):
    c=center(t);radius=60.
    if kind=='center_px':c[:,0]+=value
    elif kind=='velocity_px_s':c[:,0]+=value*(t-(t[0]+t[-1])/2)
    elif kind=='lag_ms':c=center(t+value/1000)
    elif kind=='radius_px':radius+=value
    return c,radius

def main(out):
    verify_lock();out=Path(out);out.mkdir(parents=True,exist_ok=False);results=[]
    save(out/'protocol.json',dict(axes=AXES,seeds=[0,1],rpm=247,window_ms=20,cases=CASES,common_inner_mask=.8,
         caveat='Ideal visible surface point tracks. No DVS response, real texture, quantization, blur or illumination. Thresholds are exploratory only.'))
    for axid,axis in enumerate(AXES):
        for seed in [0,1]:
            uv,t,truth=generate(axis,seed);hyps=[hypothesis(t,k,v) for k,v in CASES]
            keep=np.logical_and.reduce([np.linalg.norm(uv-c,axis=1)<=.8*r for c,r in hyps]);indices=np.flatnonzero(keep)
            save(out/f'axis{axid}_seed{seed}_input.json',dict(total=len(t),common_count=len(indices),common_indices=indices.tolist(),truth_omega=truth.tolist()))
            for (kind,value),(c,r) in zip(CASES,hyps):
                ix=indices[::2] if kind=='drop_fraction' else indices
                rel=uv[ix]-c[ix];start=perf_counter();est=cmax.estimate(cmax.lift(rel[:,0],rel[:,1],0,0,r),t[ix]-.1,max_rpm=3500)
                angle=None
                if est['valid']:angle=float(np.degrees(np.arccos(np.clip(np.dot(est['axis_camera'],truth/np.linalg.norm(truth)),-1,1))))
                results.append(dict(axis=axid,seed=seed,kind=kind,value=value,events=len(ix),valid=est['valid'],rpm=est['spin_rpm'],
                  rpm_abs_error=abs(est['spin_rpm']-247) if est['valid'] else None,axis_error_deg=angle,reason=est['failure_reason'],elapsed_s=perf_counter()-start))
            print('axis',axid,'seed',seed,flush=True)
    save(out/'results.json',results)
    grouped=[]
    for kind,value in CASES:
        rr=[r for r in results if r['kind']==kind and r['value']==value];good=[r for r in rr if r['valid']]
        grouped.append(dict(kind=kind,value=value,valid=len(good),planned=len(rr),rpm_abs_error_median=float(np.median([r['rpm_abs_error'] for r in good])) if good else None,
            axis_error_median_deg=float(np.median([r['axis_error_deg'] for r in good])) if good else None))
    save(out/'summary.json',grouped)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();main(a.out)
