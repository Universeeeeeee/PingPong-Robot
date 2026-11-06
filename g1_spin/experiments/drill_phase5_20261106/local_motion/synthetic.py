"""Synthetic-only observations, oracle metric calibration and same-support2x2 fits."""
import json
from pathlib import Path
import hashlib
import numpy as np
from scipy.spatial.transform import Rotation
from models import representation,fit,residual,equal_unit_distances,canonical

OUT=Path(__file__).resolve().parent;P5=OUT.parent;P4=P5.parent/'drill_phase4_20261106'
SCENES=['base','reverse','second','translation','mixed','zero','tail','flicker_shuffle']

def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def select_patches(xy,times):
    T=(times<.006)&((xy[:,0]+2*xy[:,1])%5!=0);cells=xy[T]//12;unique,count=np.unique(cells,axis=0,return_counts=True)
    choices=[]
    for i in sorted(range(len(unique)),key=lambda i:(-count[i],int(unique[i,0]),int(unique[i,1]))):
        if count[i]<20:continue
        center=unique[i]*12+6
        if any(np.linalg.norm(center-np.array(c['center']))<16 for c in choices):continue
        choices.append(dict(patch=len(choices),cell=unique[i].tolist(),center=center.tolist(),box=np.r_[unique[i]*12,unique[i]*12+12].tolist(),selection_T_events=int(count[i])))
        if len(choices)==6:break
    return choices
def box_mask(xy,box):return (xy[:,0]>=box[0])&(xy[:,0]<box[2])&(xy[:,1]>=box[1])&(xy[:,1]<box[3])

def inputs():
    (OUT/'inputs').mkdir(exist_ok=True);(OUT/'evaluator_truth').mkdir(exist_ok=True);meta=[]
    for name in SCENES:
        if name=='flicker_shuffle':
            y=np.arange(348,373);x=np.array([639,640,641]);xy=np.tile(np.array([(a,b) for a in x for b in y]),(16,1))
            time=np.repeat(np.arange(16)*.0005+.00025,len(x)*len(y));pol=np.repeat(np.arange(16)%2,len(x)*len(y)).astype(np.uint8)
            rng=np.random.default_rng(7)
            for p in [0,1]:
                ix=np.flatnonzero(pol==p);time[ix]=time[ix][rng.permutation(len(ix))]
            truth=dict(edge_ids=np.zeros(len(time),np.int8),source_time_s=time.copy(),tail_offset_s=np.zeros(len(time)),omega=np.zeros(3),velocity=np.zeros(2),normals=np.array([[1.,0,0]]),anchors=np.array([[0.,0,-1.]]))
            geo=dict(center_px=[640.,360.],radius_px=120.,reference_us=15000)
        else:
            source='base' if name=='tail' else name;d=np.load(P4/'inputs'/f'{source}.npz');t=np.load(P4/'truth'/f'{source}.npz');xy=np.column_stack([d['x'],d['y']]).astype(int);time=(d['t_us']-10000)*1e-6;pol=d['polarity']
            geo=json.loads((P4/'inputs'/f'{source}_geometry.json').read_text());edge=t['edge_ids'];source_time=t['exact_time_s'];offset=np.zeros(len(time))
            if name=='tail':
                arrays=[]
                for delta in [0.,.0005,.001]:
                    keep=(time+delta>=0)&(time+delta<.01)
                    arrays.append((xy[keep],time[keep]+delta,pol[keep],edge[keep],source_time[keep],np.full(keep.sum(),delta)))
                xy=np.concatenate([a[0] for a in arrays]);time=np.concatenate([a[1] for a in arrays]);pol=np.concatenate([a[2] for a in arrays]);edge=np.concatenate([a[3] for a in arrays]);source_time=np.concatenate([a[4] for a in arrays]);offset=np.concatenate([a[5] for a in arrays])
            truth=dict(edge_ids=edge,source_time_s=source_time,tail_offset_s=offset,omega=t['omega'],velocity=t['velocity'],normals=t['normals'],anchors=t['anchors'])
        order=np.argsort(time,kind='stable');xy=xy[order];time=time[order];pol=pol[order]
        for key in ['edge_ids','source_time_s','tail_offset_s']:truth[key]=truth[key][order]
        np.savez_compressed(OUT/'inputs'/f'{name}.npz',xy=xy,time_s=time,polarity=pol,ids=np.arange(len(time)))
        write(OUT/'inputs'/f'{name}_geometry.json',geo);np.savez_compressed(OUT/'evaluator_truth'/f'{name}.npz',**truth)
        meta.append(dict(scene=name,events=len(time),positive_oracle_calibration=name in ['base','reverse','second','translation','mixed'],
            extra_control=name in ['tail','flicker_shuffle'],patches=select_patches(xy,time),events_sha256=digest(OUT/'inputs'/f'{name}.npz')))
    write(OUT/'input_manifest.json',meta);return meta

def metric(model,xy,times,pol,mask):
    distance=residual(model,xy[mask],times[mask]);out=equal_unit_distances(xy[mask],times[mask],pol[mask],distance)
    out['raw_event_ids']=np.flatnonzero(mask).tolist();out['raw_distances_px']=distance.tolist();return out

def oracle_model(xy,times,training,truth,geometry,edge,t0):
    center=np.array(geometry['center_px']);radius=geometry['radius_px'];omega=truth['omega'];velocity=truth['velocity'];normal_material=truth['normals'][edge]
    # Evaluate true instantaneous edge normal at the material point closest to patch T centroid.
    local=(xy[training].mean(axis=0)-center-(t0-.005)*velocity)/radius
    if local@local>=1:return None
    point=np.r_[local,-np.sqrt(1-local@local)];normal3=Rotation.from_rotvec(omega*(t0-.005)).apply(normal_material)
    J=np.vstack([np.eye(2),-point[:2]/point[2]])
    gradient=J.T@normal3/radius;norm=np.linalg.norm(gradient)
    if norm<1e-12:return None
    normal=gradient/norm;full_velocity=radius*np.cross(omega,point)[:2]+velocity;speed=float(normal@full_velocity)
    normal,speed,_=canonical(normal,speed,0.);x0=xy[training].mean(axis=0);intercept=float(np.mean((xy[training]-x0)@normal-speed*(times[training]-t0)))
    return dict(available=True,normal=normal.tolist(),speed_px_s=speed,intercept_px=intercept,center=x0.tolist(),reference_s=t0,
        scope='Evaluator true instantaneous normal/speed; intercept estimated only from same-edge T, labels never runner inputs')

def main():
    if (OUT/'oracle_calibration.json').exists():raise FileExistsError('Preserve one fixed calibration pass')
    meta=inputs();oracles=[];models=[];positive=[]
    for record in meta:
        name=record['scene'];d=np.load(OUT/'inputs'/f'{name}.npz');truth=np.load(OUT/'evaluator_truth'/f'{name}.npz');geo=json.loads((OUT/'inputs'/f'{name}_geometry.json').read_text())
        xy=d['xy'];times=d['time_s'];pol=d['polarity'];ids=d['ids'];Hpixel=(xy[:,0]+2*xy[:,1])%5==0
        for patch in record['patches']:
            patch_mask=box_mask(xy,patch['box'])
            for segment in range(3):
                low=segment*.002;high=low+.002;t0=(low+high)/2
                for polarity in [0,1]:
                    mask=patch_mask&(pol==polarity);T=mask&(times>=low-1e-12)&(times<high-1e-12)&~Hpixel
                    H=mask&(times>=low-1e-12)&(times<high-1e-12)&Hpixel
                    future=mask&(times>=high-1e-12)&(times<high+.001-1e-12)&(Hpixel if segment<2 else True)
                    identity=dict(scene=name,patch=patch['patch'],segment=segment,polarity=polarity,start_s=low,end_s=high,training_ids=ids[T].tolist(),H_ids=ids[H].tolist(),future_ids=ids[future].tolist())
                    if T.sum():
                        for edge in np.unique(truth['edge_ids'][T]):
                            edgeT=T&(truth['edge_ids']==edge);edgeH=H&(truth['edge_ids']==edge);edgeF=future&(truth['edge_ids']==edge)
                            oracle=oracle_model(xy,times,edgeT,truth,geo,int(edge),t0)
                            if oracle is None:continue
                            h=metric(oracle,xy,times,pol,edgeH);h['evaluable']=h['unique_pixels']>=4
                            f=metric(oracle,xy,times,pol,edgeF)
                            eligible=len(np.unique(xy[edgeT],axis=0))>=12
                            row=dict(**identity,edge_eval_only=int(edge),T_unique_pixels=len(np.unique(xy[edgeT],axis=0)),model=oracle,H=h,future=f,
                                positive_calibration=record['positive_oracle_calibration'],support_eligible=eligible)
                            oracles.append(row)
                            if record['positive_oracle_calibration'] and eligible:
                                for target,score in [('H',h),('future',f)]:
                                    if score['evaluable']:positive.append(dict(scene=name,patch=patch['patch'],segment=segment,polarity=polarity,edge=int(edge),target=target,p95=score['p95'],units=score['units']))
                    for rep in ['R_all','R_latest']:
                        training=np.flatnonzero(T);chosen=training[representation(xy[training],times[training],ids[training],rep)]
                        for modelkind in ['M_time','M_spatial']:
                            model=fit(xy[chosen],times[chosen],modelkind)
                            h=metric(model,xy,times,pol,H) if model['available'] else None
                            if h is not None:h['evaluable']=h['unique_pixels']>=4
                            models.append(dict(**identity,representation=rep,model_kind=modelkind,support_ids=ids[chosen].tolist(),model=model,
                                H=h,future=metric(model,xy,times,pol,future) if model['available'] else None,
                                status='synthetic_diagnostic_only_no_final_real_motion_qualification'))
        print('synthetic local',name,'events',len(times),'patches',len(record['patches']),flush=True)
    assert positive,'No eligible oracle observations; report unsupported metric rather than invent tau'
    worst=max(positive,key=lambda r:r['p95']);tau=max(.5,1.25*worst['p95']);calibration=dict(status='passed' if tau<3 else 'resolution_failed',tau_px=tau,minimum_tau_px=.5,max_allowed_tau_exclusive_px=3.,
        clean_positive_oracle_scores=len(positive),worst_positive_oracle=worst,all_positive_P95=positive,
        tail_and_flicker_role='Observed negative/misspecification controls; not allowed to inflate clean positive tau',
        real_execution_authorized=False,scope='Synthetic metric gate only; if tau>=3 do not increase thresholds or run real models.')
    write(OUT/'oracle_observations.json',oracles);write(OUT/'model_observations.json',models);write(OUT/'oracle_calibration.json',calibration)
    print('ORACLE',calibration['status'],'tau',tau,'worst',worst)

if __name__=='__main__':main()
