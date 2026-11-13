"""Independent metric units, oracle differential and same-support2x2 checks."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from synthetic import OUT,P4,digest,write

def main():
    models=json.loads((OUT/'model_observations.json').read_text());oracles=json.loads((OUT/'oracle_observations.json').read_text());manifest=json.loads((OUT/'input_manifest.json').read_text())
    support={};max_oracle_normal=0.;max_oracle_speed=0.;metric_errors=[];positive=[]
    inputs={}
    for r in manifest:
        path=OUT/'inputs'/f'{r["scene"]}.npz';assert digest(path)==r['events_sha256'];inputs[r['scene']]=np.load(path)
    for r in models:
        key=(r['scene'],r['patch'],r['segment'],r['polarity'],r['representation']);support.setdefault(key,[]).append(r['support_ids'])
        d=inputs[r['scene']];ix=np.array(r['support_ids'],int);xy=d['xy'];times=d['time_s'];pol=d['polarity']
        assert np.all((times[ix]>=r['start_s']-1e-12)&(times[ix]<r['end_s']-1e-12)) and np.all(pol[ix]==r['polarity'])
        assert np.all((xy[ix,0]+2*xy[ix,1])%5!=0)
        if r['representation']=='R_latest':
            assert len(np.unique(xy[ix],axis=0))==len(ix)
            training=np.array(r['training_ids'],int)
            for i in ix:
                same=training[np.all(xy[training]==xy[i],axis=1)];expected=max(same,key=lambda j:(times[j],int(j)));assert i==expected
        if r['model']['available']:
            m=r['model'];x=xy[ix]-m['center'];dt=times[ix]-m['reference_s'];res=x@np.array(m['normal'])-m['speed_px_s']*dt-m['intercept_px']
            assert np.isclose(np.sqrt(np.mean(res*res)),m['normal_residual_RMS_px'],atol=1e-10)
            if r['model_kind']=='M_spatial':
                assert np.isclose(np.mean(res*res),m['eigenvalues_px2'][0],atol=1e-10)
            for target in ['H','future']:
                score=r[target];groups={}
                for key2,value in zip(score['raw_event_ids'],score['raw_distances_px']):
                    assert key2 not in set(r['training_ids']);pixel=xy[key2];binidx=int(np.floor(times[key2]/.0005+1e-10))
                    groups.setdefault((int(pixel[0]),int(pixel[1]),int(pol[key2]),binidx),[]).append(value)
                means=[sum(v)/len(v) for v in groups.values()]
                if means:
                    err=abs(float(np.percentile(means,95))-score['p95']);metric_errors.append(err);assert err<1e-12
    assert all(len(v)==2 and v[0]==v[1] for v in support.values())
    for r in oracles:
        truth=np.load(OUT/'evaluator_truth'/f'{r["scene"]}.npz');geo=json.loads((OUT/'inputs'/f'{r["scene"]}_geometry.json').read_text());omega=truth['omega'];velocity=truth['velocity'];normal=truth['normals'][r['edge_eval_only']];radius=geo['radius_px'];center=np.array(geo['center_px']);m=r['model'];point=np.array(m['center']);time=m['reference_s']
        def function(pixel,t):
            xy=(pixel-center-(t-.005)*velocity)/radius;s=np.r_[xy,-np.sqrt(1-xy@xy)]
            return float(Rotation.from_rotvec(-omega*(t-.005)).apply(s)@normal)
        eps=1e-4;g=np.array([(function(point+np.eye(2)[i]*eps,time)-function(point-np.eye(2)[i]*eps,time))/(2*eps) for i in range(2)])
        ft=(function(point,time+1e-8)-function(point,time-1e-8))/(2e-8);n=g/np.linalg.norm(g);speed=-ft/np.linalg.norm(g)
        if n[0]<-1e-12 or (abs(n[0])<=1e-12 and n[1]<0):n=-n;speed=-speed
        ne=float(np.max(abs(n-m['normal'])));ve=abs(speed-m['speed_px_s']);max_oracle_normal=max(max_oracle_normal,ne);max_oracle_speed=max(max_oracle_speed,ve)
        assert ne<1e-7 and ve<1e-3
        if r['positive_calibration'] and r['support_eligible']:
            for target in ['H','future']:
                if r[target]['evaluable']:positive.append(r[target]['p95'])
    calibration=json.loads((OUT/'oracle_calibration.json').read_text());assert len(positive)==calibration['clean_positive_oracle_scores'] and np.isclose(max(.5,1.25*max(positive)),calibration['tau_px'])
    assert calibration['tau_px']<3 and json.loads((OUT/'model_verification.json').read_text())['status']=='passed'
    prior=json.loads((P4/'manifest_task13.json').read_text())
    for rel,h in prior['files'].items():assert digest(P4/rel)==h['sha256'],rel
    source=json.loads((OUT.parent/'manifest_raw_structure.json').read_text())
    for rel,h in source['files'].items():assert digest(OUT.parent/rel)==h['sha256'],rel
    result=dict(status='passed',completed_date='2026-11-07',synthetic_scenes=len(manifest),model_records=len(models),oracle_records=len(oracles),
        same_representation_model_support_equal=True,all_training_T_only=True,latest_same_segment_causal=True,equal_pixel_polarity_bin_P95_recomputed=True,
        max_metric_difference=max(metric_errors,default=0.),max_oracle_normal_finite_difference_error=max_oracle_normal,
        max_oracle_speed_finite_difference_error_px_s=max_oracle_speed,tau_px=calibration['tau_px'],clean_positive_calibration_scores=len(positive),
        raw_structure_manifest_unchanged=True,stage4_task13_manifest_unchanged=True,no_real_model_execution=True,
        limitation='Intermediate synthetic measurement gate; final T-only training qualification/bootstrap/permutation and real protocol are separate.')
    write(OUT/'verification_synthetic.json',result);print(result)

if __name__=='__main__':main()
