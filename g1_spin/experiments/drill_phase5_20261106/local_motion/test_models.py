"""Analytic unit/support checks for2x2 local fits and equal observation scoring."""
import json
import numpy as np
from pathlib import Path
from models import representation,fit,residual,equal_unit_distances

OUT=Path(__file__).resolve().parent

def main():
    n=np.array([.6,.8]);tangent=np.array([-.8,.6]);speed=240.
    times=np.repeat(np.linspace(0,.002,17),13);along=np.tile(np.linspace(-5,5,13),17)
    xy=np.array([640.,360.])+speed*(times-.001)[:,None]*n+along[:,None]*tangent
    records=[]
    for name in ['M_time','M_spatial']:
        model=fit(xy,times,name);assert model['available']
        assert np.allclose(model['normal'],n,atol=1e-9) and np.isclose(model['speed_px_s'],speed,atol=1e-8)
        assert np.max(residual(model,xy,times))<1e-9
        shifted=fit(xy+[13.,-27.],times+.03,name);assert np.allclose(shifted['normal'],n,atol=1e-9) and np.isclose(shifted['speed_px_s'],speed,atol=1e-8)
        scaled=fit(xy,times*1000,name);assert np.isclose(scaled['speed_px_s'],speed/1000,atol=1e-9) and np.allclose(scaled['normal'],n,atol=1e-9)
        records.append(dict(model=name,normal=model['normal'],speed_px_s=model['speed_px_s'],RMS_px=model['normal_residual_RMS_px']))
    # Fixed line is identifiable for pixel-space regression with speed0, despite undefined inverse time-plane.
    stationary=np.column_stack([np.full(len(times),640.),np.tile(np.linspace(354,366,13),17)])
    model=fit(stationary,times,'M_spatial');assert model['available'] and np.allclose(model['normal'],[1,0]) and model['speed_px_s']==0
    assert not fit(stationary,times,'M_time')['available']
    assert not fit(xy,np.zeros(len(times)),'M_spatial')['available']
    symmetric=np.tile(np.array([[-1,-1],[-1,1],[1,-1],[1,1]],float),(20,1));independent=np.repeat(np.linspace(0,.002,20),4)
    assert fit(symmetric,independent,'M_spatial')['failure_reason']=='normal_direction_not_unique'
    # Covariance lower eigenvalue equals the direct minimized pixel residual variance.
    corrupted=xy.copy();corrupted[:,0]+=.02*np.sin(np.arange(len(xy)))
    check=fit(corrupted,times,'M_spatial');assert np.isclose(check['normal_residual_RMS_px']**2,check['eigenvalues_px2'][0],atol=1e-12)
    # Latest uses every input pixel once, greatest timestamp and then greatest physicalID; no fit-dependent selection.
    pixels=np.array([[1,1],[1,1],[1,1],[2,2]]);ticks=np.array([.001,.0015,.0015,.001]);ids=np.array([4,5,9,3])
    chosen=representation(pixels,ticks,ids,'R_latest');assert np.array_equal(chosen,[2,3])
    assert np.array_equal(representation(pixels,ticks,ids,'R_all'),np.arange(4))
    # Repeating one(pixel,pol,bin)100 times cannot overwhelm three other equally scored units.
    px=np.r_[np.tile([0,0],(100,1)),[[1,0],[2,0],[3,0]]];tt=np.r_[np.full(100,.0001),[.0001,.0006,.0006]];pol=np.zeros(103,int);dist=np.r_[np.zeros(100),[2.,2.,2.]]
    score=equal_unit_distances(px,tt,pol,dist);assert score['units']==4 and score['time_bins']==2 and np.isclose(score['mean'],1.5) and score['median']==2
    report=dict(status='passed',tests=['time_and_pixel_fits_recover_exact_motion','translation_and_units_preserved','stationary_v0_identifiable_in_pixel_model',
        'time_variance_and_direction_degeneracy','spatial_covariance_eigenvalue_equals_pixel_MSE','same_support_representation_latest_tie','pixel_polarity_bin_equal_weight_no_event_count_bias'],
        analytic_models=records,stationary_model=model,equal_unit_test=score,scope='Synthetic measurement/code checks only, not real parameter selection.')
    (OUT/'model_verification.json').write_text(json.dumps(report,indent=2));print('PASSED',len(report['tests']),'model/measurement test groups')

if __name__=='__main__':main()
