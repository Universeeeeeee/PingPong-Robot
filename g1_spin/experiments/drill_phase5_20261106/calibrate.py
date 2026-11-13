"""Measurement positives and exact repetition lower bounds, without a spin solver."""
import numpy as np
from diagnose import OUT,decomposition,write

def main():
    x,y=np.meshgrid(np.arange(13),np.arange(13));xy=np.column_stack([x.ravel(),y.ravel()]);gradient=np.array([100e-6,150e-6]);times=xy@gradient+.001
    exact=decomposition(xy,times,np.zeros(len(xy),int));X=np.column_stack([xy,np.ones(len(xy))]);coef=np.linalg.lstsq(X,times,rcond=None)[0]
    residual=float(np.sqrt(np.mean((X@coef-times)**2)));assert np.allclose(coef[:2],gradient,atol=1e-15) and residual<1e-15
    assert exact['irreducible_RMS_us']==0 and not exact['impossible_under_old_gate'];assert exact['maximum_single_value_R2']==1
    delta=np.array([-200e-6,0,200e-6]);repeated_xy=np.repeat(xy,3,axis=0);repeated_t=np.repeat(times,3)+np.tile(delta,len(xy));pol=np.zeros(len(repeated_t),int)
    repeated=decomposition(repeated_xy,repeated_t,pol);expected=float(np.sqrt(np.mean(delta*delta)))
    assert np.isclose(repeated['irreducible_RMS_us']*1e-6,expected,atol=1e-15)
    # Any known predictor must have RMS at least this lower bound; per-pixel means attain it exactly.
    predicted=np.repeat(times,3);measured=float(np.sqrt(np.mean((repeated_t-predicted)**2)));assert np.isclose(measured,expected)
    assert repeated['impossible_under_old_gate']
    # Opposite polarities each get their own timestamp: combined grouping must not merge them.
    px=np.tile(xy,(2,1));pt=np.r_[times,times+.001];pp=np.r_[np.zeros(len(xy)),np.ones(len(xy))]
    perpol=decomposition(px,pt,pp);merged=decomposition(px,pt)
    assert perpol['irreducible_RMS_us']==0 and np.isclose(merged['irreducible_RMS_us'],500.)
    # A broad case is equivalent to R²>=.9775; tiny spread must respect20us floor instead.
    assert np.isclose(1-.15**2,.9775)
    tiny=decomposition(np.zeros((3,2)),np.array([0,1e-6,2e-6]),np.zeros(3))
    assert tiny['gate_branch']=='20us_floor' and not tiny['impossible_under_old_gate']
    report=dict(status='passed',completed_date='2026-11-07',tests=['exact_plane_gradient_RMS_and_R2_positive','repeated_pixel_irreducible_variance_exact',
        'mean_predictor_attains_lower_bound','polarity_groups_not_silently_merged','R2_gate_and20us_floor_branches'],
        exact_plane=dict(events=len(xy),gradient=gradient.tolist(),fitted_gradient=coef[:2].tolist(),plane_RMS_us=residual*1e6,lower_bound_us=exact['irreducible_RMS_us']),
        repeated_pixel=dict(events=len(repeated_t),expected_lower_us=expected*1e6,measured_lower_us=repeated['irreducible_RMS_us'],gate_us=repeated['old_gate_us'],impossible=repeated['impossible_under_old_gate']),
        separate_polarity_lower_us=perpol['irreducible_RMS_us'],merged_polarity_lower_us=merged['irreducible_RMS_us'],floor_case=tiny,
        caveat='Calibrates measurement algebra, not DVS or physical spin accuracy.')
    write(OUT/'measurement_calibration.json',report);print(report)

if __name__=='__main__':main()
