"""Deterministic measurement/cause-order controls; no new scene or real tuning."""
import numpy as np
from helpers import OUT,write_json
from frontends import extract
from prediction import predict,future,compare_future
from pipeline import accepted_prediction

def main():
    tests=[]
    # Latest physical pixel has moved outside lane radius: old in-radius event must not return.
    xy=np.r_[np.array([[3.,0.],[5.,0.],[0.,0.]]),np.column_stack([np.arange(97)+100,np.zeros(97)])]
    raw=np.r_[np.array([[10,10],[10,10],[0,0]]),np.column_stack([np.arange(97)+100,np.zeros(97)])].astype(int)
    time=np.r_[[.001,.002,.003],np.full(97,.003)];ids=np.arange(100);pol=np.zeros(100,int)
    _,meta=extract(xy,time,pol,raw,ids,[0,0],1000,'F1');record=next(r for r in meta['seed_diagnostics'] if r['seed_index']==2)
    assert 0 not in record['candidate_physical_ids'] and 1 not in record['candidate_physical_ids'] and 2 in record['candidate_physical_ids']
    tests.append('surface_latest_global_then_lane_radius_no_old_fallback')
    invalid=predict(np.array([[2.,0.]]),np.array([.005]),.006,np.zeros(2),np.zeros(2),1.,[0,0,0])
    assert invalid['status']=='geometry_invalid';tests.append('invalid_source_not_clamped_to_sphere')
    # The first V fragment maps the front pole exactly to the back pole; projected xy would match falsely.
    sources=np.zeros((10,2));src_t=np.full(10,.0055);sourcepol=np.zeros(10,int);sourceids=np.arange(10)
    target=np.zeros((10,2));vt=np.full(10,.00625);vp=np.zeros(10,int);vids=np.arange(100,110)
    result=future(sources,src_t,sourcepol,sourceids,target,vt,vp,vids,np.zeros(2),np.zeros(2),1.,[0,np.pi/.00075,0],[-5,-5,5,5])
    first=result['segments'][0];assert first['invisible_source_count']==10 and first['pred_to_target']['median']==3 and first['target_to_pred']['median']==3
    tests.append('backface_predictions_remain_three_pixel_penalties')
    # 100 coincident source points and only one matching target: count-pooled median would hide wrong reverse direction.
    source=np.zeros((100,2));st=np.full(100,.0055);sp=np.zeros(100,int);si=np.arange(100)
    obs=np.r_[np.zeros((1,2)),np.tile([3.,0.],(9,1))];ot=np.full(10,.00625)
    imbalance=future(source,st,sp,si,obs,ot,np.zeros(10,int),np.arange(200,210),np.zeros(2),np.zeros(2),10.,[0,0,0],[-5,-5,5,5])
    assert imbalance['statistics']['median']==0 and imbalance['target_to_pred']['median']==3
    assert np.isclose(imbalance['equal_direction_mean_chamfer'],1.35)
    baseline={**imbalance,'evaluable_segments':8,'pred_to_target':{**imbalance['pred_to_target'],'median':3.},'target_to_pred':{**imbalance['target_to_pred'],'median':3.}}
    model={**imbalance,'evaluable_segments':8};comparison=compare_future(model,baseline)
    assert not comparison['engineering_pass'];tests.append('two_directions_separate_despite_pooled_zero_median')
    same=compare_future(baseline,baseline);assert not same['engineering_pass']
    assert not accepted_prediction({'decision':{'measurement_valid':True}},{'engineering_pass':True},same)
    tests.append('translation_only_tie_cannot_support_nonzero_spin')
    assert not accepted_prediction({'decision':{'measurement_valid':False}},{'engineering_pass':True},{'engineering_pass':True})
    tests.append('rejected_LOO_candidate_cannot_count_as_success')
    report=dict(status='passed',tests=tests,scope='Small deterministic cause-order/prediction path controls; no extra real window or tuned threshold.')
    write_json(OUT/'path_verification.json',report);print(report)

if __name__=='__main__':main()
