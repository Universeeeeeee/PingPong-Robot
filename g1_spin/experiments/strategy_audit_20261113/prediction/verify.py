"""Independent lift, Rodrigues projection, brute-force distances and reversal algebra."""
import importlib.util,json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('experiment',OUT/'run.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
results=json.loads((OUT/'results.json').read_text());files={name:(path,start,duration/1e6) for name,path,start,duration in r.FILES}

def rotate(points,dt,omega):
    omega=np.array(omega);speed=np.linalg.norm(omega)
    if speed<1e-15:return points.copy()
    axis=omega/speed;theta=np.asarray(dt)*speed;return points*np.cos(theta)[:,None]+np.cross(axis,points)*np.sin(theta)[:,None]+np.outer(points@axis,axis)*(1-np.cos(theta))[:,None]


def brute_nearest(a,b):
    if not len(a):return np.empty(0)
    if not len(b):return np.full(len(a),3.)
    return np.concatenate([np.minimum(cdist(a[i:i+64],b).min(axis=1),3) for i in range(0,len(a),64)])
checks=[];summary=[]
for item in results:
    if item['case'].startswith('synthetic'):
        data,_=r.synthetic([0,1,0] if item['case']=='synthetic_y' else [.6,.7,.4]);xy,t,pol,cs,radius=data;duration=.05
    else:
        path,start,duration=files[item['case']];d=np.load(path)
        if 'events' in d:
            e=d['events'];xy=np.column_stack([e['x'],e['y']]).astype(float);t=(e['t'].astype(np.int64)-int(d['origin_us'])-start)*1e-6;pol=e['p'];cs=d['linear50'];radius=48.
        else:
            xy=d['xy'];t=(d['time_us']-start)*1e-6;pol=d['polarity'];cs=np.broadcast_to(d['center'],xy.shape);radius=float(d['radius_px'])
    normalized=(xy-cs)/radius;pts=np.column_stack([normalized,-np.sqrt(1-np.sum(normalized**2,axis=1))]);cut=.6*duration;train=t<cut;source=train&(t>=cut-min(.005,cut/3))
    assert train.sum()==item['train_events'] and (~train).sum()==item['holdout_events'];assert not np.any(source&~train)
    estimates={f['kind']:f['estimate'] for f in item['fits']};cand={k:np.array(v['omega_rad_s'] if v['valid'] else v['quality'].get('rejected_candidate_omega_rad_s')) if v['valid'] or v['quality'].get('rejected_candidate_omega_rad_s') is not None else None for k,v in estimates.items()};cand['zero']=np.zeros(3);cand['opposite']=-cand['real'] if cand['real'] is not None else None
    max_difference=0
    for kind,stored in item['evaluation'].items():
        if stored is None:continue
        omega=cand[kind];scores=[]
        for j in range(4):
            lo=cut+j*duration*.1;hi=cut+(j+1)*duration*.1;mid=(lo+hi)/2;target=(t>=lo)&(t<hi)
            if source.sum()<10 or target.sum()<10:continue
            future=rotate(pts[source],mid-t[source],omega);pred=future[:,:2]*radius+np.array([np.interp(mid,t,cs[:,dim]) for dim in range(2)]);visible=future[:,2]<0;f=[];b=[]
            for p in [0,1]:
                q=pol[source]==p;actual=xy[target][pol[target]==p];seen=visible[q];prediction=pred[q];forward=np.full(len(prediction),3.);forward[seen]=brute_nearest(prediction[seen],actual);back=brute_nearest(actual,prediction[seen]);f.extend(forward);b.extend(back)
            scores.append((np.mean(f)+np.mean(b))/2)
        actual=float(np.mean(scores));difference=abs(actual-stored['future_balanced_mean']);assert difference<1e-9,(item['case'],kind,difference);max_difference=max(max_difference,difference)
    # Correct time reversal requires reversing timestamps as well as point order.
    omega=cand['real'];reversal_error=None
    if omega is not None:
        tt=t[train];pp=pts[train];ref=(tt[0]+tt[-1])/2;revtime=tt[0]+tt[-1]-tt[::-1]
        forward=rotate(pp,-(tt-ref),omega);reverse=rotate(pp[::-1],-(revtime-ref),-omega)
        reversal_error=float(np.max(abs(forward[::-1]-reverse)));assert reversal_error<1e-12
    evaluation=item['evaluation'];base=evaluation['real'];competitors=[evaluation[k] for k in ['zero','opposite','shuffle101','shuffle102','shuffle103']]
    native=estimates['real']['valid'];prediction_win=bool(base and base['evaluable_slices']>=3 and all(v and v['evaluable_slices']>=3 and base['future_balanced_mean']<v['future_balanced_mean'] for v in competitors));contrast_win=bool(base and all(v and base['iwe_holdout_cosine']>v['iwe_holdout_cosine'] for v in competitors))
    gain=1-base['future_balanced_mean']/evaluation['zero']['future_balanced_mean'] if base else None
    summary.append(dict(case=item['case'],native_valid=native,rpm=estimates['real']['spin_rpm'],reason=estimates['real']['failure_reason'],future_distance_px=base['future_balanced_mean'] if base else None,zero_distance_px=evaluation['zero']['future_balanced_mean'],future_gain_vs_zero=gain,prediction_beats_all_controls=prediction_win,contrast_beats_all_controls=contrast_win,predeclared_joint_pass=native and prediction_win and contrast_win,shuffle_native_valid=sum(estimates[k]['valid'] for k in ['shuffle101','shuffle102','shuffle103'])))
    checks.append(dict(case=item['case'],train_holdout_disjoint=True,brute_force_metric_max_error=max_difference,correct_time_reversal_point_error=reversal_error))
for filename,h in json.loads((OUT/'freeze.json').read_text()).items():assert r.digest(r.ROOT/filename)==h
r.save('verification.json',dict(checks=checks,files_unchanged=True,discarded_diagnostic='reversed_training in original run only reversed points at nonuniform fixed times; it is not exact time reversal, excluded from conclusions. Correct event-time reversal algebra independently checked here.'))
r.save('summary.json',summary)
print(json.dumps(summary,indent=2))
