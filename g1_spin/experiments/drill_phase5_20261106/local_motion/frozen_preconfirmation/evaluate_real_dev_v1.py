"""Evaluate frozen T models; heldouts never select constraints for the global systems."""
import hashlib
import json
from pathlib import Path
import numpy as np
from models import representation,fit,residual,equal_unit_distances,canonical

OUT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False))

def static_line(xy,times):
    if len(xy)<2:return None
    center=xy.mean(axis=0);X=xy-center;values,vectors=np.linalg.eigh(X.T@X/len(X));n,speed,_=canonical(vectors[:,0],0.,0.)
    return dict(available=True,normal=n.tolist(),speed_px_s=speed,center=center.tolist(),reference_s=float(times.mean()),intercept_px=0.,
        eigenvalues=values.tolist(),direction_unique=bool(values[1]-values[0]>1e-8*max(1.,values.sum())),
        role='Best stationary spatial line on exact same T representation; no H/future optimization')

def scores(model,static,xy,times,pol,is_future,tau):
    moving=equal_unit_distances(xy,times,pol,residual(model,xy,times));fixed=equal_unit_distances(xy,times,pol,residual(static,xy,times))
    support=bool(moving['unique_pixels']>=4 and (not is_future or moving['time_bins']>=2))
    gain=None if moving['mean'] is None else (fixed['mean']-moving['mean'])/max(fixed['mean'],tau)
    return dict(moving=moving,static=fixed,support_sufficient=support,gain=gain,
        distance_and_gain_pass=bool(support and moving['p95']<=tau and gain>=.2))

def permutations(Txy,Ttime,Tids,Hxy,Htime,Hpol,Fxy,Ftime,Fpol,rep,kind,identity,tau):
    seed=int.from_bytes(hashlib.sha256(identity.encode()).digest()[:8],'little');rng=np.random.default_rng(seed);rows=[]
    for draw in range(19):
        tt=Ttime[rng.permutation(len(Ttime))];ht=Htime[rng.permutation(len(Htime))];ft=Ftime[rng.permutation(len(Ftime))]
        chosen=representation(Txy,tt,Tids,rep);model=fit(Txy[chosen],tt[chosen],kind);stat=None;h=None;f=None
        if model['available']:
            fixed=static_line(Txy[chosen],tt[chosen])
            if fixed is not None:
                h=scores(model,fixed,Hxy,ht,Hpol,False,tau);f=scores(model,fixed,Fxy,ft,Fpol,True,tau)
                if h['support_sufficient'] and f['support_sufficient']:stat=min(h['gain'],f['gain'])
        rows.append(dict(draw=draw,model_available=model['available'],support_ids=Tids[chosen].tolist(),statistic=stat,
            H_gain=None if h is None else h['gain'],future_gain=None if f is None else f['gain']))
    return dict(seed=seed,draws=19,valid_statistics=sum(r['statistic'] is not None for r in rows),records=rows,
        shuffle_boundary='T,H,next1ms future separately; each observation already same polarity; xy/event sets unchanged, R_latest reselected on shuffled T')

def evaluate(partition='development'):
    real=OUT/'real';out=real/f'{partition}_evaluation.json'
    if out.exists():raise FileExistsError('Preserve heldout evaluation')
    config=json.loads((OUT/'real_evaluation_freeze.json').read_text())
    for rel,h in config['files'].items():assert digest(OUT/rel)==h,rel
    manifest=json.loads((real/f'{partition}_inputs.json').read_text());rows=json.loads((real/f'{partition}_training_qualification.json').read_text());result=[]
    for r in rows:
        case=next(c for c in manifest if c['case']==r['case']);patch=next(p for p in case['patches'] if p['patch']==r['patch']);train=np.load(OUT/patch['T_file']);held=np.load(OUT/patch['evaluation_file'])
        model=r['model'];low=r['start_s'];high=r['end_s'];pol=r['polarity'];tau=config['tau_px']
        TM=(train['time_s']>=low-1e-12)&(train['time_s']<high-1e-12)&(train['polarity']==pol)
        HM=(held['H_time_s']>=low-1e-12)&(held['H_time_s']<high-1e-12)&(held['H_polarity']==pol)
        if r['segment']<2:
            FM=(held['H_time_s']>=high-1e-12)&(held['H_time_s']<high+.001-1e-12)&(held['H_polarity']==pol);prefix='H'
        else:
            FM=(held['V_time_s']>=high-1e-12)&(held['V_time_s']<high+.001-1e-12)&(held['V_polarity']==pol);prefix='V'
        Hxy,Htime,Hpol=held['H_xy'][HM],held['H_time_s'][HM],held['H_polarity'][HM]
        Fxy,Ftime,Fpol=held[prefix+'_xy'][FM],held[prefix+'_time_s'][FM],held[prefix+'_polarity'][FM]
        record=dict(observation_id=r['observation_id'],case=r['case'],patch=r['patch'],segment=r['segment'],polarity=pol,representation=r['representation'],model_kind=r['model_kind'],
            T_qualification=r['T_qualification'],H_ids=held['H_ids'][HM].tolist(),shortfuture_ids=held[prefix+'_ids'][FM].tolist(),shortfuture_source=prefix,
            H=None,shortfuture=None,permutation=None,nonzero_motion_evidence=False,stationary_constraint=False,status='model_unidentifiable')
        if model['available']:
            lookup={int(v):i for i,v in enumerate(train['ids'])};selected=np.array([lookup[i] for i in r['support_ids']],int);fixed=static_line(train['xy'][selected],train['time_s'][selected])
            h=scores(model,fixed,Hxy,Htime,Hpol,False,tau);f=scores(model,fixed,Fxy,Ftime,Fpol,True,tau)
            identity=f'{r["case"]}/{r["patch"]}/{r["segment"]}/{pol}/{r["representation"]}'
            perm=permutations(train['xy'][TM],train['time_s'][TM],train['ids'][TM],Hxy,Htime,Hpol,Fxy,Ftime,Fpol,r['representation'],r['model_kind'],identity,tau)
            valid=[x['statistic'] for x in perm['records'] if x['statistic'] is not None]
            statistic=min(h['gain'],f['gain']) if h['support_sufficient'] and f['support_sufficient'] else None
            perm_pass=bool(statistic is not None and len(valid)>=15 and statistic>max(valid));resolved=abs(model['speed_px_s'])*.001>=tau
            q=r['T_qualification'];candidate=q['training_qualified'] and q['training_motion_class']=='candidate_nonzero'
            nonzero=bool(candidate and resolved and h['distance_and_gain_pass'] and f['distance_and_gain_pass'] and perm_pass)
            stationary=bool(q['training_qualified'] and q['training_motion_class']=='static_compatible' and h['support_sufficient'] and f['support_sufficient'] and h['static']['p95']<=tau and f['static']['p95']<=tau)
            status='supported_nonzero' if nonzero else ('stationary_constraint' if stationary else ('insufficient_heldout_support' if not h['support_sufficient'] or not f['support_sufficient'] else 'not_supported_by_fixed_metric'))
            record.update(H=h,shortfuture=f,static_training_model=fixed,permutation=perm,actual_statistic=statistic,permutation_pass=perm_pass,
                predicted_displacement_1ms_px=abs(model['speed_px_s'])*.001,above_tau_motion_resolution=bool(resolved),
                nonzero_motion_evidence=nonzero,stationary_constraint=stationary,status=status)
        result.append(record)
    write(out,result);write(real/f'{partition}_evaluation_summary.json',dict(records=len(result),nonzero=sum(r['nonzero_motion_evidence'] for r in result),
        stationary=sum(r['stationary_constraint'] for r in result),all_constraints_frozen_before_evaluation=True,
        note='H/future/permutation evaluate evidence only; never remove or refit T constraints before global systems.'))
    print('EVALUATED',partition,len(result),'nonzero',sum(r['nonzero_motion_evidence'] for r in result),'stationary',sum(r['stationary_constraint'] for r in result))

if __name__=='__main__':evaluate()
