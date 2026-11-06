"""Bounded incremental confirmation verification and failure summary; no model refits."""
import json
from collections import Counter
import numpy as np
from synthetic import OUT,P4,digest,write
from models import residual
from verify_real import equal_mean_p95

def checked_freeze(name):
    for rel,h in json.loads((OUT/name).read_text())['files'].items():
        if digest(OUT/rel)==h:continue
        if rel=='qualify_training.py':amend=json.loads((OUT/'geometry_path_amendment.json').read_text())
        elif rel=='evaluate_real.py':amend=json.loads((OUT/'confirmation_evaluator_hash_amendment.json').read_text())
        else:raise AssertionError((name,rel))
        assert amend['old_source_sha256']==h and digest(OUT/amend['old_source'])==h
        if rel=='evaluate_real.py' and amend['new_source_sha256']!=digest(OUT/rel):
            correction=json.loads((OUT/'confirmation_input_code_amendment.json').read_text())
            assert correction['old_source_sha256']==amend['new_source_sha256'] and correction['new_source_sha256']==digest(OUT/rel)
            assert digest(OUT/correction['old_evaluator_source'])==correction['old_source_sha256']
        else:assert amend['new_source_sha256']==digest(OUT/rel)

def quantiles(a):
    a=[v for v in a if v is not None];return dict(n=len(a),median=float(np.median(a)) if a else None,p95=float(np.percentile(a,95)) if a else None,max=float(max(a)) if a else None)

def main():
    real=OUT/'real';train=json.loads((real/'confirmation_training_qualification.json').read_text());evaluation=json.loads((real/'confirmation_evaluation.json').read_text());manifest=json.loads((real/'confirmation_inputs.json').read_text());systems=json.loads((real/'confirmation_T_global_systems.json').read_text())
    assert len(train)==len(evaluation)==384 and len(systems)==12
    for name in ['confirmation_freeze.json','real/confirmation_training_freeze.json','confirmation_corrected_evaluation_freeze.json','development_selection_freeze.json']:checked_freeze(name)
    old_sources=json.loads((P4/'real_frontend/input_manifest.json').read_text());byid={r['observation_id']:r for r in train};support={};metric_checks=0
    for case in manifest:
        old=next(r for r in old_sources if r['case']==case['case']);assert digest(case['source_file'])==case['source_sha256']==old['input_sha256']
        assert case['radius_px']==old['radius_px'] and case['center']==old['center']
    for e in evaluation:
        r=byid[e['observation_id']];case=next(c for c in manifest if c['case']==r['case']);patch=next(p for p in case['patches'] if p['patch']==r['patch']);T=np.load(OUT/patch['T_file']);held=np.load(OUT/patch['evaluation_file'])
        assert digest(OUT/patch['T_file'])==patch['T_sha256'] and digest(OUT/patch['evaluation_file'])==patch['evaluation_sha256']
        tset=set(T['ids'].tolist());hset=set(held['H_ids'].tolist());vset=set(held['V_ids'].tolist())
        assert not tset&hset and not tset&vset and not hset&vset and set(r['support_ids'])<=tset
        assert np.all(T['time_s']<.006) and np.all((T['xy'][:,0]+2*T['xy'][:,1])%5!=0)
        assert set(e['H_ids'])<=hset and set(e['shortfuture_ids'])<=(hset if r['segment']<2 else vset)
        key=(r['case'],r['patch'],r['segment'],r['polarity'],r['representation']);support.setdefault(key,[]).append(r['support_ids'])
        assert r['bootstrap']['draws']==32
        if e['H'] is not None:
            for field,prefix in [('H','H'),('shortfuture','H' if r['segment']<2 else 'V')]:
                ids=e['H_ids'] if field=='H' else e['shortfuture_ids'];where={int(x):i for i,x in enumerate(held[prefix+'_ids'])};ix=np.array([where[x] for x in ids],int)
                xy=held[prefix+'_xy'][ix];times=held[prefix+'_time_s'][ix];pol=held[prefix+'_polarity'][ix]
                for score,model in [('moving',r['model']),('static',e['static_training_model'])]:
                    mean,p95=equal_mean_p95(xy,times,pol,residual(model,xy,times));saved=e[field][score]
                    assert (mean is None and saved['mean'] is None) or np.isclose(mean,saved['mean'],atol=1e-12)
                    assert (p95 is None and saved['p95'] is None) or np.isclose(p95,saved['p95'],atol=1e-12);metric_checks+=1
            assert len(e['permutation']['records'])==19 and e['permutation']['draws']==19
    assert all(len(v)==2 and v[0]==v[1] for v in support.values())
    for s in systems:
        expected=[r['observation_id'] for r in train if r['case']==s['case'] and r['representation']==s['representation'] and r['model_kind']==s['model_kind'] and r['T_qualification']['training_qualified']]
        assert expected==s['observation_ids']
        for leave in s['leave_patch_out']:assert all(byid[i]['patch']!=leave['heldout_patch'] for i in leave['training_observation_ids'])
    groups=[]
    for case in manifest:
        for rep in ['R_all','R_latest']:
            for kind in ['M_time','M_spatial']:
                rows=[r for r in evaluation if r['case']==case['case'] and r['representation']==rep and r['model_kind']==kind];trs=[byid[r['observation_id']] for r in rows];serial=Counter();marginal=Counter()
                for e,r in zip(rows,trs):
                    q=r['T_qualification'];h=e['H'];f=e['shortfuture'];gates=[('model_algebra',r['model']['available']),('T_training_qualified',q['training_qualified']),('T_candidate_nonzero',q['training_qualified'] and q['training_motion_class']=='candidate_nonzero'),
                        ('H_support',bool(h and h['support_sufficient'])),('future_support',bool(f and f['support_sufficient'])),('H_P95',bool(h and h['moving']['p95'] is not None and h['moving']['p95']<=.5)),('future_P95',bool(f and f['moving']['p95'] is not None and f['moving']['p95']<=.5)),
                        ('H_gain',bool(h and h['gain'] is not None and h['gain']>=.2)),('future_gain',bool(f and f['gain'] is not None and f['gain']>=.2)),('permutation',e.get('permutation_pass',False)),('motion_resolution',e.get('above_tau_motion_resolution',False))]
                    for name,passed in gates:
                        if passed:marginal[name]+=1
                    serial[next((name for name,passed in gates if not passed),'supported_nonzero')]+=1
                system=next(s for s in systems if s['case']==case['case'] and s['representation']==rep and s['model_kind']==kind)
                groups.append(dict(case=case['case'],representation=rep,model_kind=kind,records=len(rows),T_qualified=system['T_qualified_count'],T_patches=system['independent_spatial_patches'],global_T_sufficient=system['global_training_sufficient'],
                    T_motion_classes=dict(Counter(r['T_qualification']['training_motion_class'] for r in trs)),marginal_passes=dict(marginal),serial_first_failure=dict(serial),
                    H_P95_px=quantiles([r['H']['moving']['p95'] for r in rows if r['H']]),future_P95_px=quantiles([r['shortfuture']['moving']['p95'] for r in rows if r['shortfuture']]),
                    nonzero=sum(r['nonzero_motion_evidence'] for r in rows),stationary=sum(r['stationary_constraint'] for r in rows),supported_nonzero_patches=len(set(r['patch'] for r in rows if r['nonzero_motion_evidence']))))
    write(real/'confirmation_failure_attribution.json',groups)
    assert not any(r['nonzero_motion_evidence'] for r in evaluation) and not any(r['stationary_constraint'] for r in evaluation)
    assert json.loads((real/'selection_development.json').read_text())['selected_combination'] is None
    correction=json.loads((OUT/'confirmation_input_correction.json').read_text());assert correction['status']=='passed'
    oldtrain=json.loads((OUT/'prereview_clipped/confirmation_training_qualification.json').read_text());oldeval=json.loads((OUT/'prereview_clipped/confirmation_evaluation.json').read_text())
    assert [r for r in train if r['case']=='run4_02495000']==[r for r in oldtrain if r['case']=='run4_02495000']
    assert [r for r in evaluation if r['case']=='run4_02495000']==[r for r in oldeval if r['case']=='run4_02495000']
    assert digest(real/'selection_development.json')==correction['selection_none_sha256']
    verification=dict(status='passed',date='2026-11-07',confirmation_windows=3,selected_patches=[r['patch_count'] for r in manifest],confirmation_local_model_records=384,
        development_records=360,total_real_model_records=744,same_representation_2fits_same_support=True,T_H_V_ids_disjoint=True,
        independent_equal_unit_mean_P95_checks=metric_checks,bootstrap_draws=32,timepermutation_draws=19,global_constraints_T_only=True,
        prior_dev_freezes_and_selection_none_unchanged=True,code_hash_plumbing_amendment='confirmation_evaluator_hash_amendment.json',
        complete_decoded_raw_patch_correction='confirmation_input_correction.json',clipped_confirmation_preserved='prereview_clipped/',
        corrected_run3_run5_rerun_only=True,run4_training_and_evaluation_reused_exactly=True,
        nonzero_motion_evidence=0,stationary_evidence=0,no_CMax_run=True,no108_run=True,source_geometry_hashes_preserved=True,
        stopping_reason='No>=3patch heldout-supported nonzero motion in any confirmation combination; frozen development none remains, no new fit or combination selection.')
    write(real/'verification_confirmation.json',verification);print(verification)

if __name__=='__main__':main()
