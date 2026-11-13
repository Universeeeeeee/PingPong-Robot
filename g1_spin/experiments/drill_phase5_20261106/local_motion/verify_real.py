"""Incremental dev checkpoint verification; no new real models or thresholds."""
import json
from pathlib import Path
import numpy as np
from synthetic import OUT,P4,digest,write
from models import residual

def equal_mean_p95(xy,times,pol,distance):
    groups={}
    for x,t,p,d in zip(xy,times,pol,distance):groups.setdefault((int(x[0]),int(x[1]),int(p),int(np.floor(t/.0005+1e-10))),[]).append(float(d))
    values=[sum(v)/len(v) for v in groups.values()]
    return (float(np.mean(values)),float(np.percentile(values,95))) if values else (None,None)

def main():
    real=OUT/'real';manifest=json.loads((real/'development_inputs.json').read_text());train=json.loads((real/'development_training_qualification.json').read_text());evaluation=json.loads((real/'development_evaluation.json').read_text());systems=json.loads((real/'development_T_global_systems.json').read_text())
    assert len(train)==len(evaluation)==360 and len(systems)==16
    for name in ['pre_real_freeze.json','real/development_training_freeze.json','real/development_qualification_freeze.json','real_evaluation_freeze.json']:
        for rel,h in json.loads((OUT/name).read_text())['files'].items():
            if rel=='qualify_training.py' and digest(OUT/rel)!=h:
                amendment=json.loads((OUT/'geometry_path_amendment.json').read_text())
                assert amendment['status']=='passed' and amendment['old_source_sha256']==h and amendment['new_source_sha256']==digest(OUT/rel)
                assert digest(OUT/amendment['old_source'])==h
            else:assert digest(OUT/rel)==h,rel
    lookup={r['observation_id']:r for r in train};metrics=0;constraints=0;support={}
    for e in evaluation:
        r=lookup[e['observation_id']];case=next(c for c in manifest if c['case']==r['case']);patch=next(p for p in case['patches'] if p['patch']==r['patch']);T=np.load(OUT/patch['T_file']);held=np.load(OUT/patch['evaluation_file']);Tids=set(T['ids'].tolist());Hids=set(held['H_ids'].tolist());Vids=set(held['V_ids'].tolist())
        assert digest(OUT/patch['T_file'])==patch['T_sha256'] and digest(OUT/patch['evaluation_file'])==patch['evaluation_sha256']
        assert set(r['support_ids'])<=Tids and not Tids&Hids and not Tids&Vids
        assert set(e['H_ids'])<=Hids and set(e['shortfuture_ids'])<=(Hids if r['segment']<2 else Vids)
        assert set(r['support_ids']).isdisjoint(e['H_ids']) and set(r['support_ids']).isdisjoint(e['shortfuture_ids'])
        key=(r['case'],r['patch'],r['segment'],r['polarity'],r['representation']);support.setdefault(key,[]).append(r['support_ids'])
        bootstrap=r['bootstrap'];assert bootstrap['draws']==32
        if r['model']['available']:
            assert len(bootstrap['samples'])==32
            valid=[x for x in bootstrap['samples'] if x['available']];assert len(valid)==bootstrap['valid_draws']
        if e['H'] is not None:
            for target,source_prefix in [('H','H'),('shortfuture','H' if r['segment']<2 else 'V')]:
                score=e[target];ids=e['H_ids'] if target=='H' else e['shortfuture_ids'];where={int(v):i for i,v in enumerate(held[source_prefix+'_ids'])};ix=np.array([where[i] for i in ids],int)
                xy=held[source_prefix+'_xy'][ix];times=held[source_prefix+'_time_s'][ix];pol=held[source_prefix+'_polarity'][ix]
                for modelname,output in [('moving',r['model']),('static',e['static_training_model'])]:
                    m,p=equal_mean_p95(xy,times,pol,residual(output,xy,times));assert m==score[modelname]['mean'] or (m is not None and np.isclose(m,score[modelname]['mean'],atol=1e-12));assert p==score[modelname]['p95'] or (p is not None and np.isclose(p,score[modelname]['p95'],atol=1e-12));metrics+=1
            assert e['permutation']['draws']==19 and len(e['permutation']['records'])==19
    assert all(len(x)==2 and x[0]==x[1] for x in support.values())
    for s in systems:
        expected=[r['observation_id'] for r in train if r['case']==s['case'] and r['representation']==s['representation'] and r['model_kind']==s['model_kind'] and r['T_qualification']['training_qualified']]
        assert expected==s['observation_ids'];constraints+=len(expected)
        for leave in s['leave_patch_out']:
            assert all(lookup[i]['patch']!=leave['heldout_patch'] for i in leave['training_observation_ids'])
    prediction=json.loads((real/'development_global_LOPO_predictions.json').read_text())
    for s in prediction:
        frozen=next(g for g in systems if g['case']==s['case'] and g['representation']==s['representation'] and g['model_kind']==s['model_kind'])
        for patch in s['patches']:
            expected=next(l for l in frozen['leave_patch_out'] if l['heldout_patch']==patch['heldout_patch']);assert patch['training_observation_ids']==expected['training_observation_ids']
            assert not patch['accepted_rotation_explanation']
            assert all(lookup[i]['patch']!=patch['heldout_patch'] for i in patch['training_observation_ids'])
    assert not any(r['nonzero_motion_evidence'] for r in evaluation)
    selection=json.loads((real/'selection_development.json').read_text());assert selection['selected_combination'] is None
    for rel,h in json.loads((OUT.parent/'manifest_raw_structure.json').read_text())['files'].items():assert digest(OUT.parent/rel)==h['sha256'],rel
    for rel,h in json.loads((P4/'manifest_task13.json').read_text())['files'].items():assert digest(P4/rel)==h['sha256'],rel
    report=dict(status='passed',date='2026-11-07',dev_model_observations=360,global_combination_systems=16,
        T_only_selected_global_constraint_rows=constraints,H_future_never_prune_global_systems=True,same_representation_2fits_same_support=True,
        independent_equal_unit_mean_and_P95_checks=metrics,T_H_V_physical_IDs_disjoint=True,
        actual_bootstrap_draws=32,actual_timepermutation_draws=19,global_LOPO_heldpatch_removed=True,
        accepted_global_rotation_claims=0,nonzero_local_support=0,stationary_local_support=0,selection_none=True,
        pre_real_T_training_qualification_and_evaluation_outputs_unchanged=True,old_dev_code_freeze_preserved=True,
        future_geometry_invalid_branch_amendment='geometry_path_amendment.json',raw_checkpoint_and_stage4_task13_unchanged=True,
        confirmation_not_evaluated=True,no_CMax_or108_run=True,scope='Developer local2x2 diagnostic checkpoint; not real accuracy or negative-control success.')
    write(real/'verification_development.json',report);print(report)

if __name__=='__main__':main()
