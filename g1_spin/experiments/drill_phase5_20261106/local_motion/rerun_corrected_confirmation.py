"""Rerun only corrected raw patch cases; run4 T/eval results hash-reused unchanged."""
import json
import numpy as np
from synthetic import OUT,write,digest
from models import representation,fit
from train_real import bootstrap
from qualify_training import classification,build_system
from evaluate_real import evaluate

def main():
    real=OUT/'real';prior=OUT/'prereview_clipped'
    if (real/'confirmation_training.json').exists():raise FileExistsError('Preserve corrected confirmation')
    correction=json.loads((OUT/'confirmation_input_correction.json').read_text());assert correction['status']=='passed'
    selection_hash=digest(real/'selection_development.json');assert selection_hash==correction['selection_none_sha256']
    manifest=json.loads((real/'confirmation_inputs.json').read_text());oldrows=json.loads((prior/'confirmation_training_qualification.json').read_text());rows=[];reuse=[]
    for case in manifest:
        if case['case']=='run4_02495000':
            old=[r for r in oldrows if r['case']==case['case']]
            assert all(r['observation_id']==len(rows)+i for i,r in enumerate(old));rows.extend(old);reuse.append(dict(case=case['case'],rows=len(old),all_model_bootstrap_values_reused=True));print('HASH REUSED T',case['case'],len(old),flush=True);continue
        for p in case['patches']:
            d=np.load(OUT/p['T_file']);xy=d['xy'];times=d['time_s'];pol=d['polarity'];ids=d['ids']
            for segment in range(3):
                low=segment*.002;high=low+.002
                for polarity in [0,1]:
                    index=np.flatnonzero((times>=low-1e-12)&(times<high-1e-12)&(pol==polarity))
                    for rep in ['R_all','R_latest']:
                        picked=index[representation(xy[index],times[index],ids[index],rep)];identity=f'{case["case"]}/patch{p["patch"]}/segment{segment}/pol{polarity}/{rep}'
                        for kind in ['M_time','M_spatial']:
                            model=fit(xy[picked],times[picked],kind);boot=bootstrap(xy[picked],times[picked],kind,identity,model)
                            r=dict(observation_id=len(rows),case=case['case'],patch=p['patch'],box=p['box'],segment=segment,polarity=polarity,representation=rep,model_kind=kind,
                                start_s=low,end_s=high,T_file=p['T_file'],T_sha256=p['T_sha256'],all_segment_T_ids=ids[index].tolist(),support_ids=ids[picked].tolist(),model=model,bootstrap=boot)
                            r['T_qualification']=classification(r);rows.append(r)
        print('CORRECTED T TRAIN',case['case'],flush=True)
    assert len(rows)==384
    systems=[];oldsystems=json.loads((prior/'confirmation_T_global_systems.json').read_text())
    for case in manifest:
        for rep in ['R_all','R_latest']:
            for kind in ['M_time','M_spatial']:
                if case['case']=='run4_02495000':systems.append(next(s for s in oldsystems if s['case']==case['case'] and s['representation']==rep and s['model_kind']==kind));continue
                local=[r for r in rows if r['case']==case['case'] and r['representation']==rep and r['model_kind']==kind]
                systems.append(dict(case=case['case'],representation=rep,model_kind=kind,**build_system(local,np.array(case['center']),case['radius_px'])))
    write(real/'confirmation_training.json',rows);write(real/'confirmation_training_qualification.json',rows);write(real/'confirmation_T_global_systems.json',systems)
    write(real/'confirmation_training_freeze.json',dict(stage='Corrected complete decoded raw patches T-only, before corrected H/future; run4 reused',records=len(rows),
        files={r:digest(OUT/r) for r in ['repair_confirmation_inputs.py','rerun_corrected_confirmation.py','confirmation_input_correction.json','models.py','train_real.py','qualify_training.py','real/confirmation_inputs.json','real/confirmation_training.json','real/confirmation_training_qualification.json','real/confirmation_T_global_systems.json']},
        selection_development_sha256=selection_hash,H_future_evaluation_performed=False,reused_run4=reuse))
    config=json.loads((OUT/'confirmation_evaluation_freeze.json').read_text());config['stage']='Corrected raw patch T freeze before corrected evaluation; same tau/bootstrap/permutation/metrics'
    config['files']={r:digest(OUT/r) for r in ['tau_freeze.json','evaluate_real.py','models.py','qualify_training.py','real/confirmation_training_freeze.json','real/confirmation_training_qualification.json','real/confirmation_T_global_systems.json']}
    write(OUT/'confirmation_corrected_evaluation_freeze.json',config)
    # Frozen evaluator supports case subsets; run4 evaluation is reused exactly, not re-scored.
    evaluate('confirmation',case_filter=['run3_01495000','run5_05095000'],configuration='confirmation_corrected_evaluation_freeze.json')
    new=json.loads((real/'confirmation_evaluation.json').read_text());old=json.loads((prior/'confirmation_evaluation.json').read_text());run4=[r for r in old if r['case']=='run4_02495000']
    allrows=sorted(new+run4,key=lambda r:r['observation_id']);assert len(allrows)==384
    write(real/'confirmation_evaluation.json',allrows);write(real/'confirmation_evaluation_summary.json',dict(records=384,nonzero=sum(r['nonzero_motion_evidence'] for r in allrows),stationary=sum(r['stationary_constraint'] for r in allrows),
        run4_evaluation_rows_hash_reused=len(run4),raw_patch_input_correction=True,old_clipped_results='prereview_clipped/'))
    assert digest(real/'selection_development.json')==selection_hash
    print('CORRECTED CONFIRMATION COMPLETE',len(allrows),'nonzero',sum(r['nonzero_motion_evidence'] for r in allrows),'stationary',sum(r['stationary_constraint'] for r in allrows),flush=True)

if __name__=='__main__':main()
