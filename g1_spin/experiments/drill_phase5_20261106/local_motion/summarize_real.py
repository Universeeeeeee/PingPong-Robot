"""Failure attribution is descriptive, never changes frozen training constraints."""
import json
from pathlib import Path
from collections import Counter
import numpy as np

OUT=Path(__file__).resolve().parent
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False))
def dist(a):
    a=[v for v in a if v is not None];return dict(n=len(a),median=float(np.median(a)) if a else None,p95=float(np.percentile(a,95)) if a else None,max=float(max(a)) if a else None)

def main():
    real=OUT/'real';train=json.loads((real/'development_training_qualification.json').read_text());evaluation=json.loads((real/'development_evaluation.json').read_text());lookup={r['observation_id']:r for r in train};out=[]
    for case in ['run1_01495000','run2_02095000','run5_04095000','run5_01525000']:
        for rep in ['R_all','R_latest']:
            for kind in ['M_time','M_spatial']:
                rows=[r for r in evaluation if r['case']==case and r['representation']==rep and r['model_kind']==kind];serial=Counter();marginal=Counter();Tclasses=Counter();Hp=[];Fp=[];Hg=[];Fg=[]
                for r in rows:
                    t=lookup[r['observation_id']];q=t['T_qualification'];Tclasses[q['training_motion_class']]+=1
                    h=r['H'];f=r['shortfuture'];has_model=t['model']['available'];tq=q['training_qualified'];tn=q['training_motion_class']=='candidate_nonzero' and tq
                    Hsupport=bool(h and h['support_sufficient']);Fsupport=bool(f and f['support_sufficient']);Hd=bool(h and h['moving']['p95'] is not None and h['moving']['p95']<=.5);Fd=bool(f and f['moving']['p95'] is not None and f['moving']['p95']<=.5)
                    Hgain=bool(h and h['gain'] is not None and h['gain']>=.2);Fgain=bool(f and f['gain'] is not None and f['gain']>=.2);resolved=r.get('above_tau_motion_resolution',False);perm=r.get('permutation_pass',False)
                    gates=[('model_algebra',has_model),('T_training_qualified',tq),('T_candidate_nonzero',tn),('H_support',Hsupport),('future_support',Fsupport),('H_distance_P95',Hd),('future_distance_P95',Fd),('H_gain',Hgain),('future_gain',Fgain),('permutation',perm),('displacement_resolution',resolved)]
                    for name,passed in gates:
                        if passed:marginal[name]+=1
                    failed=next((name for name,passed in gates if not passed),None);serial[failed or 'supported_nonzero']+=1
                    if h:Hp.append(h['moving']['p95']);Hg.append(h['gain'])
                    if f:Fp.append(f['moving']['p95']);Fg.append(f['gain'])
                out.append(dict(case=case,representation=rep,model_kind=kind,records=len(rows),T_motion_classes=dict(Tclasses),
                    marginal_passes={name:marginal[name] for name in ['model_algebra','T_training_qualified','T_candidate_nonzero','H_support','future_support','H_distance_P95','future_distance_P95','H_gain','future_gain','permutation','displacement_resolution']},
                    serial_first_failure=dict(serial),H_P95_px=dist(Hp),future_P95_px=dist(Fp),H_gain=dist(Hg),future_gain=dist(Fg),
                    supported_nonzero=sum(r['nonzero_motion_evidence'] for r in rows),supported_spatial_patches=len(set(r['patch'] for r in rows if r['nonzero_motion_evidence'])),stationary=sum(r['stationary_constraint'] for r in rows)))
    write(real/'development_failure_attribution.json',out)
    systems=json.loads((real/'development_T_global_systems.json').read_text());qualified=[]
    for rep in ['R_all','R_latest']:
        for kind in ['M_time','M_spatial']:
            windows=[]
            for case in ['run1_01495000','run2_02095000']:
                evidence=next(r for r in out if r['case']==case and r['representation']==rep and r['model_kind']==kind);system=next(s for s in systems if s['case']==case and s['representation']==rep and s['model_kind']==kind)
                windows.append(dict(case=case,T_global_sufficient=system['global_training_sufficient'],supported_nonzero_patches=evidence['supported_spatial_patches'],qualified=False))
            qualified.append(dict(representation=rep,model_kind=kind,qualified=False,windows=windows))
    write(real/'selection_development.json',dict(status='none',selected_combination=None,combinations=qualified,
        reason='Neither rotation development window has>=3 spatial patches with heldout/permutation-supported nonzero motion; run2 also lacks12 T-qualified constraints for every combination.',
        H_future_never_filtered_global_constraints=True,negative_support_lack_not_success=True,confirmation_evaluated=False,regression108_run=False,
        prospective_action='Astra checkpoint decides fixed confirmation accounting; no new patch, parameter or representation selection.'))
    print('Saved marginal/serial failure attribution for16cells and selection=none')

if __name__=='__main__':main()
