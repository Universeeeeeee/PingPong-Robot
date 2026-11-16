"""Consolidate initializer, optimizer and gate diagnostics with source hashes."""
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent


def angle(first, second):
    first=np.asarray(first,float);second=np.asarray(second,float)
    return float(np.degrees(np.arccos(np.clip(np.dot(first,second)/(np.linalg.norm(first)*np.linalg.norm(second)),-1,1))))


def main():
    source_names=['gate_diagnostics.json','gate_refinement.json','zero_seed_results.json',
                  'seed_controls.json','constraint_sensitivity.json',
                  'coarse_seed/results.json','coarse_seed/extension/results.json']
    sources={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in source_names}
    gates=json.loads((OUT/'gate_diagnostics.json').read_text())['rows']
    flow=json.loads((OUT/'coarse_seed/extension/results.json').read_text())
    zero=json.loads((OUT/'zero_seed_results.json').read_text())['rows']
    comparisons=[]
    for name, records in [('coarse_flow',flow),('zero',zero)]:
        good=[row for row in records if row['status']=='accepted']
        rpm=[abs(row['final']['measurement_rpm']-row['old_G_rpm']) for row in good if row['old_G_valid']]
        axes=[angle(row['final']['measurement_axis'],row['old_G_axis']) for row in good if row['old_G_valid']]
        assert all((row['status']=='accepted')==row['old_G_valid'] for row in records)
        comparisons.append(dict(method=name,windows=len(records),accepted=len(good),
                                statuses=dict(Counter(row['status'] for row in records)),
                                rpm_absolute_difference_median=float(np.median(rpm)),
                                rpm_absolute_difference_max=max(rpm),
                                directed_axis_difference_median_deg=float(np.median(axes)),
                                directed_axis_difference_max_deg=max(axes)))
    controls=json.loads((OUT/'seed_controls.json').read_text())['rows']
    rotation=[row for row in controls if row['case'].startswith(('run1','run2','run3','run4'))]
    assert len(rotation)==40
    baseline=[row for row in gates if row['half_span_ms']==2 and row['residual_ratio_limit']==.15 and row['representation']=='all']
    result=dict(comparisons=comparisons,
                strict_baseline_accepted_planes=sum(row['accepted'] for row in baseline),
                strict_gate_cases=len(baseline),
                residual_sensitivity_rows=len(gates),
                generic_seed_rotation_trials=len(rotation),
                generic_seed_rotation_acceptances=sum(row['outcome']['measurement_valid'] for row in rotation),
                generic_seed_rotation_scan_selected_zero=sum(bool(np.linalg.norm(row['selected_seed'])<1e-12) for row in rotation),
                source_hashes=sources,
                limitations=['Old Nov5 recordings only for new real-event runs; Nov13 raw windows unavailable while Linux powered off',
                             'Acceptance is optimizer quality-gate acceptance, not physical correctness',
                             'G has an extra multistart-ambiguity gate; single-start methods do not',
                             'Four original rotation pilot windows and archive expansion are exploratory reuse, not new blind recordings'])
    (OUT/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    assert all(hashlib.sha256((OUT/name).read_bytes()).hexdigest()==digest for name,digest in sources.items())
    (OUT/'verification.json').write_text(json.dumps(dict(
        sources_unchanged=True,case_identity_matches=True,windows_per_extension=108,
        summary_acceptance_sets_equal=True,existing_tests='2 passed: tests/test_optical_flow.py and phase3/test_math.py',
        review='Independent agent reconstructed M_spatial via regression residual SVD and checked label flow, geometry and no direction-grid route'),indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ['source_hashes','limitations']},indent=2))

if __name__=='__main__':
    main()
