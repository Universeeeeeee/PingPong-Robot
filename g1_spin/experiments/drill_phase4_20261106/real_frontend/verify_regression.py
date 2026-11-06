"""Incremental final checks; reuses already verified frontend and prediction interfaces."""
import gzip
import json
import numpy as np
from helpers import OUT,ROOT,P3,PHASE4,digest,hash_arrays,write_json
from pipeline import split

def read(path):return json.loads(gzip.decompress(path.read_bytes()))

def main():
    p=OUT/'regression108';summary=json.loads((p/'summary.json').read_text());states=json.loads((p/'states648.json').read_text());g=read(p/'G_reference.json.gz')
    assert len(summary)==324 and len(states)==648 and len(g)==108
    assert len({(r['case'],r['frontend'],r['method']) for r in states})==648
    assert len({r['case'] for r in g})==108
    known={r['case']:r for r in json.loads((OUT/'input_manifest.json').read_text())}
    audit={r['case']:r for r in json.loads((p/'input_audit.json').read_text())};old={r['window']:r for r in json.loads((P3/'window_audit.json').read_text())}
    reused=json.loads((p/'reused_eight_windows.json').read_text());assert len(reused)==24 and len({r['case'] for r in reused})==8
    reusechecks=[];splitchecks=[];sources={}
    for name,r in audit.items():
        source=ROOT/r['source_input'];assert digest(source)==r['source_input_sha256'];sources[r['source_input']]=r['source_input_sha256']
    for name in sorted({r['case'] for r in reused}):
        r=known[name];source=ROOT/r['source_window_input'];d=np.load(source);saved=np.load(OUT/r['input_file'])
        assert digest(source)==r['source_input_sha256']==audit[name]['source_input_sha256']
        assert float(saved['radius'])==float(old[name]['radius_px'])==r['radius_px']
        assert np.array_equal(d['common_ids'],saved['ids']) and np.array_equal(d['center'],saved['center'])
        assert digest(OUT/r['input_file'])==r['input_sha256'] and r['input_hash']==audit[name]['input_hash']
        reusechecks.append(dict(case=name,radius_px=r['radius_px'],source_sha256=r['source_input_sha256'],passed=True))
    for entry in summary:
        if entry.get('status')=='geometry_unavailable':
            assert old[entry['case']]['geometry_status']=='unavailable'
            assert all(r['failure_reason']=='geometry_unavailable' and not r['measurement_valid'] for r in entry['routes'].values());continue
        detail=OUT/entry['detail'];assert digest(detail)==entry['detail_sha256'];d=read(detail)
        T=set(d['T_ids']);H=set(d['H_ids']);V=set(d['V_ids']);assert not T&H and not T&V and not H&V
        source=np.load(ROOT/audit[entry['case']]['source_input']);assert T|H|V==set(source['common_ids'].tolist())
        if entry.get('reused'):
            b=np.load(OUT/known[entry['case']]['input_file']);m=split(b['xy'],b['times'])
            for key in ['T','H','V']:assert set(d[key+'_ids'])==set(b['ids'][m[key]].tolist())
        for record in d['extraction_metadata']['seed_diagnostics']:
            for key in ['raw_neighbor_physical_ids','candidate_physical_ids','fit_physical_ids','dropped_T_physical_ids']:assert set(record[key])<=T
        for row in d['holdout_H']['all_plane_H_pairs']:assert set(row['holdout_ids'])<=H
        for name,route in entry['routes'].items():
            assert route['LOO_status']=='not_evaluated_regression_global_only'
            assert route['qualification']['LOO_evaluable'] is None and route['qualification']['LOO_passed'] is None and route['qualification']['pass_rotation'] is None
            actual=d['routes'][name];assert actual['search'] is None and not actual['initializer']['available']
            assert not route['measurement_valid'] and route['measurement_omega'] is None and route['actual_NM_calls']==0
        splitchecks.append(dict(case=entry['case'],frontend=entry['frontend'],T=len(T),H=len(H),V=len(V),fit_T_only=True))
    for name in ['parameter_freeze.json','confirmation_regression_freeze.json']:
        for rel,h in json.loads((OUT/name).read_text())['files'].items():assert digest(OUT/rel)==h,rel
    for rel,h in json.loads((OUT/'input_source_hashes.json').read_text()).items():assert digest(ROOT/rel)==h,rel
    for rel,record in json.loads((PHASE4/'manifest_task13.json').read_text())['files'].items():assert digest(PHASE4/rel)==record['sha256'],rel
    report=dict(status='passed',regression_windows=108,frontend_method_states=648,shared_G_references=108,
        reused_cases=len(reusechecks),reused_frontend_cases=24,reuse_radius_and_source_sha_checks=reusechecks,
        available_input_frontend_cases_verified=len(splitchecks),geometry_missing_frontend_cases=60,T_H_V_disjoint_and_complete=True,
        all_fit_support_T_only=True,all_real_AS_initialization_failed=True,no_real_AS_search_or_G_fallback=True,
        regression_LOO_explicitly_not_evaluated_and_null=True,parameter_and_selection_hashes_unchanged=True,
        task13_manifest_unchanged=True,source_inputs_unchanged=len(sources),checks=splitchecks,
        G_native_valid=sum(bool(r['decision'] and r['decision']['measurement_valid']) for r in g),
        G_common_quality_pass=sum(bool(r['decision'] and r['decision']['common_quality']['passed']) for r in g),
        maximum_T_events=max(r['T'] for r in audit.values()),scope='Incremental evidence; prior synthetic/development/confirmation verifications cover frontend causality and prediction paths. Historical108 is not independent testing.')
    write_json(p/'verification.json',report);print({k:v for k,v in report.items() if k not in ['checks','reuse_radius_and_source_sha_checks']})

if __name__=='__main__':main()
