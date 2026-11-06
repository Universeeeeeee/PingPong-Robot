"""Checkpoint verification: split IDs, frontend causality, holdouts and fixed synthetic controls."""
import argparse
import gzip
import json
import numpy as np
from helpers import OUT,PHASE4,ROOT,digest,write_json,flow
from run import bundle_synthetic,bundle_real
from pipeline import split
from prediction import in_box

def read_detail(path):return json.loads(gzip.decompress(path.read_bytes()))

def verify(kind):
    frozen=json.loads((OUT/'parameter_freeze.json').read_text())
    for rel,h in frozen['files'].items():assert digest(OUT/rel)==h,rel
    summary=json.loads((OUT/kind/'summary.json').read_text());checks=[];truth_errors=[];root_tests=json.loads((OUT/'path_verification.json').read_text());assert root_tests['status']=='passed'
    manifest=json.loads((OUT/'input_manifest.json').read_text())
    for entry in summary:
        if entry.get('status')=='geometry_lane_unavailable':checks.append(dict(case=entry['case'],frontend=entry['frontend'],lane=entry['lane'],status='preserved_geometry_failure'));continue
        if kind=='synthetic':bundle,_=bundle_synthetic(entry['case'])
        else:bundle=bundle_real(next(r for r in manifest if r['case']==entry['case']))
        path=OUT/entry['detail'];assert digest(path)==entry['detail_sha256'];result=read_detail(path)
        masks=split(bundle['xy'],bundle['times']);T,H,V=[bundle['ids'][masks[k]] for k in ['T','H','V']]
        assert np.array_equal(T,result['T_ids']) and np.array_equal(H,result['H_ids']) and np.array_equal(V,result['V_ids'])
        assert not np.intersect1d(T,H).size and not np.intersect1d(T,V).size and not np.intersect1d(H,V).size
        raw=bundle['xy'][masks['T']];times=bundle['times'][masks['T']];pol=bundle['polarity'][masks['T']];velocity=np.array(result['lane']['velocity']);xy=raw-(times-.005)[:,None]*velocity
        index={int(i):k for k,i in enumerate(T)}
        for record in result['extraction_metadata']['seed_diagnostics']:
            for key in ['raw_neighbor_physical_ids','candidate_physical_ids','fit_physical_ids','dropped_T_physical_ids']:
                assert set(record[key])<=set(T.tolist())
            seed=record['seed_index'];candidate=np.array([index[i] for i in record['candidate_physical_ids']],int)
            assert np.all(pol[candidate]==pol[seed]) and np.all(np.sum((xy[candidate]-xy[seed])**2,axis=1)<=16+1e-10)
            if result['frontend'] in ['F1','F2']:
                assert np.all(times[candidate]<=times[seed]) and np.all(times[candidate]>=times[seed]-.002)
                assert len(np.unique(raw[candidate],axis=0))==len(candidate)
                for i in candidate:
                    history=np.flatnonzero(np.all(raw==raw[i],axis=1)&(pol==pol[seed])&(times<=times[seed])&(times>=times[seed]-.002))
                    expected=max(history,key=lambda k:(times[k],int(T[k])));assert i==expected
        basepath=path.name.replace('.json.gz','_planes.npz');planes=np.load(path.parent/basepath)
        if entry['frontend']=='F0':
            old,_=flow.extract(xy,times,pol,bundle['center'],bundle['radius'])
            assert np.allclose(old['gradients'],planes['gradients'],atol=1e-12) and np.array_equal(old['seed_indices'],planes['seed_indices'])
        Hpred=result['holdout_H'];assert set(Hpred['predicted_H_ids'])<=set(H.tolist())
        for pair in Hpred['all_plane_H_pairs']:assert set(pair['holdout_ids'])<=set(H.tolist())
        for route in result['routes'].values():
            if route['search'] is None:assert not route['initializer']['available'] or len(T)<100
            else:
                assert route['search']['actual_NM_calls']==1 and not route['search']['G_fallback_called']
                for k,segment in enumerate(route['future_V']['segments']):
                    assert set(segment['source_ids'])<=set(T.tolist()) and set(segment['target_ids'])<=set(V.tolist())
                    assert segment['target_ids']==result['zero_future']['segments'][k]['target_ids']==result['translation_only_future']['segments'][k]['target_ids']
                for leave in route['leave_region_out']:
                    removed=set(leave['removed_T_ids']);remaining=set(leave['remaining_T_ids']);source=set(leave['evaluator_only_source_ids'])
                    assert not source&remaining and source<=removed and removed|remaining==set(T.tolist())
                    assert not source&set(H.tolist()) and not source&set(V.tolist())
                    if leave['accepted_and_predictively_supported']:assert leave['fitted']['decision']['measurement_valid']
            if route['decision'] is not None and not route['decision']['measurement_valid']:assert route['decision']['measurement_omega'] is None
        checks.append(dict(case=entry['case'],frontend=entry['frontend'],lane=entry['lane'],T=len(T),H=len(H),V=len(V),
            planes=entry['accepted_planes'],split_and_surface_causality_passed=True))
        if kind=='synthetic' and entry['lane']=='B1' and entry['case'] in ['base','reverse','second']:
            truth=np.load(PHASE4/'truth'/f'{entry["case"]}.npz')['omega'];route=result['routes']['S_comp'];init=np.array(route['initializer']['raw_omega_rad_s']);dec=route['decision']
            assert route['initializer']['available'] and dec and dec['measurement_valid']
            final=np.array(dec['measurement_omega']);angle=lambda a:float(np.degrees(np.arccos(np.clip(a@truth/np.linalg.norm(a)/np.linalg.norm(truth),-1,1))))
            error=dict(case=entry['case'],frontend=entry['frontend'],initial_axis_error_deg=angle(init),final_axis_error_deg=angle(final),
                initial_speed_error=float(np.linalg.norm(init)/np.linalg.norm(truth)-1),final_speed_error=float(np.linalg.norm(final)/np.linalg.norm(truth)-1),
                V_pass=route['comparison_V']['engineering_pass'],translation_comparison_pass=route['model_vs_translation_only']['engineering_pass'])
            assert error['initial_axis_error_deg']<=10 and abs(error['initial_speed_error'])<=.15 and error['final_axis_error_deg']<=5 and abs(error['final_speed_error'])<=.1
            assert error['V_pass'] and error['translation_comparison_pass'];truth_errors.append(error)
    for rel,h in json.loads((OUT/'input_source_hashes.json').read_text()).items():assert digest(ROOT/rel)==h,rel
    for rel,record in json.loads((PHASE4/'manifest_task13.json').read_text())['files'].items():assert digest(PHASE4/rel)==record['sha256'],rel
    report=dict(status='passed',stage=kind,states=len(summary),checks=checks,synthetic_S_errors=truth_errors,
        F0_matches_original_extractor=True,T_H_V_ID_isolation_passed=True,global_latest_before_lane_radius_passed=True,
        all_frontend_fit_ids_are_T=True,all_future_targets_are_V=True,zero_translation_model_same_target_ids=True,
        task13_manifest_unchanged=True,input_sources_unchanged=True,measurement_path_controls=root_tests,
        caveat='Synthetic truth errors verify clean-edge code paths only; H/V from same recording are diagnostic, not independent repetitions.')
    write_json(OUT/kind/'verification.json',report);print('PASSED',kind,len(checks),'states;',len(truth_errors),'synthetic full-chain cases')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['synthetic','development']);args=parser.parse_args();verify(args.stage)
