"""Single historical108 regression, all3F, global H/V only; reuse8 B1 cases."""
import gzip
import json
import numpy as np
from helpers import OUT,ROOT,P3,methods,digest,write_json,table,hash_arrays,load_archive
from pipeline import process_case,split
from run import save_detail,compact

def main():
    assert json.loads((OUT/'selection_development.json').read_text())['selected_frontend'] is None
    assert json.loads((OUT/'confirmation/verification.json').read_text())['status']=='passed'
    directory=OUT/'regression108';directory.mkdir(exist_ok=True)
    if (directory/'summary.json').exists():raise FileExistsError('Preserve single108 regression')
    audit=json.loads((P3/'window_audit.json').read_text());known=json.loads((OUT/'input_manifest.json').read_text());known={r['case']:r for r in known}
    previous={};previousG={}
    for partition in ['development','confirmation']:
        for r in json.loads((OUT/partition/'summary.json').read_text()):
            if r['lane']=='B1':previous[(r['case'],r['frontend'])]=r
        for r in json.loads(gzip.decompress((OUT/partition/'G_reference.json.gz').read_bytes())):
            if r['lane']=='B1':previousG[r['case']]=r
    summary=[];Grefs=[];reuse=[];inputs=[]
    for run in [1,2,3,4,5]:
        events,origin=load_archive(run)
        for w in [w for w in audit if w['run']==run]:
            case=w['window'];lane=dict(name='B1',available=True,velocity=[0.,0.],source='historical stage3 fixed B1')
            if w['geometry_status']=='unavailable':
                Grefs.append(dict(case=case,lane='B1',status='geometry_unavailable',search=None,decision=None))
                for f in ['F0','F1','F2']:
                    summary.append(dict(case=case,frontend=f,lane='B1',status='geometry_unavailable',source_stage3_condition=w['condition'],
                        routes={m:dict(measurement_valid=False,measurement_omega=None,candidate_omega=None,failure_reason='geometry_unavailable',
                            LOO_status='not_evaluated_regression_global_only') for m in ['A_literal','S_comp']}))
                continue
            source=P3/w['input_file'];old=np.load(source);ids=old['common_ids'];raw=events[ids];xy=np.column_stack([raw['x'],raw['y']]).astype(float);times=(raw['t'].astype(np.int64)-origin-w['start_us'])*1e-6
            bundle=dict(case=case,source='historical108_complete_common_ids',xy=xy,times=times,polarity=raw['p'],ids=ids,center=old['center'],radius=float(w['radius_px']))
            input_hash=hash_arrays(ids,xy,times,raw['p'],old['center']);T=split(xy,times)['T']
            inputs.append(dict(case=case,run=run,condition=w['condition'],source_input=str(source.relative_to(ROOT)),source_input_sha256=digest(source),
                common_ids_hash=hash_arrays(ids),input_hash=input_hash,complete_events=len(ids),T=int(T.sum()),source_mask=w['source_mask']))
            if case in known:
                r=known[case];assert r['input_hash']==input_hash and digest(OUT/r['input_file'])==r['input_sha256']
                Grefs.append({**previousG[case], 'reused':True})
                for f in ['F0','F1','F2']:
                    entry=previous[(case,f)];assert digest(OUT/entry['detail'])==entry['detail_sha256']
                    copy=json.loads(json.dumps(entry));copy.update(reused=True,reused_from=entry['detail'],source_stage3_condition=w['condition'],LOO_status='not_evaluated_regression_global_only')
                    for route in copy['routes'].values():
                        route['qualification']['LOO_evaluable']=None;route['qualification']['LOO_passed']=None;route['qualification']['pass_rotation']=None
                        route['LOO_status']='not_evaluated_regression_global_only'
                    summary.append(copy);reuse.append(dict(case=case,frontend=f,input_hash=input_hash,detail_sha256=entry['detail_sha256']))
                print('regression reused',case,flush=True);continue
            relative=xy[T]-old['center'];points=methods.cmax.lift(relative[:,0],relative[:,1],0,0,w['radius_px']);search=methods.g_search(points,times[T]);decision=methods.decision(points,times[T],search,'G')
            Grefs.append(dict(case=case,lane='B1',reused=False,T_ids=ids[T].tolist(),input_hash=hash_arrays(ids[T],points,times[T]),search=search,decision=decision))
            for f in ['F0','F1','F2']:
                result,planes,rawplanes=process_case(bundle,lane,f,include_LOO=False);name=f'{case}_B1_{f}'
                result['LOO_status']='not_evaluated_regression_global_only'
                for route in result['routes'].values():
                    route['LOO_status']='not_evaluated_regression_global_only';route['qualification']['LOO_evaluable']=None;route['qualification']['LOO_passed']=None;route['qualification']['pass_rotation']=None
                path=directory/f'{name}.json.gz';save_detail(path,result);np.savez_compressed(directory/f'{name}_planes.npz',**planes);np.savez_compressed(directory/f'{name}_raw_planes.npz',**rawplanes)
                entry=compact(result,path);entry.update(reused=False,source_stage3_condition=w['condition'],LOO_status='not_evaluated_regression_global_only')
                for route in entry['routes'].values():route['LOO_status']='not_evaluated_regression_global_only'
                summary.append(entry)
            print('regression',case,'T',int(T.sum()),'planes',[r['accepted_planes'] for r in summary[-3:]],flush=True)
        del events
    assert len(summary)==324 and len(Grefs)==108 and len(reuse)==24
    write_json(directory/'summary.json',summary);table(directory/'summary.csv',summary);save_detail(directory/'G_reference.json.gz',Grefs)
    write_json(directory/'input_audit.json',inputs);write_json(directory/'reused_eight_windows.json',reuse)
    states=[]
    for entry in summary:
        for method,route in entry['routes'].items():states.append(dict(case=entry['case'],frontend=entry['frontend'],method=method,
            condition=entry['source_stage3_condition'],reused=entry.get('reused',False),**route))
    assert len(states)==648;write_json(directory/'states648.json',states);table(directory/'states648.csv',states)

if __name__=='__main__':main()
