"""Correct clipped confirmation pixels to full decoded raw patches; fixed choices/parameters."""
import json
from pathlib import Path
import numpy as np
from synthetic import OUT,P4,digest,write

ROOT=OUT.parents[3];ARCHIVE=ROOT/'g1_spin/data/archive_20261105/decoded'
NAMES={3:'DSP03F100025_20261105_201143',4:'DSP03F100025_20261105_203829',5:'DSP03F100025_20261105_215030'}

def main():
    real=OUT/'real';prior=OUT/'prereview_clipped';old=json.loads((prior/'confirmation_inputs.json').read_text());original=json.loads((P4/'real_frontend/input_manifest.json').read_text());corrected=[];audit=[];sources={}
    for folder in ['confirmation_T','confirmation_evaluation_only']:(real/folder).mkdir(exist_ok=True)
    for case in old:
        run=int(case['case'][3]);record=next(r for r in original if r['case']==case['case']);folder=ARCHIVE/NAMES[run];m=json.loads((folder/'manifest.json').read_text())
        events=np.concatenate([np.load(folder/name) for name in m['parts']]);origin=int(events['t'][0]);start=origin+record['start_us'];patches=[]
        sources[str(folder/'manifest.json')]=digest(folder/'manifest.json')
        for name in m['parts']:sources[str(folder/name)]=digest(folder/name)
        for patch in case['patches']:
            box=patch['box'];keep=(events['t']>=start)&(events['t']<start+10000)&(events['x']>=box[0])&(events['x']<box[2])&(events['y']>=box[1])&(events['y']<box[3])
            ids=np.flatnonzero(keep);raw=events[ids];xy=np.column_stack([raw['x'],raw['y']]).astype(float);times=(raw['t'].astype(np.int64)-start)*1e-6;pol=raw['p'];hp=(xy[:,0].astype(int)+2*xy[:,1].astype(int))%5==0
            T=(times<.006)&~hp;H=(times<.006)&hp;V=times>=.006;name=f'{case["case"]}_patch{patch["patch"]}'
            oldT=np.load(prior/'confirmation_T'/f'{name}.npz');oldEval=np.load(prior/'confirmation_evaluation_only'/f'{name}.npz');oldids=np.sort(np.r_[oldT['ids'],oldEval['H_ids'],oldEval['V_ids']]);assert np.all(np.isin(oldids,ids))
            tpath=real/'confirmation_T'/f'{name}.npz';epath=real/'confirmation_evaluation_only'/f'{name}.npz'
            reused=case['case']=='run4_02495000'
            if reused:
                assert np.array_equal(oldT['ids'],ids[T]) and np.array_equal(oldT['xy'],xy[T]) and np.array_equal(oldT['time_s'],times[T]) and np.array_equal(oldT['polarity'],pol[T])
                for key,mask in [('H',H),('V',V)]:
                    assert np.array_equal(oldEval[key+'_ids'],ids[mask]) and np.array_equal(oldEval[key+'_xy'],xy[mask]) and np.array_equal(oldEval[key+'_time_s'],times[mask]) and np.array_equal(oldEval[key+'_polarity'],pol[mask])
                import shutil
                shutil.copy2(prior/'confirmation_T'/f'{name}.npz',tpath);shutil.copy2(prior/'confirmation_evaluation_only'/f'{name}.npz',epath)
            else:
                np.savez_compressed(tpath,xy=xy[T],time_s=times[T],polarity=pol[T],ids=ids[T])
                np.savez_compressed(epath,H_xy=xy[H],H_time_s=times[H],H_polarity=pol[H],H_ids=ids[H],V_xy=xy[V],V_time_s=times[V],V_polarity=pol[V],V_ids=ids[V])
            audit.append(dict(case=case['case'],patch=patch['patch'],fixed_box=box,old_events=len(oldids),corrected_events=len(ids),added_event_ids=np.setdiff1d(ids,oldids).tolist(),
                old_T_events=len(oldT['ids']),new_T_events=int(T.sum()),old_H_events=len(oldEval['H_ids']),new_H_events=int(H.sum()),old_V_events=len(oldEval['V_ids']),new_V_events=int(V.sum()),
                raw_boolean_ID_reconstruction=True,reused_run4_unchanged=reused,old_T_sha256=digest(prior/'confirmation_T'/f'{name}.npz'),new_T_sha256=digest(tpath),
                old_evaluation_sha256=digest(prior/'confirmation_evaluation_only'/f'{name}.npz'),new_evaluation_sha256=digest(epath)))
            patches.append(dict(**{k:v for k,v in patch.items() if k not in ['T_sha256','evaluation_sha256','T_events','H_events','V_events']},T_sha256=digest(tpath),evaluation_sha256=digest(epath),T_events=int(T.sum()),H_events=int(H.sum()),V_events=int(V.sum())))
        corrected.append(dict(**{k:v for k,v in case.items() if k!='patches'},patches=patches,patch_event_source='Full official decoded integer-pixel10ms raw patch; old common support used only to choose fixed box'))
        del events
    expected={('run3_01495000',2):249,('run3_01495000',3):259,('run3_01495000',4):233,('run5_05095000',0):100,('run5_05095000',1):117,('run5_05095000',3):79}
    for (case,patch),n in expected.items():assert next(r for r in audit if r['case']==case and r['patch']==patch)['corrected_events']==n
    write(real/'confirmation_inputs.json',corrected);write(OUT/'confirmation_input_correction.json',dict(status='passed',reason='Previous confirmation used clipped stage4 common_ids to fill patches, unlike dev full decoded raw events; fixed raw observation bug only.',
        old_results='prereview_clipped/',fixed_patch_boxes_preserved=True,geometry_tau_bootstrap_permutation_parameters_unchanged=True,run4_reused_after_exact_ID_content_and_hash_check=True,
        patch_audit=audit,raw_source_hashes=sources,required_reproduction_counts=expected.__len__(),selection_none_sha256=digest(real/'selection_development.json')))
    print('RAW INPUT REPAIR PASSED',[(r['case'],r['patch'],r['old_events'],r['corrected_events']) for r in audit],flush=True)

if __name__=='__main__':main()
