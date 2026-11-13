"""Fixed3-window confirmation with unchanged development parameters and no selection."""
import json
from pathlib import Path
import numpy as np
from synthetic import OUT,P4,digest,write,select_patches
from models import representation,fit
from train_real import bootstrap
from qualify_training import classification,build_system
from evaluate_real import evaluate

CASES=['run3_01495000','run4_02495000','run5_05095000']

def main():
    real=OUT/'real'
    if (real/'confirmation_training.json').exists():raise FileExistsError('Preserve confirmation training')
    selection_hash=digest(real/'selection_development.json');assert json.loads((real/'selection_development.json').read_text())['selected_combination'] is None
    assert json.loads((OUT/'geometry_path_amendment.json').read_text())['status']=='passed'
    source=json.loads((P4/'real_frontend/input_manifest.json').read_text());manifest=[];rows=[]
    (real/'confirmation_T').mkdir(exist_ok=True);(real/'confirmation_evaluation_only').mkdir(exist_ok=True)
    for case in CASES:
        src=next(r for r in source if r['case']==case);sourcefile=P4/'real_frontend'/src['input_file'];assert digest(sourcefile)==src['input_sha256'];d=np.load(sourcefile)
        xy=d['xy'];times=d['times'];pol=d['polarity'];ids=d['ids']
        # Same count-only grid rule as raw-structure stage: selection does not use H or residual values.
        cells=np.floor(xy[times<.006]/12).astype(int);unique,count=np.unique(cells,axis=0,return_counts=True);chosen=[]
        for i in sorted(range(len(unique)),key=lambda i:(-count[i],int(unique[i,0]),int(unique[i,1]))):
            if count[i]<20:continue
            center=unique[i]*12+6
            if any(np.linalg.norm(center-np.array(p['grid_center']))<16 for p in chosen):continue
            box=np.r_[unique[i]*12,unique[i]*12+12].tolist();chosen.append(dict(patch=len(chosen),grid_center=center.tolist(),box=box,selection_first6ms_events=int(count[i])))
            if len(chosen)==6:break
        patches=[]
        for p in chosen:
            box=p['box'];mask=(xy[:,0]>=box[0])&(xy[:,0]<box[2])&(xy[:,1]>=box[1])&(xy[:,1]<box[3]);px=xy[mask];t=times[mask];pl=pol[mask];physical=ids[mask]
            Hpix=(px[:,0].astype(int)+2*px[:,1].astype(int))%5==0;T=(t<.006)&~Hpix;H=(t<.006)&Hpix;V=t>=.006
            name=f'{case}_patch{p["patch"]}';tp=real/'confirmation_T'/f'{name}.npz';ep=real/'confirmation_evaluation_only'/f'{name}.npz'
            np.savez_compressed(tp,xy=px[T],time_s=t[T],polarity=pl[T],ids=physical[T])
            np.savez_compressed(ep,H_xy=px[H],H_time_s=t[H],H_polarity=pl[H],H_ids=physical[H],V_xy=px[V],V_time_s=t[V],V_polarity=pl[V],V_ids=physical[V])
            patch=dict(**p,T_file=str(tp.relative_to(OUT)),T_sha256=digest(tp),evaluation_file=str(ep.relative_to(OUT)),evaluation_sha256=digest(ep),T_events=int(T.sum()),H_events=int(H.sum()),V_events=int(V.sum()));patches.append(patch)
        manifest.append(dict(case=case,center=src['center'],radius_px=src['radius_px'],patches=patches,patch_count=len(patches),source_file=str(sourcefile),source_sha256=digest(sourcefile),geometry_velocity_px_s=[0.,0.]))
    write(real/'confirmation_inputs.json',manifest)
    # Trainer opens only the separately saved T file, not H/V arrays.
    for case in manifest:
        for p in case['patches']:
            d=np.load(OUT/p['T_file']);xy=d['xy'];times=d['time_s'];pol=d['polarity'];ids=d['ids']
            for segment in range(3):
                low=segment*.002;high=low+.002
                for polarity in [0,1]:
                    tm=(times>=low-1e-12)&(times<high-1e-12)&(pol==polarity);local=np.flatnonzero(tm)
                    for rep in ['R_all','R_latest']:
                        picked=local[representation(xy[local],times[local],ids[local],rep)];identity=f'{case["case"]}/patch{p["patch"]}/segment{segment}/pol{polarity}/{rep}'
                        for kind in ['M_time','M_spatial']:
                            model=fit(xy[picked],times[picked],kind);boot=bootstrap(xy[picked],times[picked],kind,identity,model)
                            r=dict(observation_id=len(rows),case=case['case'],patch=p['patch'],box=p['box'],segment=segment,polarity=polarity,representation=rep,model_kind=kind,
                                start_s=low,end_s=high,T_file=p['T_file'],T_sha256=p['T_sha256'],all_segment_T_ids=ids[local].tolist(),support_ids=ids[picked].tolist(),model=model,bootstrap=boot)
                            r['T_qualification']=classification(r);rows.append(r)
        print('confirmation T trained',case['case'],'patches',len(case['patches']),flush=True)
    write(real/'confirmation_training.json',rows);write(real/'confirmation_training_qualification.json',rows)
    systems=[]
    for case in manifest:
        for rep in ['R_all','R_latest']:
            for kind in ['M_time','M_spatial']:
                local=[r for r in rows if r['case']==case['case'] and r['representation']==rep and r['model_kind']==kind]
                systems.append(dict(case=case['case'],representation=rep,model_kind=kind,**build_system(local,np.array(case['center']),case['radius_px'])))
    write(real/'confirmation_T_global_systems.json',systems)
    write(real/'confirmation_training_freeze.json',dict(stage='All confirmation T fits/bootstrap/qualification/global systems before evaluation',rows=len(rows),
        files={r:digest(OUT/r) for r in ['confirm_real.py','models.py','train_real.py','qualify_training.py','real/confirmation_inputs.json','real/confirmation_training.json','real/confirmation_training_qualification.json','real/confirmation_T_global_systems.json']},
        selection_development_sha256=selection_hash,H_future_evaluation_performed=False))
    print('FROZEN confirmation T',len(rows),'global sufficient',sum(s['global_training_sufficient'] for s in systems),flush=True)
    evaluate('confirmation');assert digest(real/'selection_development.json')==selection_hash

if __name__=='__main__':main()
