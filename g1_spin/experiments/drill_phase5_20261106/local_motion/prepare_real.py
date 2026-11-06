"""Fixed development patches split once; trainer receives only T NPZ files."""
import json
import hashlib
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent;P5=OUT.parent;P4=P5.parent/'drill_phase4_20261106'
DEVELOPMENT=['run1_01495000','run2_02095000','run5_04095000','run5_01525000']

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False))

def main():
    target=OUT/'real';target.mkdir(exist_ok=True)
    if (target/'development_inputs.json').exists():raise FileExistsError('Preserve fixed development split')
    for folder in ['T','evaluation_only']:(target/folder).mkdir(exist_ok=True)
    selections=json.loads((P5/'patch_selection.json').read_text());rows=[];sources={}
    for case in DEVELOPMENT:
        selection=next(r for r in selections if r['case']==case);patches=[]
        for p in selection['patches']:
            name=f'{case}_patch{p["patch"]}';source=P5/'patches'/f'{name}.npz';d=np.load(source)
            xy=d['xy'];times=d['time_s'];pol=d['polarity'];ids=d['ids'];hp=(xy[:,0]+2*xy[:,1])%5==0
            T=(times<.006)&~hp;H=(times<.006)&hp;V=times>=.006
            tpath=target/'T'/f'{name}.npz';epath=target/'evaluation_only'/f'{name}.npz'
            np.savez_compressed(tpath,xy=xy[T],time_s=times[T],polarity=pol[T],ids=ids[T])
            np.savez_compressed(epath,H_xy=xy[H],H_time_s=times[H],H_polarity=pol[H],H_ids=ids[H],
                V_xy=xy[V],V_time_s=times[V],V_polarity=pol[V],V_ids=ids[V])
            patches.append(dict(patch=p['patch'],box=p['box'],source_file=str(source.relative_to(P5)),source_sha256=digest(source),
                T_file=str(tpath.relative_to(OUT)),T_sha256=digest(tpath),T_events=int(T.sum()),
                evaluation_file=str(epath.relative_to(OUT)),evaluation_sha256=digest(epath),H_events=int(H.sum()),V_events=int(V.sum())))
            sources[str(source.relative_to(P5))]=digest(source)
        rows.append(dict(case=case,center=selection['center'],radius_px=selection['radius_px'],patches=patches,patch_count=len(patches),
            status='fixed_development_input' if patches else 'no_count_eligible_patch_activity_control',
            geometry_velocity_px_s=[0.,0.],geometry_source='stage4 frozen B1; no new geometry fitting'))
    write(target/'development_inputs.json',rows)
    config=dict(stage='Before any new real2x2 training output',date='2026-11-07',tau_px=.5,
        development_cases=DEVELOPMENT,confirmation_cases=['run3_01495000','run4_02495000','run5_05095000'],
        representations=['R_all','R_latest'],models=['M_time','M_spatial'],segments_s=[[0,.002],[.002,.004],[.004,.006]],
        polarities=[0,1],minimum_T_unique_pixels=12,bootstrap_draws=32,normal_spread95_max_deg=15,
        future_and_H_not_training_filters=True,real_geometry_lane='B1 fixed center, cdot0',
        bootstrap_sampling='Raw-pixel clusters, sample unique pixels with replacement, preserve all events in drawn clusters; deterministic SHA256 identity seed, same draws for both fits.',
        global_readiness='At least4 spatial patches,12 training-qualified local observations, rank3 condition<=100 and each leave-patch-out rank3 before global chain; no evaluation filtering.',
        sources=sources,files={'prepare_real.py':digest(OUT/'prepare_real.py'),'models.py':digest(OUT/'models.py'),'tau_freeze.json':digest(OUT/'tau_freeze.json'),
            'real/development_inputs.json':digest(target/'development_inputs.json')},
        qualification_status='Bootstrap and raw statistics saved first; stationary/nonzero/permutation/selection supplemental rules separate before evaluator.')
    write(OUT/'pre_real_freeze.json',config);print('Prepared T-only patches',sum(len(r['patches']) for r in rows),'cases',len(rows),'confirmation not opened')

if __name__=='__main__':main()
