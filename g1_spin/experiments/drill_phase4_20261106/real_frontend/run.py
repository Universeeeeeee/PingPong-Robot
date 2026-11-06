"""Bounded candidate checks and development only; confirmation requires next checkpoint."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from helpers import OUT,PHASE4,methods,digest,write_json,table,hash_arrays
from pipeline import process_case,split

def save_detail(path,result):
    content=json.dumps(result,allow_nan=False,separators=(',',':')).encode()
    path.write_bytes(gzip.compress(content,mtime=0))

def bundle_synthetic(name):
    d=np.load(PHASE4/'inputs'/f'{name}.npz');g=json.loads((PHASE4/'inputs'/f'{name}_geometry.json').read_text())
    return dict(case=name,source='task1_fixed_integer_edge_crossings',xy=np.column_stack([d['x'],d['y']]).astype(float),times=(d['t_us']-10000)*1e-6,
        polarity=d['polarity'],ids=np.arange(len(d['x'])),center=np.array(g['center_px']),radius=float(g['radius_px'])),[
        dict(name='B1',available=True,velocity=[0.,0.],source='fixed center'),dict(name='secant',available=True,velocity=g['center_velocity_px_s'],source='given synthetic trajectory')]

def bundle_real(row):
    d=np.load(OUT/row['input_file']);return dict(case=row['case'],source=row['source'],xy=d['xy'],times=d['times'],polarity=d['polarity'],ids=d['ids'],center=d['center'],radius=float(d['radius']))

def compact(result,path):
    routes={}
    for name,r in result['routes'].items():
        init=r['initializer'];dec=r['decision'];search=r['search']
        routes[name]=dict(initializer=init,measurement_valid=bool(dec and dec['measurement_valid']),failure_reason=r['status'],
            candidate_omega=None if search is None else search['optimizer_candidate_omega'],measurement_omega=None if dec is None else dec['measurement_omega'],
            common_quality=None if dec is None else dec['common_quality'],native_scalar_status=init['native_scalar_status'],
            future_comparison=r['comparison_V'],model_vs_translation_only=r['model_vs_translation_only'],qualification=r['qualification'],pipeline_s=r['pipeline_s'],flow_s=r['flow_s'],
            actual_NM_calls=0 if search is None else search['actual_NM_calls'])
    return dict(case=result['case'],frontend=result['frontend'],lane=result['lane']['name'],counts=result['counts'],input_hash=result['coordinate_input_hash'],
        accepted_planes=result['extraction_metadata']['accepted_planes'],seed_exit_reasons=result['extraction_metadata']['failures'],
        normal_support=result['normal_support'],H_statistics=result['holdout_H']['event_statistics'],H_pass=result['holdout_H']['engineering_pass'],
        H_total=result['holdout_H']['H_total'],H_predicted=result['holdout_H']['unique_H_events'],H_total_prediction_coverage=result['holdout_H']['total_H_prediction_coverage'],
        H_regions=result['holdout_H']['regions_with_at_least5_H'],translation_vs_zero=result['translation_vs_zero'],
        S_raw_RHS={k:v for k,v in result['S_raw_RHS'].items() if k!='metadata'},routes=routes,detail=str(path.relative_to(OUT)),detail_sha256=digest(path))

def execute(kind):
    directory=OUT/kind;directory.mkdir(exist_ok=True)
    if (directory/'summary.json').exists():raise FileExistsError('Preserve completed bounded pass')
    frozen=json.loads((OUT/'parameter_freeze.json').read_text())
    for rel,h in frozen['files'].items():assert digest(OUT/rel)==h,rel
    jobs=[]
    if kind=='synthetic':
        for scene in ['base','reverse','second','zero','translation','mixed']:
            b,lanes=bundle_synthetic(scene);jobs.append((b,lanes))
    elif kind=='development':
        assert (OUT/'synthetic/verification.json').exists(),'Synthetic branch checks required first'
        for row in json.loads((OUT/'input_manifest.json').read_text()):
            if row['partition']=='development':jobs.append((bundle_real(row),row['lanes']))
    else:raise ValueError('Confirmation/regression not in this checkpoint')
    summaries=[];Grefs=[]
    for bundle,lanes in jobs:
        masks=split(bundle['xy'],bundle['times']);T=masks['T']
        for lane in lanes:
            if not lane['available']:
                for frontend in ['F0','F1','F2']:summaries.append(dict(case=bundle['case'],frontend=frontend,lane=lane['name'],status='geometry_lane_unavailable',reason=lane['failure_reason']))
                continue
            velocity=np.array(lane['velocity']);xy=bundle['xy'][T]-bundle['center']-(bundle['times'][T]-.005)[:,None]*velocity
            points=methods.cmax.lift(xy[:,0],xy[:,1],0,0,bundle['radius']);search=methods.g_search(points,bundle['times'][T]);decision=methods.decision(points,bundle['times'][T],search,'G')
            Grefs.append(dict(case=bundle['case'],lane=lane['name'],T_ids=bundle['ids'][T].tolist(),input_hash=hash_arrays(bundle['ids'][T],points,bundle['times'][T]),search=search,decision=decision,
                note='Same T/lane reference only; never fallback and not frontend selection. Frozen G may internally cap5000 events.'))
            for frontend in ['F0','F1','F2']:
                result,planes,raw=process_case(bundle,lane,frontend,include_LOO=kind=='development')
                name=f'{bundle["case"]}_{lane["name"]}_{frontend}';path=directory/f'{name}.json.gz';save_detail(path,result)
                np.savez_compressed(directory/f'{name}_planes.npz',**planes);np.savez_compressed(directory/f'{name}_raw_planes.npz',**raw)
                entry=compact(result,path);summaries.append(entry)
                print(kind,name,'planes',entry['accepted_planes'],'H',entry['H_statistics']['median'],
                    'S',entry['routes']['S_comp']['failure_reason'],'V',entry['routes']['S_comp']['future_comparison'],flush=True)
    write_json(directory/'summary.json',summaries);save_detail(directory/'G_reference.json.gz',Grefs)
    table(directory/'summary.csv',summaries)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['synthetic','development']);args=parser.parse_args();execute(args.stage)
