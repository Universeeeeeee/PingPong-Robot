"""Bounded27 support-only comparisons at original T physical seeds; no searches."""
import gzip
import json
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree
from helpers import OUT,flow,write_json,table,digest
from frontends import stable_surface,cluster
from run import bundle_real
from pipeline import split

def statistics(a):
    a=[x for x in a if x is not None];return dict(n=len(a),min=float(np.min(a)) if a else None,median=float(np.median(a)) if a else None,max=float(np.max(a)) if a else None)

def extract_at_seeds(bundle,mask,physical_seeds,frontend):
    xy=bundle['xy'][mask];times=bundle['times'][mask];pol=bundle['polarity'][mask];ids=bundle['ids'][mask];lookup={int(v):i for i,v in enumerate(ids)};tree=cKDTree(xy)
    rows=[];points=[];gradients=[]
    for physical in physical_seeds:
        seed=lookup[physical];indices=np.asarray(tree.query_ball_point(xy[seed],4),int)
        indices=indices[(abs(times[indices]-times[seed])<=.002)&(pol[indices]==pol[seed])];original=indices.copy()
        if frontend in ['F1','F2']:
            history=np.flatnonzero((pol==pol[seed])&(times<=times[seed])&(times>=times[seed]-.002));latest=stable_surface(history,xy.astype(int),times,ids,seed)
            indices=latest[np.sum((xy[latest]-xy[seed])**2,axis=1)<=16]
        candidate=indices.copy();diagnostic=None
        if frontend=='F2':indices,diagnostic=cluster(indices,xy,times,seed)
        n=len(indices);rank=None;ratio=None;std=None;rms=None;gate=None;gradient=None;reason='too_few_neighbors'
        if n>=2:
            local=xy[indices];centered=local-local.mean(axis=0);s=np.linalg.svd(centered,compute_uv=False);ratio=float(s[1]/s[0]) if len(s)>1 and s[0]>1e-6 else 0.
            dt=times[indices]-times[seed];std=float(np.std(dt));gate=max(20e-6,.15*std);X=np.column_stack([centered,np.ones(n)]);coef,_,rank,_=np.linalg.lstsq(X,dt,rcond=None);rank=int(rank)
            rms=float(np.sqrt(np.mean((X@coef-dt)**2)));gradient=float(np.linalg.norm(coef[:2]))
            if n>=10:
                sx,sy=(local.mean(axis=0)-bundle['center'])/bundle['radius']
                if ratio<.08:reason='spatial_plane_degenerate'
                elif std<1e-6:reason='insufficient_local_time_span'
                elif rank<3 or rms>gate:reason='plane_residual_or_rank'
                elif sx*sx+sy*sy>=.9**2:reason='outside_flow_inner_sphere'
                elif gradient<1e-8:reason='vanishing_gradient'
                else:reason='accepted';points.append([sx,sy,-np.sqrt(1-sx*sx-sy*sy)]);gradients.append(coef[:2])
        rows.append(dict(seed_physical_id=physical,source_event_count=len(ids),source_unique_pixels=len(np.unique(xy,axis=0)),raw_neighbors=len(original),
            raw_unique_pixels=len(np.unique(xy[original],axis=0)),latest_neighbors=len(candidate),cluster_size=n,
            triplet_hypotheses=None if diagnostic is None else diagnostic.get('hypotheses'),valid_triplets=None if diagnostic is None else diagnostic.get('valid_hypotheses'),
            fit_physical_ids=ids[indices].tolist(),design_rank=rank,spatial_ratio=ratio,time_std_s=std,RMS_s=rms,gate_s=gate,RMS_over_gate=None if rms is None else rms/gate,status=reason))
    init={method:flow.initialize(np.asarray(points).reshape(-1,3),np.asarray(gradients).reshape(-1,2),bundle['radius'],method) for method in ['A','S']}
    return rows,init

def main():
    path=OUT/'split_diagnostic';path.mkdir(exist_ok=True)
    if (path/'summary.json').exists():raise FileExistsError('Preserve bounded27 comparisons')
    manifest=json.loads((OUT/'input_manifest.json').read_text());summaries=[];allrows=[]
    for name in ['run1_01495000','run2_02095000','run5_04095000']:
        bundle=bundle_real(next(r for r in manifest if r['case']==name));masks=split(bundle['xy'],bundle['times'])
        detail=json.loads(gzip.decompress((OUT/'development'/f'{name}_B1_F0.json.gz').read_bytes()))
        seeds=[r['seed_physical_id'] for r in detail['extraction_metadata']['seed_diagnostics']]
        for support,mask in [('T',masks['T']),('T_union_H',masks['T']|masks['H']),('full10ms',np.ones(len(bundle['ids']),bool))]:
            for frontend in ['F0','F1','F2']:
                rows,init=extract_at_seeds(bundle,mask,seeds,frontend)
                for r in rows:r.update(case=name,support=support,frontend=frontend)
                if support=='T':
                    saved=json.loads(gzip.decompress((OUT/'development'/f'{name}_B1_{frontend}.json.gz').read_bytes()))
                    assert [r['status'] for r in rows]==[r['status'] for r in saved['extraction_metadata']['seed_diagnostics']]
                summaries.append(dict(case=name,support=support,frontend=frontend,seeds=len(seeds),physical_seed_ids=seeds,
                    source_events=int(mask.sum()),source_unique_pixels=len(np.unique(bundle['xy'][mask],axis=0)),
                    accepted_planes=sum(r['status']=='accepted' for r in rows),exit_reasons=dict(Counter(r['status'] for r in rows)),
                    raw_neighbors=statistics([r['raw_neighbors'] for r in rows]),latest_neighbors=statistics([r['latest_neighbors'] for r in rows]),
                    cluster_size=statistics([r['cluster_size'] for r in rows]),RMS_over_gate=statistics([r['RMS_over_gate'] for r in rows]),initializers=init,
                    role='Measurement support diagnostic only; full10ms uses H/V and is not independent evaluation or a replacement frontend.'))
                allrows+=rows;print(name,support,frontend,'accepted',summaries[-1]['accepted_planes'],flush=True)
    assert len(summaries)==27
    table(path/'seeds.csv',allrows);write_json(path/'summary.json',summaries)
    write_json(path/'verification.json',dict(status='passed',combinations=27,same_original_T_physical_seeds=True,T_exit_reasons_match_main_results=True,
        no_search=True,no_frontend_selection_change=True,rows=len(allrows),selection_file_sha256=digest(OUT/'selection_development.json')))

if __name__=='__main__':main()
