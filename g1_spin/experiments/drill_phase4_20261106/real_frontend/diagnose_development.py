"""Reconstruct saved T fit supports for seed diagnostics; no new estimates or thresholds."""
import gzip
import json
from collections import Counter
import numpy as np
from helpers import OUT,write_json,table
from run import bundle_real
from pipeline import split

def summarize(values):
    values=[v for v in values if v is not None]
    return dict(n=len(values),min=float(np.min(values)) if values else None,median=float(np.median(values)) if values else None,max=float(np.max(values)) if values else None)

def main():
    summaries=json.loads((OUT/'development/summary.json').read_text());manifest=json.loads((OUT/'input_manifest.json').read_text());rows=[];groups=[]
    for summary in summaries:
        if summary.get('status'):continue
        bundle=bundle_real(next(r for r in manifest if r['case']==summary['case']));T=split(bundle['xy'],bundle['times'])['T']
        detail=json.loads(gzip.decompress((OUT/summary['detail']).read_bytes()));velocity=np.array(detail['lane']['velocity']);raw=bundle['xy'][T];t=bundle['times'][T];xy=raw-(t-.005)[:,None]*velocity;ids=bundle['ids'][T];lookup={int(p):i for i,p in enumerate(ids)}
        local=[]
        for seed in detail['extraction_metadata']['seed_diagnostics']:
            indices=np.array([lookup[i] for i in seed['fit_physical_ids']],int);n=len(indices)
            rank=None;ratio=None;std=None;rms=None;gate=None
            if n>=2:
                coords=xy[indices];centered=coords-coords.mean(axis=0);sing=np.linalg.svd(centered,compute_uv=False)
                ratio=float(sing[1]/sing[0]) if len(sing)>1 and sing[0]>1e-12 else 0.
                X=np.column_stack([centered,np.ones(n)]);rank=int(np.linalg.matrix_rank(X));dt=t[indices]-t[seed['seed_index']];std=float(np.std(dt));gate=max(20e-6,.15*std)
                if rank==3:
                    coeff=np.linalg.lstsq(X,dt,rcond=None)[0];rms=float(np.sqrt(np.mean((X@coeff-dt)**2)))
            triplet=seed['triplet_diagnostic'] or {}
            row=dict(case=summary['case'],lane=summary['lane'],frontend=summary['frontend'],seed_physical_id=seed['seed_physical_id'],
                raw_neighbors=seed['raw_neighbors'],available_latest_neighbors=seed['frontend_neighbors'],cluster_size=n,
                hypotheses=triplet.get('hypotheses'),valid_hypotheses=triplet.get('valid_hypotheses'),selected_cluster_support=triplet.get('selected_support'),
                final_design_rank=rank,final_spatial_ratio=ratio,time_std_s=std,final_plane_RMS_s=rms,residual_gate_s=gate,
                residual_over_gate=None if rms is None else rms/gate,original_exit=seed['status'],
                rank_residual_scope='Independent audit of saved training support, including early-exit small supports; not a new acceptance decision')
            rows.append(row);local.append(row)
        groups.append(dict(case=summary['case'],lane=summary['lane'],frontend=summary['frontend'],T=summary['counts']['T'],H=summary['counts']['H'],V=summary['counts']['V'],
            seeds=len(local),exit_reasons=dict(Counter(r['original_exit'] for r in local)),
            raw_neighbors=summarize([r['raw_neighbors'] for r in local]),latest_neighbors=summarize([r['available_latest_neighbors'] for r in local]),
            latest_at_least10=sum(r['available_latest_neighbors']>=10 for r in local),clusters_at_least10=sum(r['cluster_size']>=10 for r in local),
            cluster_size=summarize([r['cluster_size'] for r in local]),valid_triplet_hypotheses=summarize([r['valid_hypotheses'] for r in local]),
            final_ranks=dict(Counter(str(r['final_design_rank']) for r in local)),final_spatial_ratio=summarize([r['final_spatial_ratio'] for r in local]),
            final_RMS_us=summarize([r['final_plane_RMS_s']*1e6 if r['final_plane_RMS_s'] is not None else None for r in local]),
            residual_over_gate=summarize([r['residual_over_gate'] for r in local])))
    table(OUT/'development/seed_diagnostics.csv',rows);write_json(OUT/'development/seed_diagnostics_summary.json',groups)
    print('Reconstructed',len(rows),'saved seed supports; no extraction/search or real parameter changes')

if __name__=='__main__':main()
