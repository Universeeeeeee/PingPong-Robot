"""Three bounded, fixed frontends. Inputs contain T events only."""
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree

def stable_surface(indices,raw_pixels,times,physical_ids,seed):
    indices=indices[(times[indices]<=times[seed])&(times[indices]>=times[seed]-.002)]
    if not len(indices):return indices
    order=np.lexsort((physical_ids[indices],times[indices],raw_pixels[indices,1],raw_pixels[indices,0]))
    ordered=indices[order];pixels=raw_pixels[ordered];last=np.r_[np.any(pixels[1:]!=pixels[:-1],axis=1),True]
    return np.sort(ordered[last])

def cluster(indices,xy,times,seed):
    others=indices[indices!=seed]
    if len(others)<2:return np.empty(0,int),dict(hypotheses=0,reason='no_triplets')
    a,b=np.triu_indices(len(others),1);pairs=np.column_stack([others[a],others[b]])
    if len(pairs)>64:pairs=pairs[np.linspace(0,len(pairs)-1,64).astype(int)]
    best=None;valid_hypotheses=0
    for ordinal,(i,j) in enumerate(pairs):
        basis=xy[[i,j]]-xy[seed];det=float(np.linalg.det(basis))
        if abs(det)/2<2:continue
        gradient=np.linalg.solve(basis,times[[i,j]]-times[seed]);norm=float(np.linalg.norm(gradient))
        if norm<1e-8:continue
        distances=abs(times[indices]-times[seed]-(xy[indices]-xy[seed])@gradient)/norm
        keep=distances<=.75;support=int(keep.sum());median=float(np.median(distances));valid_hypotheses+=1
        key=(support,-median,-ordinal)
        if best is None or key>best[0]:best=(key,indices[keep],dict(hypotheses=len(pairs),valid_hypotheses=valid_hypotheses,selected_triplet=[int(seed),int(i),int(j)],
            selected_support=support,median_distance_all_neighbors_px=median,ordinal=ordinal))
    if best is None:return np.empty(0,int),dict(hypotheses=len(pairs),valid_hypotheses=0,reason='no_valid_triplet')
    best[2]['valid_hypotheses']=valid_hypotheses
    return best[1],best[2]

def extract(xy,times,polarity,raw_pixels,physical_ids,center,radius,frontend):
    xy=np.asarray(xy,float);times=np.asarray(times,float);raw_pixels=np.asarray(raw_pixels);physical_ids=np.asarray(physical_ids)
    failures=Counter();records=[];points=[];gradients=[];centroids=[];seed_indices=[];coeff=[];rmses=[];stds=[];neighbor_ids=[];offsets=[0];dropped_ids=[];dropped_offsets=[0]
    tree=cKDTree(xy);seeds=np.linspace(0,len(xy)-1,min(300,len(xy)),dtype=int) if len(xy)>=100 else np.empty(0,int)
    for seed in seeds:
        indices=np.asarray(tree.query_ball_point(xy[seed],4),int)
        indices=indices[(abs(times[indices]-times[seed])<=.002)&(polarity[indices]==polarity[seed])]
        original=indices.copy();extra=None
        if frontend in ['F1','F2']:
            history=np.flatnonzero((polarity==polarity[seed])&(times<=times[seed])&(times>=times[seed]-.002))
            latest=stable_surface(history,raw_pixels,times,physical_ids,seed)
            indices=latest[np.sum((xy[latest]-xy[seed])**2,axis=1)<=16]
        candidate=indices.copy()
        if frontend=='F2':indices,extra=cluster(indices,xy,times,seed)
        record=dict(seed_index=int(seed),seed_physical_id=int(physical_ids[seed]),raw_neighbors=len(original),frontend_neighbors=len(candidate),cluster_neighbors=len(indices),
            raw_neighbor_physical_ids=physical_ids[original].tolist(),candidate_physical_ids=physical_ids[candidate].tolist(),fit_physical_ids=physical_ids[indices].tolist(),
            dropped_T_physical_ids=physical_ids[np.setdiff1d(original,indices)].tolist(),triplet_diagnostic=extra,status=None,
            spatial_ratio=None,time_std_s=None,plane_rms_s=None,plane_gate_s=None,gradient_norm=None,accepted_plane_index=None)
        reason=None
        if len(indices)<10:reason='too_few_neighbors'
        else:
            coords=xy[indices];centered=coords-coords.mean(axis=0);singular=np.linalg.svd(centered,compute_uv=False)
            ratio=float(singular[1]/singular[0]) if singular[0]>1e-6 else 0.;record['spatial_ratio']=ratio
            dt=times[indices]-times[seed];std=float(np.std(dt));record['time_std_s']=std
            if ratio<.08:reason='spatial_plane_degenerate'
            elif std<1e-6:reason='insufficient_local_time_span'
            else:
                design=np.column_stack([centered,np.ones(len(indices))]);fit,_,rank,_=np.linalg.lstsq(design,dt,rcond=None)
                rms=float(np.sqrt(np.mean((design@fit-dt)**2)));gate=max(20e-6,.15*std);norm=float(np.linalg.norm(fit[:2]))
                record.update(plane_rms_s=rms,plane_gate_s=gate,gradient_norm=norm)
                sx,sy=(coords.mean(axis=0)-center)/radius
                if rank<3 or rms>gate:reason='plane_residual_or_rank'
                elif sx*sx+sy*sy>=.9**2:reason='outside_flow_inner_sphere'
                elif norm<1e-8:reason='vanishing_gradient'
                else:
                    record['accepted_plane_index']=len(points);points.append([sx,sy,-np.sqrt(1-sx*sx-sy*sy)]);gradients.append(fit[:2]);centroids.append(coords.mean(axis=0));seed_indices.append(seed)
                    coeff.append(fit);rmses.append(rms);stds.append(std);neighbor_ids.extend(indices.tolist());offsets.append(len(neighbor_ids))
                    dropped_ids.extend(np.setdiff1d(original,indices).tolist());dropped_offsets.append(len(dropped_ids))
        record['status']=reason or 'accepted';records.append(record)
        if reason:failures[reason]+=1
    if len(xy)<100:failures['insufficient_events']=1
    data=dict(points=np.array(points).reshape(-1,3),gradients=np.array(gradients).reshape(-1,2),pixel_centers=np.array(centroids).reshape(-1,2),seed_indices=np.array(seed_indices,int),
        coefficients=np.array(coeff).reshape(-1,3),rms_s=np.array(rmses),time_std_s=np.array(stds),neighbor_offsets=np.array(offsets,int),neighbor_indices=np.array(neighbor_ids,int),
        dropped_offsets=np.array(dropped_offsets,int),dropped_indices=np.array(dropped_ids,int))
    return data,dict(frontend=frontend,accepted_planes=len(points),seeds_attempted=len(seeds),failures=dict(failures),seed_diagnostics=records)
