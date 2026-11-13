"""Heldout event evaluation; never used by plane extraction or optimizer."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

def stats(values):
    a=np.asarray(values,float)
    return dict(count=len(a),median=float(np.median(a)) if len(a) else None,mean=float(np.mean(a)) if len(a) else None,
        p95=float(np.percentile(a,95)) if len(a) else None,coverage_1p5=float(np.mean(a<=1.5)) if len(a) else None)

def regions(raw_xy,center):
    cells=np.floor((raw_xy-center)/24).astype(int);unique,counts=np.unique(cells,axis=0,return_counts=True)
    order=sorted(range(len(unique)),key=lambda i:(-counts[i],int(unique[i,0]),int(unique[i,1])));chosen=[]
    for i in order:
        if counts[i]<20:continue
        candidate=unique[i];middle=center+(candidate+.5)*24
        if any(np.linalg.norm(middle-np.array(r['grid_center']))<24-1e-9 for r in chosen):continue
        mask=np.all(cells==candidate,axis=1);pixels=raw_xy[mask]
        chosen.append(dict(region=len(chosen),grid_cell=candidate.tolist(),grid_center=middle.tolist(),count=int(mask.sum()),
            bbox=[float(pixels[:,0].min()),float(pixels[:,1].min()),float(pixels[:,0].max()),float(pixels[:,1].max())]))
        if len(chosen)==6:break
    return chosen

def in_box(xy,box,pad=0):return (xy[:,0]>=box[0]-pad)&(xy[:,0]<=box[2]+pad)&(xy[:,1]>=box[1]-pad)&(xy[:,1]<=box[3]+pad)

def roi_from_source(raw_xy):
    if not len(raw_xy):return None
    return [max(0.,float(raw_xy[:,0].min())-15),max(0.,float(raw_xy[:,1].min())-15),min(1279.,float(raw_xy[:,0].max())+15),min(719.,float(raw_xy[:,1].max())+15)]

def normal_support(planes,training_raw,selected_regions):
    counts=[];directions=[]
    for r in selected_regions:
        seedxy=training_raw[planes['seed_indices']];mask=in_box(seedxy,r['bbox']);counts.append(dict(region=r['region'],accepted_planes=int(mask.sum())))
        if mask.any():
            gradient=planes['gradients'][mask];gradient=gradient/np.linalg.norm(gradient,axis=1)[:,None]
            # Median orientation tensor avoids directed sign cancellation, without fitting truth labels.
            tensor=gradient.T@gradient;_,vectors=np.linalg.eigh(tensor);directions.append(vectors[:,-1])
    angles=[float(np.degrees(np.arccos(np.clip(abs(a@b),-1,1)))) for i,a in enumerate(directions) for b in directions[i+1:]]
    regions_count=sum(c['accepted_planes']>0 for c in counts)
    return dict(plane_region_counts=counts,regions_with_planes=regions_count,max_unoriented_normal_separation_deg=max(angles,default=None),
        sufficient=regions_count>=3 and max(angles,default=0)>=20)

def local_holdout(planes,Txy,Ttime,Tpol,Hxy,Htime,Hpol,Hids):
    tree=cKDTree(Hxy);rows=[];by_event={}
    for k,pixel in enumerate(planes['pixel_centers']):
        seed=planes['seed_indices'][k];indices=np.asarray(tree.query_ball_point(pixel,4),int)
        indices=indices[(abs(Htime[indices]-Ttime[seed])<=.002)&(Hpol[indices]==Tpol[seed])]
        gradient=planes['gradients'][k];norm=np.linalg.norm(gradient)
        distance=abs(Htime[indices]-Ttime[seed]-(Hxy[indices]-pixel)@gradient-planes['coefficients'][k,2])/norm
        for idx,value in zip(indices,distance):by_event.setdefault(int(Hids[idx]),[]).append(float(value))
        rows.append(dict(plane=k,holdout_ids=Hids[indices].tolist(),normal_distance_px=distance.tolist(),statistics=stats(distance)))
    values=[v for row in rows for v in row['normal_distance_px']];event_values=[float(np.median(v)) for v in by_event.values()]
    statistics=stats(event_values)
    return dict(status='evaluable' if event_values else 'no_local_H_support',all_plane_H_pairs=rows,pair_statistics=stats(values),
        event_statistics=statistics,unique_H_events=len(by_event),H_total=len(Hids),
        predicted_H_ids=sorted(by_event),total_H_prediction_coverage=len(by_event)/len(Hids) if len(Hids) else None,
        engineering_pass=bool(len(by_event)>=20 and len(by_event)>=.2*len(Hids) and statistics['median']<=1.5 and statistics['coverage_1p5']>=.5),
        selection='All same-polarity H within4px centroid and±2ms seed; no residual filtering; median per H across all predictions, never best prediction.')

def predict(source_xy,source_times,target_time,center,velocity,radius,omega):
    local=(source_xy-center-(source_times-.005)[:,None]*velocity)/radius
    if not len(local):return dict(status='ok',xy=np.empty((0,2)),visible=np.empty(0,bool),invalid_source_indices=[])
    q=np.sum(local*local,axis=1);invalid=~np.isfinite(q)|(q>=1)
    if invalid.any():return dict(status='geometry_invalid',xy=None,visible=None,invalid_source_indices=np.flatnonzero(invalid).tolist())
    points=np.column_stack([local,-np.sqrt(1-q)])
    warped=Rotation.from_rotvec((target_time-source_times)[:,None]*np.asarray(omega)).apply(points)
    return dict(status='ok',xy=warped[:,:2]*radius+center+(target_time-.005)*velocity,visible=warped[:,2]<0,invalid_source_indices=[])

def future(source_xy,source_time,source_pol,source_ids,Vxy,Vtime,Vpol,Vids,center,velocity,radius,omega,roi):
    all_dist=[];all_forward=[];all_reverse=[];segments=[]
    for k in range(8):
        low=.006+k*.0005;high=low+.0005;target=(low+high)/2
        keep=(Vtime>=low-1e-12)&(Vtime<high-1e-12)&in_box(Vxy,roi) if roi is not None else np.zeros(len(Vxy),bool)
        targets=Vxy[keep];pol=Vpol[keep];ids=Vids[keep]
        segment=dict(index=k,start_s=low,end_s=high,target_time_s=target,source_ids=source_ids.tolist(),target_ids=ids.tolist(),roi=roi,
            source_count=len(source_xy),target_count=len(targets),status='insufficient_source_or_target',statistics=stats([]))
        if len(source_xy)>=10 and len(targets)>=10:
            propagated=predict(source_xy,source_time,target,center,velocity,radius,omega)
            if propagated['status']!='ok':
                segment.update(status='geometry_invalid',invalid_source_ids=source_ids[propagated['invalid_source_indices']].tolist());segments.append(segment);continue
            prediction=propagated['xy'];visible=propagated['visible'];forward=[];reverse=[]
            distances=[];directions=[]
            for polarity in [0,1]:
                matching=source_pol==polarity;pred=prediction[matching];seen=visible[matching];observed=targets[pol==polarity]
                d1=np.full(len(pred),3.)
                if seen.any() and len(observed):d1[seen]=np.minimum(cKDTree(observed).query(pred[seen])[0],3.)
                d2=np.minimum(cKDTree(pred[seen]).query(observed)[0],3.) if seen.any() and len(observed) else np.full(len(observed),3.)
                distances.extend(d1.tolist()+d2.tolist());directions.append(dict(polarity=polarity,pred_to_target=stats(d1),target_to_pred=stats(d2)))
                forward.extend(d1.tolist());reverse.extend(d2.tolist())
            segment.update(status='evaluable',statistics=stats(distances),directions=directions,predicted_xy=prediction.tolist(),
                invisible_source_ids=source_ids[~visible].tolist(),invisible_source_count=int((~visible).sum()),
                pred_to_target=stats(forward),target_to_pred=stats(reverse),equal_direction_mean_chamfer=(np.mean(forward)+np.mean(reverse))/2)
            all_dist.extend(distances);all_forward.extend(forward);all_reverse.extend(reverse)
        segments.append(segment)
    return dict(segments=segments,evaluable_segments=sum(s['status']=='evaluable' for s in segments),statistics=stats(all_dist),
        pred_to_target=stats(all_forward),target_to_pred=stats(all_reverse),geometry_invalid_segments=sum(s['status']=='geometry_invalid' for s in segments),
        equal_direction_mean_chamfer=float((np.mean(all_forward)+np.mean(all_reverse))/2) if all_forward and all_reverse else None,
        method='Polarity bidirectional nearest-neighbor Chamfer, all sources/targets in fixed ROI, cap3px; no V correspondence fitting')

def compare_future(model,zero):
    directions={}
    for name in ['pred_to_target','target_to_pred']:
        median=model[name]['median'];baseline=zero[name]['median'];improvement=None if median is None or baseline is None or baseline<=1e-12 else 1-median/baseline
        directions[name]=dict(model_median_px=median,baseline_median_px=baseline,relative_improvement=improvement,
            passed=bool(improvement is not None and median<=1.5 and improvement>=.2))
    median=model['statistics']['median'];baseline=zero['statistics']['median'];improvement=None if median is None or baseline is None or baseline<=1e-12 else 1-median/baseline
    evaluable=model['evaluable_segments']>=4 and zero['evaluable_segments']>=4 and not model['geometry_invalid_segments'] and not zero['geometry_invalid_segments']
    return dict(evaluable=evaluable,model_distance_median_px=median,zero_distance_median_px=baseline,relative_zero_improvement=improvement,
        model_equal_direction_mean_chamfer=model['equal_direction_mean_chamfer'],baseline_equal_direction_mean_chamfer=zero['equal_direction_mean_chamfer'],
        direction_comparisons=directions,engineering_pass=bool(evaluable and all(v['passed'] for v in directions.values())))
