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
        engineering_pass=bool(event_values and statistics['median']<=1.5 and statistics['coverage_1p5']>=.5),
        selection='All same-polarity H within4px centroid and±2ms seed; no residual filtering; median per H across all predictions, never best prediction.')

def predict(source_xy,source_times,target_time,center,velocity,radius,omega):
    local=(source_xy-center-(source_times-.005)[:,None]*velocity)/radius
    if not len(local):return np.empty((0,2))
    points=np.column_stack([local,-np.sqrt(np.maximum(0,1-np.sum(local*local,axis=1)))])
    warped=Rotation.from_rotvec((target_time-source_times)[:,None]*np.asarray(omega)).apply(points)
    return warped[:,:2]*radius+center+(target_time-.005)*velocity

def future(source_xy,source_time,source_pol,source_ids,Vxy,Vtime,Vpol,Vids,center,velocity,radius,omega,roi):
    all_dist=[];segments=[]
    for k in range(8):
        low=.006+k*.0005;high=low+.0005;target=(low+high)/2
        keep=(Vtime>=low-1e-12)&(Vtime<high-1e-12)&in_box(Vxy,roi) if roi is not None else np.zeros(len(Vxy),bool)
        targets=Vxy[keep];pol=Vpol[keep];ids=Vids[keep]
        segment=dict(index=k,start_s=low,end_s=high,target_time_s=target,source_ids=source_ids.tolist(),target_ids=ids.tolist(),roi=roi,
            source_count=len(source_xy),target_count=len(targets),status='insufficient_source_or_target',statistics=stats([]))
        if len(source_xy)>=10 and len(targets)>=10:
            prediction=predict(source_xy,source_time,target,center,velocity,radius,omega)
            distances=[];directions=[]
            for polarity in [0,1]:
                pred=prediction[source_pol==polarity];observed=targets[pol==polarity]
                d1=np.minimum(cKDTree(observed).query(pred)[0],3.) if len(pred) and len(observed) else np.full(len(pred),3.)
                d2=np.minimum(cKDTree(pred).query(observed)[0],3.) if len(pred) and len(observed) else np.full(len(observed),3.)
                distances.extend(d1.tolist()+d2.tolist());directions.append(dict(polarity=polarity,pred_to_target=stats(d1),target_to_pred=stats(d2)))
            segment.update(status='evaluable',statistics=stats(distances),directions=directions,predicted_xy=prediction.tolist())
            all_dist.extend(distances)
        segments.append(segment)
    return dict(segments=segments,evaluable_segments=sum(s['status']=='evaluable' for s in segments),statistics=stats(all_dist),
        method='Polarity bidirectional nearest-neighbor Chamfer, all sources/targets in fixed ROI, cap3px; no V correspondence fitting')

def compare_future(model,zero):
    median=model['statistics']['median'];baseline=zero['statistics']['median'];improvement=None if median is None or baseline is None or baseline<=1e-12 else 1-median/baseline
    evaluable=model['evaluable_segments']>=4 and zero['evaluable_segments']>=4 and improvement is not None
    return dict(evaluable=evaluable,model_distance_median_px=median,zero_distance_median_px=baseline,relative_zero_improvement=improvement,
        engineering_pass=bool(evaluable and median<=1.5 and improvement>=.2))
