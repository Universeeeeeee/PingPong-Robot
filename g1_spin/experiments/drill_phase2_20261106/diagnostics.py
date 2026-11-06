"""Raw geometry evidence and velocity coupling; never chooses geometry by spin."""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from run import OUT,OLD,SPEED,load,select,picture,table,write_json,cmax,digest

def shaft(raw,center,radius,direction):
    """External line in a fixed direction corridor, with independent support gates."""
    direction=np.asarray(direction,float);direction/=np.linalg.norm(direction)
    perpendicular=np.array([-direction[1],direction[0]])
    xy=np.column_stack([raw['x'],raw['y']]).astype(float)-center
    along=xy@direction;across=xy@perpendicular;pts=[]
    for low in np.arange(radius+15,radius+170,12):
        keep=(along>=low)&(along<low+12)&(abs(across)<25)
        if keep.sum()>=3:pts.append([np.median(along[keep]),np.median(across[keep])])
    row=dict(status='insufficient_support',points=len(pts),span_px=None,rms_px=None,direction=None,
        slope=None,center_relative_intercept_px=None,x_at_global_y450_px=None)
    if len(pts)<6:return row
    pts=np.array(pts);span=float(np.ptp(pts[:,0]));row['span_px']=span
    if span<70:return row
    opt=least_squares(lambda p:p[0]*pts[:,0]+p[1]-pts[:,1],[0.,0.],loss='soft_l1')
    rms=float(np.sqrt(np.mean((opt.x[0]*pts[:,0]+opt.x[1]-pts[:,1])**2)))
    vector=direction+opt.x[0]*perpendicular;vector/=np.linalg.norm(vector)
    linepoint=center+opt.x[1]*perpendicular
    x450=float(linepoint[0]+(450-linepoint[1])*vector[0]/vector[1]) if abs(vector[1])>1e-6 else None
    row.update(status='accepted_projection_candidate' if rms<=3 else 'large_line_residual',rms_px=rms,direction=vector.tolist(),
        slope=float(opt.x[0]),center_relative_intercept_px=float(opt.x[1]),x_at_global_y450_px=x450)
    return row

def contour(raw,center,radius):
    """Candidate-relative opposing arc support; diagnostic, not a center detector."""
    xy=np.unique(np.column_stack([raw['x'],raw['y']]),axis=0).astype(float)
    delta=xy-center;dist=np.linalg.norm(delta,axis=1);theta=np.arctan2(delta[:,1],delta[:,0])
    bins=np.floor((theta+np.pi)/(2*np.pi)*36).astype(int)%36
    keep=(abs(dist-radius)<4)&~((theta>np.deg2rad(60))&(theta<np.deg2rad(120)))
    ids=[j for j in range(36) if np.count_nonzero(keep&(bins==j))>=2]
    pts=np.array([np.median(xy[keep&(bins==j)],axis=0) for j in ids])
    metrics=dict(support_bins=ids,left_bins=0,right_bins=0,left_center=None,right_center=None)
    if not len(pts):return metrics,pts
    for name,mask in [('left',pts[:,0]<center[0]-radius*.35),('right',pts[:,0]>center[0]+radius*.35)]:
        side=pts[mask];metrics[name+'_bins']=len(side)
        if len(side)>=6:
            fit=least_squares(lambda c:np.linalg.norm(side-c,axis=1)-radius,center,loss='soft_l1',bounds=(center-10,center+10))
            if not np.any(fit.active_mask):metrics[name+'_center']=fit.x.tolist()
    return metrics,pts

def main():
    metadata=[];coupling=[];(OUT/'diagnostics').mkdir(exist_ok=True)
    for run,starts in SPEED.items():
        e,t0=load(run);oldkey='rotation' if run==1 else 'control'
        geom=list(csv.DictReader((OLD/f'{oldkey}_geometry.csv').open()))
        gt=np.array([float(g['time_s']) for g in geom])
        for begin in starts:
            name=f'run{run}_{round(begin*1e6):08d}_ablation'
            d=np.load(OUT/'windows'/f'{name}.npz');ref=begin+.005
            context,first=select(e,t0,round(begin*1e6)-25000,60000)
            context_path=OUT/'diagnostics'/f'{name}_raw_context.npz'
            np.savez_compressed(context_path,ids=np.arange(first,first+len(context)),events=context,origin_us=t0)
            fig,axes=plt.subplots(2,11,figsize=(27,6.5));record=[]
            midpoint=d['A0'][0];box=np.r_[np.floor(midpoint-[85,85]).astype(int),np.ceil(midpoint+[85,155]).astype(int)]
            for j,offset in enumerate(range(-25,26,5)):
                nominal=ref+offset*.001;g=geom[int(np.argmin(abs(gt-nominal)))];valid=g['geometry_valid']=='True'
                center=np.array([float(g['cx']),float(g['cy'])]) if valid else None
                for i,duration in enumerate([1000,5000]):
                    q,_=select(e,t0,round(nominal*1e6)-duration//2,duration)
                    median=float(np.median(q['t']-t0)*1e-6) if len(q) else None
                    ax=axes[i,j];ax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]])
                    metrics=None;line=None;roi_median=None;support_median=None;old_observation_median=None
                    if valid:
                        ax.add_patch(plt.Circle(center,60,fill=False,color='lime',lw=.6));ax.plot(*center,'+',color='lime',ms=4)
                        metrics,pts=contour(q,center,60.)
                        if len(pts):ax.scatter(pts[:,0],pts[:,1],s=3,c='magenta')
                        line=shaft(q,center,60.,[0,1])
                        roi=q[(abs(q['x'].astype(float)-center[0])<105)&(abs(q['y'].astype(float)-center[1])<105)]
                        if len(roi):roi_median=float(np.median(roi['t']-t0)*1e-6)
                        delta=np.column_stack([q['x'],q['y']]).astype(float)-center
                        support=(abs(np.linalg.norm(delta,axis=1)-60)<4)
                        theta=np.arctan2(delta[:,1],delta[:,0]);support&=~((theta>np.deg2rad(60))&(theta<np.deg2rad(120)))
                        if support.any():support_median=float(np.median(q['t'][support]-t0)*1e-6)
                        oldq,_=select(e,t0,round(float(g['start_s'])*1e6),5000)
                        oldq=oldq[(abs(oldq['x'].astype(float)-center[0])<105)&(abs(oldq['y'].astype(float)-center[1])<105)]
                        if len(oldq):old_observation_median=float(np.median(oldq['t']-t0)*1e-6)
                    h=g['hough_radius_px'] if valid else '-';bins=len(metrics['support_bins']) if metrics else 0
                    title=f'{nominal:.3f}s {duration/1000:g}ms\nmed {median:.4f}' if median is not None else 'empty'
                    ax.set_title(title+f'\ng{float(g["time_s"]):.3f} H{h[:5]} b{bins}',fontsize=7);ax.tick_params(labelsize=6)
                    record.append(dict(nominal_s=nominal,duration_us=duration,event_median_s=median,
                        roi_event_median_s=roi_median,ring_support_event_median_s=support_median,
                        old_observation_roi_event_median_s=old_observation_median,
                        old_geometry_time_s=float(g['time_s']),old_geometry_valid=valid,old_center=None if center is None else center.tolist(),
                        hough_radius_px=float(h) if valid else None,angular_support=metrics,external_shaft=line,raw_count=len(q)))
            fig.suptitle(f'{name} | raw 1ms / 5ms, nominal vs median times | green=old candidate, magenta=ring support | no spin warp')
            fig.tight_layout();fig.savefig(OUT/'diagnostics'/f'{name}.png',dpi=140);plt.close(fig)
            metadata.append(dict(window=name,raw_context_file=str(context_path.relative_to(OUT)),raw_context_sha256=digest(context_path),
                raw_context_start_us=round(begin*1e6)-25000,raw_context_end_us=round(begin*1e6)+35000,
                raw_context_events=len(context),snapshots=record,
                warning='Opposing arcs are selected near old candidates; neither they nor shaft lines are independent center truth.'))
            dt=(d['time_us']-d['time_us'][0])*1e-6;span=float(np.ptp(dt))
            c1=d['A0'];c2=d['A1'];velocity=(c2[-1]-c2[0])/span
            near=np.array([-velocity[1]/60,velocity[0]/60,0.])
            s1=cmax.lift(*(d['xy']-c1).T,0,0,60);s2=cmax.lift(*(d['xy']-c2).T,0,0,60)
            # First-order sphere displacement versus rotational displacement, all three coordinates.
            tref=(d['time_us']-(round(ref*1e6)))*1e-6
            x,y,z=s1.T;zero=np.zeros(len(x))
            matrix=np.stack([np.stack([zero,z,-y],axis=1),np.stack([-z,zero,x],axis=1),np.stack([y,-x,zero],axis=1)],axis=1)
            matrix*=tref[:,None,None]
            predicted,_,rank,_=np.linalg.lstsq(matrix.reshape(-1,3),(s2-s1).reshape(-1),rcond=None)
            coupling.append(dict(window=name,run=run,start_us=round(begin*1e6),event_span_s=span,
                center_delta_px=(c2[-1]-c2[0]).tolist(),center_velocity_px_s=velocity.tolist(),
                near_axis_delta_omega=near.tolist(),full_sphere_linear_delta_omega=predicted.tolist(),
                full_sphere_rank=int(rank),mean_z=float(s1[:,2].mean()),
                note='Algebraic coupling prediction for a velocity hypothesis, not measured physical center motion.'))
            print('diagnostics',name,flush=True)
    write_json(OUT/'raw_diagnostics.json',metadata);write_json(OUT/'velocity_predictions.json',coupling)
    table(OUT/'velocity_predictions.csv',coupling)

if __name__=='__main__':main()
