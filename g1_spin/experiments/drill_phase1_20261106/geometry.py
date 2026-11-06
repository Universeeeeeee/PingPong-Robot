"""Geometry-only diagnostics; all candidates remain without independent truth."""
import csv
import json
from pathlib import Path
import cv2
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from explore import OUT, load, select, picture

SEED={'rotation':(642.,357.),'control':(700.,360.)}
RADIUS={'rotation':60.,'control':60.}

def fit(q,seed,radius):
    # Restrict to vicinity of prior detection; no spin score enters geometry.
    x0,y0=np.floor(np.asarray(seed)-105).astype(int)
    x1,y1=x0+210,y0+210
    q=q[(q['x']>=x0)&(q['x']<x1)&(q['y']>=y0)&(q['y']<y1)]
    im=np.zeros((210,210),np.uint8);im[q['y']-y0,q['x']-x0]=255
    image=cv2.GaussianBlur(im,(5,5),1)
    circles=cv2.HoughCircles(image,cv2.HOUGH_GRADIENT,1,50,param1=80,param2=20,
                            minRadius=round(radius-7),maxRadius=round(radius+7))
    if circles is None:return None
    xy=np.unique(np.column_stack((q['x'],q['y'])),axis=0).astype(float)
    candidates=[]
    for cc in circles[0]:
        c=cc[:2]+[x0,y0]
        # Use equal angular bins in a thin ring, not all texture events.
        delta=xy-c;dist=np.linalg.norm(delta,axis=1)
        theta=np.arctan2(delta[:,1],delta[:,0]);bins=np.floor((theta+np.pi)/(2*np.pi)*36).astype(int)%36
        near=(abs(dist-radius)<4)&~((theta>np.deg2rad(60))&(theta<np.deg2rad(120)))
        pts=np.array([np.median(xy[near&(bins==j)],axis=0) for j in range(36) if np.count_nonzero(near&(bins==j))>=2])
        if len(pts)<15:continue
        opt=least_squares(lambda center:np.linalg.norm(pts-center,axis=1)-radius,c,
                          bounds=(c-5,c+5),loss='soft_l1',f_scale=1.)
        center=opt.x;res=np.linalg.norm(pts-center,axis=1)-radius
        good=abs(res)<=2;support=int(good.sum())
        rms=float(np.sqrt(np.mean(res**2)))
        if support<15 or rms>2.5 or np.any(opt.active_mask):continue
        candidates.append((support,-rms,center,cc[2],rms,len(pts)))
    if not candidates:return None
    support,_,c,hr,rms,n=max(candidates,key=lambda x:x[:2])
    return dict(cx=float(c[0]),cy=float(c[1]),radius_px=radius,hough_radius_px=float(hr),
                ring_rms_px=rms,angular_bins=support,fit_bins=n)

def main():
    summaries={}
    for key in ['rotation','control']:
        e,t0=load(key);duration=(int(e['t'][-1])-t0)*1e-6
        seed=np.array(SEED[key]);rows=[]
        for i,s in enumerate(np.arange(.005,duration-.005,.01)):
            q=select(e,t0,s-.0025,5)
            result=fit(q,seed,RADIUS[key])
            row=dict(time_s=float(s),start_s=float(s-.0025),end_s=float(s+.0025),
                     geometry_valid=result is not None,cx=None,cy=None,radius_px=RADIUS[key],
                     hough_radius_px=None,ring_rms_px=None,angular_bins=0,fit_bins=0)
            if result:
                row.update(result);seed=np.array([result['cx'],result['cy']])
            rows.append(row)
        with (OUT/f'{key}_geometry.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
        good=[r for r in rows if r['geometry_valid']]
        summaries[key]=dict(observations=len(rows),accepted=len(good),first_accepted_s=good[0]['time_s'] if good else None,
                            radius_px=RADIUS[key],note='candidate geometry only, no independent center truth')
        print(key,summaries[key],flush=True)
        times=np.array([r['time_s'] for r in rows])
        fig,axes=plt.subplots(3,1,figsize=(13,8),sharex=True)
        for ax,k in zip(axes,['cx','cy','ring_rms_px']):
            ax.plot(times,[r[k] if r[k] is not None else np.nan for r in rows],'.-',ms=2,lw=.5);ax.set_ylabel(k)
        axes[-1].set_xlabel('Seconds from first CD event');fig.suptitle(f'{key}: fixed-radius geometry candidates; gaps remain gaps')
        fig.tight_layout();fig.savefig(OUT/f'{key}_geometry.png',dpi=140);plt.close(fig)
        review_times=([.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.]
                      if key=='rotation' else [3.5,3.6,3.7,3.8,4.,4.5,5.,5.5,6.,6.5,7.,8.])
        fig,axes=plt.subplots(3,4,figsize=(14,10))
        for ax,s in zip(axes.flat,review_times):
            row=min(rows,key=lambda r:abs(r['time_s']-s))
            center=np.array([row['cx'],row['cy']]) if row['geometry_valid'] else np.array(SEED[key])
            x0,y0=np.floor(center-95).astype(int);box=(x0,y0,x0+190,y0+190)
            q=select(e,t0,row['start_s'],5)
            ax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]])
            if row['geometry_valid']:ax.add_patch(plt.Circle(center,RADIUS[key],fill=False,color='lime',lw=.8))
            ax.set_title(f'{row["time_s"]:.3f}s | '+('candidate' if row['geometry_valid'] else 'NO GEOMETRY'),fontsize=9)
        fig.suptitle(f'{key}: geometry review before spin results')
        fig.tight_layout();fig.savefig(OUT/f'{key}_geometry_review.jpg',dpi=140);plt.close(fig)
        if key=='control':
            # Finer boundary inspection; event activity does not itself establish motion.
            fig,axes=plt.subplots(3,4,figsize=(14,10))
            for ax,s in zip(axes.flat,np.arange(3.50,3.80,.025)):
                box=(570,240,820,490);q=select(e,t0,float(s),5)
                ax.imshow(picture(q,box),extent=[570,820,490,240]);ax.set_title(f'{s:.3f}-{s+.005:.3f}s',fontsize=9)
            fig.suptitle('Control onset audit: raw 5ms samples, no geometry fitting')
            fig.tight_layout();fig.savefig(OUT/'control_onset.jpg',dpi=140);plt.close(fig)
    (OUT/'geometry_summary.json').write_text(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
