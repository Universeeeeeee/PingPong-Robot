"""Cross-camera ball-position alignment diagnostic, no calibrated spin truth."""
import json,sqlite3
from pathlib import Path
import numpy as np,cv2
from scipy.optimize import least_squares
from audit_remote import OUT,SESSION,parse_image,save

def detect(im,kind):
    if kind=='event':
        gray=cv2.cvtColor(im,cv2.COLOR_BGR2GRAY);gray=(gray>0).astype(np.uint8)*255
        x0,y0,x1,y1=180,100,850,640;minr,maxr=38,80;threshold=20
    else:
        gray=im;x0,y0,x1,y1=220,70,460,300;minr,maxr=9,22;threshold=12
    crop=cv2.GaussianBlur(gray[y0:y1,x0:x1],(5,5),1)
    circles=cv2.HoughCircles(crop,cv2.HOUGH_GRADIENT,1,30,param1=60,param2=threshold,minRadius=minr,maxRadius=maxr)
    candidates=[]
    if circles is not None:
        for cx,cy,r in circles[0]:
            cx+=x0;cy+=y0
            # Distinct angle support around the candidate perimeter.
            yy,xx=np.nonzero(gray[max(0,int(cy-r-4)):int(cy+r+5),max(0,int(cx-r-4)):int(cx+r+5)]> (0 if kind=='event' else 45))
            xx=xx+max(0,int(cx-r-4));yy=yy+max(0,int(cy-r-4))
            near=abs(np.hypot(xx-cx,yy-cy)-r)<(2 if kind=='event' else 1.5)
            bins=np.unique((np.floor((np.arctan2(yy[near]-cy,xx[near]-cx)+np.pi)*24/(2*np.pi)).astype(int))%24)
            if len(bins)>=14:candidates.append(dict(x=float(cx),y=float(cy),r=float(r),support=len(bins)))
    # No nearest expected trajectory scoring; multiple choices retain ambiguity.
    candidates.sort(key=lambda x:x['support'],reverse=True)
    return candidates

def main():
    event=[];ir=[]
    for path in sorted((OUT/'event_samples').glob('*.png')):
        t=int(path.stem)*.5;im=cv2.imread(str(path));cc=detect(im,'event');event.append(dict(t=t,candidates=cc))
    db=sqlite3.connect('file:'+str(SESSION/'realsense.db3')+'?mode=ro',uri=True);rough=json.loads((OUT/'rough_alignment.json').read_text())['corresponding_db_s'];iid=110
    for target in np.arange(rough-3,min(262.5,rough+58),.1):
        row=db.execute('select timestamp,data from messages where topic_id=? and timestamp>=? order by timestamp limit 1',(iid,round(target*1e9))).fetchone()
        if row:
            im,_=parse_image(row[1]);ir.append(dict(t=row[0]/1e9,candidates=detect(im,'ir')))
    save('ball_tracking_candidates.json',dict(event=event,ir=ir,rough=rough))
    print('single candidates event',sum(len(r['candidates'])==1 for r in event),'ir',sum(len(r['candidates'])==1 for r in ir),flush=True)
    # Exactly one circle only: avoid silently choosing from competing detections.
    ee=[r for r in event if len(r['candidates'])==1];ii=[r for r in ir if len(r['candidates'])==1]
    et=np.array([r['t'] for r in ee]);exy=np.array([[r['candidates'][0][k] for k in ['x','y']] for r in ee]);it=np.array([r['t'] for r in ii]);ixy=np.array([[r['candidates'][0][k] for k in ['x','y']] for r in ii])
    save('ball_tracking_single.json',dict(event_times=et.tolist(),event_xy=exy.tolist(),ir_times=it.tolist(),ir_xy=ixy.tolist()))
    # Joint affine spatial fit is only a diagnostic; camera depth/parallax not calibrated.
    # Chronological split: first 30s fit, last period validation, no tuning to final reference.
    scans=[]
    for lag in np.arange(rough-2,rough+2.001,.02):
        target=et+lag;right=np.searchsorted(it,target);safe=(right>0)&(right<len(it));valid=np.flatnonzero(safe)
        valid=valid[(it[right[valid]]-it[right[valid]-1])<=.25]
        train=valid[(et[valid]>=3)&(et[valid]<30)];test=valid[(et[valid]>=30)&(et[valid]<52)]
        if len(train)<12:continue
        y=np.column_stack([np.interp(target,it,ixy[:,d]) for d in range(2)])
        X=np.column_stack([exy,np.ones(len(exy))]);initial=np.linalg.lstsq(X[train],y[train],rcond=None)[0]
        fit=least_squares(lambda z:(X[train]@z.reshape(3,2)-y[train]).ravel(),initial.ravel(),loss='soft_l1',f_scale=2)
        residual=np.linalg.norm(X@fit.x.reshape(3,2)-y,axis=1)
        scans.append(dict(lag_s=float(lag),train_n=len(train),test_n=len(test),train_median_px=float(np.median(residual[train])),train_mean_px=float(np.mean(residual[train])),test_median_px=float(np.median(residual[test])) if len(test) else None,affine=fit.x.reshape(3,2).tolist()))
    save('alignment_scan.json',scans)
    if scans:
        best=min(scans,key=lambda x:x['train_mean_px']);save('alignment_diagnostic.json',dict(best=best,warning='Position-derived time/spatial fit using uncalibrated cameras; not hardware synchronization or spin GT. Temporal support changes with lag. No drift estimate.'))
        print(best,flush=True)

if __name__=='__main__':main()
