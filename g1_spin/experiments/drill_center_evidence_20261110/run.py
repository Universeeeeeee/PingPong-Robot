"""Spatially separated geometry evidence from raw events, no spin estimator calls."""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT.parent/'drill_geometry_diagnosis_20261110'))
s=importlib.util.spec_from_file_location('previous',OUT.parent/'drill_geometry_diagnosis_20261110/run.py');prev=importlib.util.module_from_spec(s);s.loader.exec_module(prev)
old=prev.old;p2=prev.p2;ROOT=prev.ROOT
CASES=prev.ENDS


def fit_circle(points,anchor,radius,minimum):
    if len(points)<minimum:return dict(valid=False,reason='insufficient_bins',bins=len(points))
    fit=least_squares(lambda c:np.linalg.norm(points-c,axis=1)-radius,anchor,bounds=(anchor-10,anchor+10),loss='soft_l1',f_scale=1.)
    delta=fit.x-points;dist=np.linalg.norm(delta,axis=1);jac=delta/dist[:,None]
    singular=np.linalg.svd(jac,compute_uv=False);condition=float(singular[0]/singular[-1])
    rms=float(np.sqrt(np.mean((dist-radius)**2)))
    reason='bound' if np.any(fit.active_mask) else 'residual' if rms>2.5 else 'condition' if condition>10 else None
    return dict(valid=reason is None,reason=reason,bins=len(points),center=fit.x.tolist(),rms_px=rms,condition=condition)


def measure(raw,origin,mid,duration,anchor,radius,direction):
    pixels=np.column_stack([raw['x'],raw['y']]).astype(float);delta=pixels-anchor
    unique=np.unique(pixels,axis=0);u=unique-anchor;angle=np.arctan2(u[:,1],u[:,0]);dist=np.linalg.norm(u,axis=1)
    d=np.array(direction);perp=np.array([-d[1],d[0]]);sector=np.abs(np.arctan2(np.sin(angle-np.arctan2(d[1],d[0])),np.cos(angle-np.arctan2(d[1],d[0]))))
    near=(abs(dist-radius)<=8)&(sector>np.pi/6);bins=np.floor((angle+np.pi)*36/(2*np.pi)).astype(int)%36
    pts=[];bin_ids=[];timing=[]
    raw_angles=np.arctan2(delta[:,1],delta[:,0]);raw_bins=np.floor((raw_angles+np.pi)*36/(2*np.pi)).astype(int)%36
    for b in range(36):
        mask=near&(bins==b)
        if mask.sum()>=2:
            pts.append(np.median(unique[mask],axis=0));bin_ids.append(b)
            rmask=(raw_bins==b)&(abs(np.linalg.norm(delta,axis=1)-radius)<=8)
            timing.append(float(np.median(raw['t'][rmask].astype(np.int64)-origin-mid)))
    points=np.asarray(pts).reshape(-1,2);sides=(points-anchor)@perp
    fits={key:fit_circle(points[mask],anchor,radius,minimum) for key,mask,minimum in [
          ('all',np.ones(len(points),bool),12),('side_a',sides<0,6),('side_b',sides>=0,6)]}
    # Spatially separate shaft support with a fixed coordinate origin and evaluation plane.
    along=delta@d;across=delta@perp;shaft_pts=[]
    for lo in np.arange(radius+15,radius+170,12):
        mask=(along>=lo)&(along<lo+12)&(abs(across)<25)
        if mask.sum()>=3:shaft_pts.append([np.median(along[mask]),np.median(across[mask])])
    shaft=dict(valid=False,reason='insufficient_support',points=len(shaft_pts))
    if len(shaft_pts)>=6:
        sp=np.asarray(shaft_pts);span=float(np.ptp(sp[:,0]))
        if span>=70:
            fit=least_squares(lambda c:c[0]*sp[:,0]+c[1]-sp[:,1],[0.,0.],loss='soft_l1')
            rms=float(np.sqrt(np.mean((fit.x[0]*sp[:,0]+fit.x[1]-sp[:,1])**2)))
            shaft.update(valid=rms<=3,reason=None if rms<=3 else 'residual',rms_px=rms,span_px=span,slope=float(fit.x[0]),
                intercept=float(fit.x[1]),offset_px=float(fit.x[0]*(radius+50)+fit.x[1]))
    return dict(fits=fits,shaft=shaft,ring_points=points.tolist(),bin_ids=bin_ids,bin_median_time_offset_us=timing)


def main():
    if (OUT/'measurements.json').exists():raise FileExistsError('Preserve output')
    paths=[OUT/'run.py',OUT/'protocol.md',p2.OUT/'geometry.json',p2.OLD/'rotation_geometry.csv',p2.OLD/'control_geometry.csv',
           p2.OUT/'run.py',old.OUT/'run.py',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py']
    hashes={str(p.relative_to(ROOT)):old.digest(p) for p in paths};old.save(OUT/'freeze.json',hashes)
    (OUT/'raw').mkdir(exist_ok=True);records=[];anchors=[]
    for run,ends in CASES.items():
        folder=p2.ARCHIVE/'decoded'/p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']]
        origin=int(parts[0]['t'][0]);ts,xy=p2.observations(run);radius=prev.prev.RADIUS[run]
        theta=np.deg2rad(35 if run==3 else 90);direction=np.array([np.cos(theta),np.sin(theta)])
        for end in ends:
            anchor=p2.centers(ts,xy,[end-10000])[0];key=f'run{run}_{end:08d}'
            anchors.append(dict(case=key,run=run,end_us=end,center=anchor.tolist(),radius=radius,direction=direction.tolist()))
            fig,axes=plt.subplots(2,9,figsize=(20,6))
            for col,mid in enumerate(range(end-50000,end+30001,10000)):
                try:center_old=p2.centers(ts,xy,[mid])[0].tolist()
                except ValueError:center_old=None
                for row,duration in enumerate([1000,5000]):
                    raw,ids=old.read_window(parts,origin,mid-duration//2,mid+duration//2)
                    # Persist full geometry ROI, never only preselected ring points.
                    pix=np.column_stack([raw['x'],raw['y']]);mask=np.linalg.norm(pix-anchor,axis=1)<=radius+190
                    raw=raw[mask];ids=ids[mask]
                    path=OUT/'raw'/f'{key}_{mid:08d}_{duration}us.npz'
                    np.savez_compressed(path,events=raw,ids=ids,origin_us=origin)
                    result=measure(raw,origin,mid,duration,anchor,radius,direction)
                    records.append(dict(case=key,run=run,end_us=end,mid_us=mid,duration_us=duration,old_center=center_old,
                                        input_file=str(path.relative_to(OUT)),input_sha256=old.digest(path),events=len(raw),**result))
                    ax=axes[row,col];box=np.r_[np.floor(anchor-90).astype(int),np.ceil(anchor+90).astype(int)]
                    ax.imshow(p2.picture(raw,box),extent=[box[0],box[2],box[3],box[1]])
                    ax.add_patch(plt.Circle(anchor,radius,fill=False,color='gray',lw=.6,ls='--'))
                    if center_old:ax.plot(*center_old,'wx',ms=4)
                    for label,color in [('side_a','lime'),('side_b','magenta')]:
                        f=result['fits'][label]
                        if f['valid']:
                            ax.add_patch(plt.Circle(f['center'],radius,fill=False,color=color,lw=.6));ax.plot(*f['center'],'+',color=color,ms=4)
                    ax.set_title(f'{mid/1e6:.3f}s / {duration/1000:g}ms\nA:{result["fits"]["side_a"]["bins"]} B:{result["fits"]["side_b"]["bins"]}',fontsize=8)
                    ax.set_xticks([]);ax.set_yticks([])
            fig.suptitle(f'{key}: raw polarity events | gray fixed prior; white old center; green/magenta separate halves',fontsize=12)
            fig.tight_layout();fig.savefig(OUT/f'{key}_raw.png',dpi=130);plt.close(fig)
            print(key,flush=True)
        del parts
    old.save(OUT/'anchors.json',anchors);old.save(OUT/'measurements.json',records)
    assert len(records)==144 and all(old.digest(ROOT/p)==h for p,h in hashes.items())
    print('144 exposures saved; frozen files unchanged')

if __name__=='__main__':main()
