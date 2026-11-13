"""Outermost locally supported event bands; no invented/projected pixels."""
import numpy as np
from scipy.spatial import cKDTree


def extract(raw, anchor, radius, direction):
    points=np.unique(np.column_stack([raw['x'],raw['y']]),axis=0).astype(float)
    delta=points-anchor;radii=np.linalg.norm(delta,axis=1);theta=np.arctan2(delta[:,1],delta[:,0])
    shaft_angle=np.arctan2(direction[1],direction[0])
    separation=np.abs(np.arctan2(np.sin(theta-shaft_angle),np.cos(theta-shaft_angle)))
    mask=(abs(radii-radius)<=12)&(separation>np.pi/6)
    points,radii,theta=points[mask],radii[mask],theta[mask]
    if not len(points):return dict(points=[],bins=[],reason='no_ring_support',boundary_bins=0)
    # Repeated events at a hot pixel cannot manufacture distinct spatial support.
    neighbors=cKDTree(points).query_ball_point(points,2.5,return_length=True)
    good=neighbors>=3
    points,radii,theta=points[good],radii[good],theta[good]
    bins=np.floor((theta+np.pi)*36/(2*np.pi)).astype(int)%36
    candidates={};boundary=0
    for b in range(36):
        ix=np.flatnonzero(bins==b)
        if len(ix)<3:continue
        chosen=None
        for upper in sorted(radii[ix],reverse=True):
            band=ix[(radii[ix]>=upper-2.5)&(radii[ix]<=upper)]
            # Support must extend tangentially, not just a radial streak.
            if len(band)>=3 and np.ptp(theta[band])>=.035:
                chosen=band;break
        if chosen is None:continue
        if radii[chosen].max()>=radius+11:
            boundary+=1;continue
        median=np.median(points[chosen],axis=0)
        actual=chosen[np.argmin(np.linalg.norm(points[chosen]-median,axis=1))]
        candidates[b]=(points[actual],float(np.median(radii[chosen])))
    # Require a neighboring angular bin at a compatible radius, without filling gaps.
    kept=[];kept_bins=[]
    for b,(p,r) in candidates.items():
        if any(n in candidates and abs(candidates[n][1]-r)<=3 for n in [(b-1)%36,(b+1)%36]):
            kept.append(p.tolist());kept_bins.append(b)
    return dict(points=kept,bins=kept_bins,reason=None if kept else 'no_supported_arc',boundary_bins=boundary)
