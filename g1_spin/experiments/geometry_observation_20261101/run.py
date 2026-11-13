"""Offline geometry observation experiment. No estimator, hardware, or network access.

Run from repository root: g1_spin/.venv/bin/python g1_spin/experiments/geometry_observation_20261101/run.py
"""
import csv
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DATA = ROOT / 'g1_spin/data/decoded/local_analysis_222403'
SEGMENTS = {'translation': (.30, .60), 'early_turn': (1.60, 1.90),
            'middle': (6.20, 6.50), 'late': (11.20, 11.50)}
ORIGIN = np.array([556., 224.])
ANGLES = np.linspace(-np.pi, 0, 121)
RADII = np.arange(30., 76., .5)
XX = ORIGIN[0] + np.cos(ANGLES[:, None]) * RADII
YY = ORIGIN[1] + np.sin(ANGLES[:, None]) * RADII
POLICIES = {'scaled_0.6': .6, 'fixed_0.6': .6, 'scaled_0.3': .3, 'scaled_1.0': 1.}
BOUNDS = ([525., 200., 35.], [590., 250., 70.])


def write_csv(name, rows):
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: round(float(v), 6) if isinstance(v, (float, np.floating)) else v
                             for k, v in row.items()})


def longest_run(mask):
    return max((len(x) for x in ''.join('1' if b else '0' for b in mask).split('0')), default=0)


def candidates(events, exposure_ms, policy):
    # Crop contains all sample locations plus >4 sigma blur padding. Coordinates stay global.
    use = events[(events['x'] >= 475) & (events['x'] < 638) &
                 (events['y'] >= 143) & (events['y'] < 237)]
    counts = np.bincount((use['y'].astype(int)-143)*163 + use['x'].astype(int)-475,
                         minlength=94*163).reshape(94, 163)
    smoothed = gaussian_filter(counts.astype(float), 1)
    profiles = map_coordinates(smoothed, [YY-143, XX-475], order=1)
    threshold = POLICIES[policy] * (1 if policy.startswith('fixed') else exposure_ms/100)
    hits = profiles >= threshold
    valid = hits.any(axis=1)
    last = len(RADII)-1 - np.argmax(hits[:, ::-1], axis=1)
    points = np.column_stack((XX[np.arange(121), last], YY[np.arange(121), last]))[valid]
    rr = RADII[last[valid]]
    return points, dict(threshold=threshold, candidate_count=int(valid.sum()),
                        ray_span_deg=float(np.ptp(ANGLES[valid])*180/np.pi) if valid.any() else 0,
                        ray_longest_run_deg=max(0, longest_run(valid)-1)*1.5,
                        radial_boundary_hits=int(np.sum((rr == RADII[0]) | (rr == RADII[-1])))), use


def fit_circle(points, radius=None, initial=None):
    if len(points) < 10:
        return None
    if initial is None:
        initial = [556., 224., 48.]
    def residual(p):
        return np.linalg.norm(points-p[:2], axis=1) - (p[2] if radius is None else radius)
    fit = least_squares(residual, initial if radius is None else initial[:2],
                        bounds=BOUNDS if radius is None else (BOUNDS[0][:2], BOUNDS[1][:2]),
                        loss='soft_l1', f_scale=1.)
    p = np.r_[fit.x, radius] if radius is not None else fit.x
    delta = points-p[:2]
    distance = np.linalg.norm(delta, axis=1)
    err = distance-p[2]
    jac = np.column_stack((-delta / distance[:, None], -np.ones(len(points))))
    singular = np.linalg.svd(jac, compute_uv=False)
    actual_singular=np.linalg.svd(jac if radius is None else jac[:,:2],compute_uv=False)
    # Unweighted geometric Jacobian diagnostic, not a statistical covariance estimate.
    colcos = float(np.dot(jac[:, 1], jac[:, 2]) /
                   (np.linalg.norm(jac[:, 1])*np.linalg.norm(jac[:, 2])))
    theta = np.arctan2(delta[:, 1], delta[:, 0])
    bins = np.floor((theta+np.pi)/(2*np.pi)*24).astype(int) % 24
    occupied = np.bincount(bins, minlength=24) > 0
    supported = abs(err) <= 2
    support_bins = np.bincount(bins[supported], minlength=24) > 0
    # A cyclic run must not be counted twice when all bins are occupied.
    longest = min(24, longest_run(np.r_[support_bins, support_bins]))*15
    # Consecutive original rays are a stricter continuity check than occupied coarse bins.
    ray_angle=np.arctan2(points[:,1]-ORIGIN[1],points[:,0]-ORIGIN[0])
    ray_angle=np.where(ray_angle>0,ray_angle-2*np.pi,ray_angle)
    ray_ids=np.rint((ray_angle+np.pi)/np.pi*120).astype(int)
    ray_supported=np.zeros(121,dtype=bool)
    ray_supported[ray_ids[supported]]=True
    return p, dict(cx=p[0], cy=p[1], radius=p[2], rms=np.sqrt(np.mean(err**2)),
                   median_abs_residual=np.median(abs(err)), points_within_1px=int(np.sum(abs(err)<=1)),
                   points_within_2px=int(supported.sum()), occupied_15deg_bins=int(occupied.sum()),
                   supported_2px_15deg_bins=int(support_bins.sum()),
                   supported_2px_longest_arc_deg=longest,
                   supported_2px_longest_original_rayrun_deg=max(0,longest_run(ray_supported)-1)*1.5,
                   jac_condition=actual_singular[0]/actual_singular[-1],
                   released_radius_jac_condition=singular[0]/singular[-1], jac_cy_r_cos=colcos,
                   fit_bound_hits=int(np.count_nonzero(fit.active_mask)), cost=fit.cost,
                   optimizer_success=bool(fit.success),
                   optimizer_valid=bool(fit.success and not np.any(fit.active_mask)))


def raw_support(events, points, p):
    xy = np.column_stack((events['x'], events['y'])).astype(float)
    if not len(xy):
        return {}
    tree = cKDTree(xy)
    n1 = np.array(tree.query_ball_point(points, 1., return_length=True))
    n2 = np.array(tree.query_ball_point(points, 2., return_length=True))
    d = np.linalg.norm(xy-p[:2], axis=1)
    theta = np.arctan2(xy[:, 1]-p[1], xy[:, 0]-p[0])
    bins = np.floor((theta+np.pi)/(2*np.pi)*24).astype(int) % 24
    near = abs(d-p[2]) <= 2
    occupied = np.bincount(bins[near], minlength=24) > 0
    outer = (d > p[2]+3) & (d <= p[2]+8)
    # Isolated points: <=1 raw event within 2 pixels. Multiple rays can share one event.
    return dict(candidate_raw_1px_supported_fraction=float(np.mean(n1>0)),
                candidate_raw_2px_event_median=float(np.median(n2)),
                candidate_raw_2px_isolated_fraction=float(np.mean(n2<=1)),
                raw_circle_2px_events=int(near.sum()), raw_circle_2px_15deg_bins=int(occupied.sum()),
                raw_lower_semicircle_2px_bins=int(occupied[12:].sum()),
                raw_outer_3to8px_events=int(outer.sum()),
                candidate_y_ge_240=int(np.sum(points[:, 1]>=240)),
                lower_occlusion_risk=True)


def main():
    manifest = json.loads((DATA/'manifest.json').read_text())
    events = np.concatenate([np.load(DATA/f) for f in manifest['parts']])
    t0 = int(events['t'][0])
    timestamps = events['t'].astype(np.int64)-t0
    rows, sensitivity, perturbations, multi = [], [], [], []
    representative = {}
    for segment, (start, end) in SEGMENTS.items():
        print(segment, flush=True)
        midpoints = np.arange(round(start*100), round(end*100)+1)/100
        for mid in midpoints:
            for duration in [20, 50, 100]:
                tstart, tend = mid-duration/2000, mid+duration/2000
                a, b = np.searchsorted(timestamps, np.rint(np.array([tstart,tend])*1e6).astype(int))
                selected = events[a:b]
                raw_ball = selected[(selected['x']>=475)&(selected['x']<638)&
                                    (selected['y']>=143)&(selected['y']<285)]
                is_rep = abs(mid-(start+end)/2) < .00001
                for policy in POLICIES:
                    points, diag, raw = candidates(selected, duration, policy)
                    for mode, radius in [('free', None), ('fixed48', 48.)]:
                        fitted = fit_circle(points, radius)
                        if fitted is None:
                            continue
                        p, fitdiag = fitted
                        row = dict(segment=segment, midpoint_s=mid, start_s=tstart, end_s=tend,
                                   duration_ms=duration, threshold_policy=policy, mode=mode,
                                   events_total=len(selected), events_sampling_crop=len(raw), **diag,
                                   **fitdiag, **raw_support(raw_ball, points, p))
                        rows.append(row)
                        if is_rep and policy == 'scaled_0.6':
                            for seed_idx, seed in enumerate([[556.,224.,48.], [545.,210.,40.], [570.,238.,62.]]):
                                pp, dd = fit_circle(points, radius, seed)
                                multi.append(dict(segment=segment, duration_ms=duration, mode=mode,
                                                  seed=seed_idx, **dd))
                            for k in range(4):
                                # Four disjoint time blocks; each perturbation deletes one whole block.
                                lo=t0+round((tstart+duration/1000*k/4)*1e6)
                                hi=t0+round((tstart+duration/1000*(k+1)/4)*1e6)
                                reduced=selected[(selected['t']<lo)|(selected['t']>=hi)]
                                pts, _, _ = candidates(reduced, duration*.75, policy)
                                result=fit_circle(pts, radius)
                                if result:
                                    pp, dd=result
                                    perturbations.append(dict(segment=segment,duration_ms=duration,mode=mode,
                                                              perturbation='delete_time_quarter',block=k,
                                                              dx=pp[0]-p[0],dy=pp[1]-p[1],dr=pp[2]-p[2],**dd))
                                # Spatial perturbation drops one 45-degree contiguous original ray sector.
                                theta=np.arctan2(points[:,1]-ORIGIN[1],points[:,0]-ORIGIN[0])
                                theta=np.where(theta>0,theta-2*np.pi,theta)
                                keep=~((theta>=-np.pi+k*np.pi/4)&(theta<-np.pi+(k+1)*np.pi/4+1e-10))
                                result=fit_circle(points[keep], radius)
                                if result:
                                    pp, dd=result
                                    perturbations.append(dict(segment=segment,duration_ms=duration,mode=mode,
                                                              perturbation='delete_45deg_sector',block=k,
                                                              dx=pp[0]-p[0],dy=pp[1]-p[1],dr=pp[2]-p[2],**dd))
                    if is_rep and policy == 'scaled_0.6':
                        for radius in [46.,48.,50.]:
                            pp,dd=fit_circle(points,radius)
                            sensitivity.append(dict(segment=segment,duration_ms=duration,reference_radius=radius,**dd))
                        representative[segment,duration]=(selected.copy(),points,
                            fit_circle(points)[0],fit_circle(points,48)[0])
    write_csv('observations.csv', rows)
    write_csv('radius_sensitivity.csv', sensitivity)
    write_csv('block_perturbations.csv', perturbations)
    write_csv('multistart.csv', multi)
    paired=[]
    for segment in SEGMENTS:
        for mode in ['free','fixed48']:
            for midpoint in sorted({r['midpoint_s'] for r in rows if r['segment']==segment}):
                rr={r['duration_ms']:r for r in rows if r['segment']==segment and r['mode']==mode and
                    r['midpoint_s']==midpoint and r['threshold_policy']=='scaled_0.6'}
                for short,long in [(20,50),(20,100),(50,100)]:
                    if short in rr and long in rr:
                        paired.append(dict(segment=segment,mode=mode,midpoint_s=midpoint,
                                           short_ms=short,long_ms=long,
                                           dx_short_minus_long=rr[short]['cx']-rr[long]['cx'],
                                           dy_short_minus_long=rr[short]['cy']-rr[long]['cy'],
                                           dr_short_minus_long=rr[short]['radius']-rr[long]['radius']))
    write_csv('paired_windows.csv',paired)
    metadata=dict(source=str(DATA),manifest_sha256=hashlib.sha256((DATA/'manifest.json').read_bytes()).hexdigest(),
                  source_t0_us=t0,events=len(events),duration_s=timestamps[-1]/1e6,
                  input_parts_sha256={f:hashlib.sha256((DATA/f).read_bytes()).hexdigest() for f in manifest['parts']},
                  segments=SEGMENTS,midpoint_step_ms=10,windows_ms=[20,50,100],
                  sampling_origin=ORIGIN.tolist(),radii_range=[30,75.5,.5],angles_deg=[-180,0,1.5],
                  thresholds={'scaled':'base threshold * duration_ms / 100','fixed':'0.6 count after sigma=1 blur'},
                  jacobian_note='jac_condition uses actual 3/free or 2/fixed columns; released_radius_jac_condition and jac_cy_r_cos always describe hypothetical free-radius geometry, not covariance.',
                  representative_only='radius46/50, three starts, time-quarter and spatial-sector deletion',
                  bounds=BOUNDS,reference_radius_note='48px is prior upper-arc estimate, not truth',
                  temporal_origin='first CD event, not SDK index start; windows half-open',
                  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  matplotlib=matplotlib.__version__,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    metadata['protocol_hypotheses']=[
        'Equal nominal midpoints isolate temporal accumulation effects from the old start-label shift.',
        'Threshold scaled with exposure compares approximately equal average count density; fixed count is a sensitivity control.',
        'Upper-arc support can come from texture, occlusion or noise; residual and convergence cannot establish silhouette identity.',
        'Fixed 48px is an explicit prior; radius46/50 tests sensitivity, not accuracy.',
        'Time-block/sector deletion and Jacobian diagnose stability and coupling, not confidence intervals.',
        'No independent truth; no absolute accuracy or RPM claim.']
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    plot(rows,representative)
    print(f'Wrote {len(rows)} observations, {len(perturbations)} perturbations',flush=True)


def plot(rows, representative):
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,segment in zip(axes,['early_turn','middle','late']):
        selected,_,_,_=representative[segment,20]
        crop=selected[(selected['x']>=490)&(selected['x']<=630)&(selected['y']>=155)&(selected['y']<=275)]
        ax.scatter(crop['x'],crop['y'],s=2,c=np.where(crop['p']>0,'#cc5522','#3377aa'),alpha=.28,rasterized=True)
        ax.set(xlim=(490,630),ylim=(275,155),aspect='equal',title=f'{segment}, 20ms, raw events only')
    fig.tight_layout();fig.savefig(OUT/'raw_events_only.png',dpi=180);plt.close(fig)
    fig, axes=plt.subplots(4,3,figsize=(15,15))
    for i,segment in enumerate(SEGMENTS):
        for j,duration in enumerate([20,50,100]):
            selected,points,free,fixed=representative[segment,duration]
            ax=axes[i,j]
            crop=selected[(selected['x']>=480)&(selected['x']<=635)&(selected['y']>=155)&(selected['y']<=277)]
            # Every raw event is drawn (no random thinning); duplicate pixels overlap.
            ax.scatter(crop['x'],crop['y'],s=.6,c=np.where(crop['p']>0,'#cc5522','#3377aa'),alpha=.18,rasterized=True)
            ax.scatter(points[:,0],points[:,1],s=9,facecolors='none',edgecolors='#e31a1c',linewidths=.6)
            for p,color in [(free,'#00bfc4'),(fixed,'#33aa33')]:
                ax.add_patch(plt.Circle(p[:2],p[2],fill=False,color=color,lw=1.1))
                ax.plot(*p[:2],'+',color=color,ms=8)
            ax.plot(*ORIGIN,'kx',ms=4)
            ax.axhline(240,color='gray',ls=':',lw=.7)
            ax.set(xlim=(480,635),ylim=(277,155),aspect='equal',
                   title=f'{segment} mid={sum(SEGMENTS[segment])/2:.2f}s, {duration}ms\nfree ({free[0]:.2f},{free[1]:.2f}) r={free[2]:.2f}; fixed48 cy={fixed[1]:.2f}')
    fig.suptitle('All raw events: ON orange / OFF blue; outer candidates red; free cyan / fixed48 green; sampling origin black\nSilhouette identity is unverified; gray y=240 marks prior occlusion concern',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(OUT/'raw_candidate_overlays.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(4,3,figsize=(15,12))
    for i,segment in enumerate(SEGMENTS):
        for j,(field,label) in enumerate([('cx','cx px'),('cy','cy px'),('radius','radius px')]):
            ax=axes[i,j]
            for duration,color in [(20,'#3366cc'),(50,'#ee7733'),(100,'#228833')]:
                for mode,style in [('free','-'),('fixed48','--')]:
                    rr=[r for r in rows if r['segment']==segment and r['duration_ms']==duration and r['threshold_policy']=='scaled_0.6' and r['mode']==mode]
                    ax.plot([r['midpoint_s'] for r in rr],[r[field] for r in rr],style,color=color,lw=1,label=f'{duration}ms {mode}')
            ax.set_ylabel(f'{segment}: {label}');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7,ncol=2)
    for ax in axes[-1]:ax.set_xlabel('Nominal midpoint (s); neighboring windows overlap')
    fig.tight_layout();fig.savefig(OUT/'center_radius_traces.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(4,3,figsize=(15,11))
    for i,segment in enumerate(SEGMENTS):
        for j,duration in enumerate([20,50,100]):
            ax=axes[i,j]
            for policy,style in [('scaled_0.3',':'),('scaled_0.6','-'),('scaled_1.0','-.'),('fixed_0.6','--')]:
                rr=[r for r in rows if r['segment']==segment and r['duration_ms']==duration and r['threshold_policy']==policy and r['mode']=='free']
                ax.plot([r['midpoint_s'] for r in rr],[r['cy'] for r in rr],style,label=policy)
            ax.set_title(f'{segment} {duration}ms free cy');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'threshold_cy_sensitivity.png',dpi=150);plt.close(fig)


if __name__ == '__main__':
    main()
