"""Fixed, offline M0/M1/M2 geometry-support ablation; see protocol.md."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

OUT = Path(__file__).resolve().parent
PREVIOUS = OUT.parent / 'geometry_observation_20261101'
spec = importlib.util.spec_from_file_location('previous_geometry', PREVIOUS/'run.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
ROOT, DATA, SEGMENTS = old.ROOT, old.DATA, old.SEGMENTS
ORIGIN, ANGLES, RADII = old.ORIGIN, old.ANGLES, old.RADII
METHODS = ['M0', 'M1', 'M2']
MASK_CENTERS = {'early_turn': [570.915304,217.762814], 'middle': [571.644594,219.409940], 'late': [559.410878,220.107632]}


def write_csv(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def digest_directory(path):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(path.iterdir()) if p.is_file()}


def crop(events):
    return events[(events['x'] >= 470) & (events['x'] < 643) &
                  (events['y'] >= 138) & (events['y'] < 242)]


def extract(events, duration, origin=ORIGIN, base=.6):
    use = crop(events)
    counts = np.bincount((use['y'].astype(int)-138)*173 + use['x'].astype(int)-470,
                         minlength=104*173).reshape(104, 173)
    xx = origin[0]+np.cos(ANGLES[:, None])*RADII
    yy = origin[1]+np.sin(ANGLES[:, None])*RADII
    profiles = map_coordinates(gaussian_filter(counts.astype(float), 1), [yy-138, xx-470], order=1)
    hits = profiles >= base*duration/100
    ids = np.flatnonzero(hits.any(axis=1))
    ix = len(RADII)-1-np.argmax(hits[ids, ::-1], axis=1)
    return ids, np.column_stack((xx[ids, ix], yy[ids, ix])), RADII[ix]


def support(events, points, lo, hi, radius=3.):
    """Count distinct spatial pixels globally and independently inside four time blocks."""
    use = crop(events)
    xy = np.column_stack((use['x'], use['y'])).astype(float)
    counts = np.zeros(len(points), dtype=int)
    distinct = counts.copy()
    blocks = counts.copy()
    if len(xy) and len(points):
        neighborhoods = cKDTree(xy).query_ball_point(points, radius)
        block_ids = np.minimum(3, ((use['t'].astype(np.int64)-lo)*4//(hi-lo)).astype(int))
        pixel_ids = use['y'].astype(np.int64)*1280+use['x'].astype(np.int64)
        for i, near in enumerate(neighborhoods):
            counts[i] = len(near)
            distinct[i] = len(np.unique(pixel_ids[near]))
            blocks[i] = sum(len(np.unique(pixel_ids[np.asarray(near, dtype=int)[block_ids[near] == k]])) >= 2 for k in range(4))
    return counts, distinct, blocks


def runs(ids, radii):
    if not len(ids):
        return []
    cuts = np.flatnonzero((np.diff(ids) != 1) | (abs(np.diff(radii)) > 2))+1
    return np.split(np.arange(len(ids)), cuts)


def fit(points, fixed):
    nan = dict(cx=np.nan, cy=np.nan, radius=np.nan, rms=np.nan, jac_cy_r_cos=np.nan,
               jac_condition=np.nan, fit_valid=False, fit_status='insufficient_points', fit_bound_hits=0)
    if len(points) < 10:
        return nan
    def residual(p):
        return np.linalg.norm(points-p[:2], axis=1)-(48 if fixed else p[2])
    result = least_squares(residual, [556., 224.] if fixed else [556., 224., 48.],
                           bounds=(old.BOUNDS[0][:2], old.BOUNDS[1][:2]) if fixed else old.BOUNDS,
                           loss='soft_l1', f_scale=1.)
    p = np.r_[result.x, 48.] if fixed else result.x
    delta = points-p[:2]
    unit = delta/np.linalg.norm(delta, axis=1)[:, None]
    jac = np.column_stack((-unit, -np.ones(len(points))))
    singular = np.linalg.svd(jac[:, :2] if fixed else jac, compute_uv=False)
    return dict(cx=p[0], cy=p[1], radius=p[2], rms=np.sqrt(np.mean(residual(result.x)**2)),
                jac_cy_r_cos=np.dot(jac[:, 1], jac[:, 2])/(np.linalg.norm(jac[:, 1])*np.linalg.norm(jac[:, 2])),
                jac_condition=singular[0]/singular[-1], fit_valid=bool(result.success and not np.any(result.active_mask)),
                fit_status='valid' if result.success and not np.any(result.active_mask) else 'optimizer_or_bounds',
                fit_bound_hits=int(np.count_nonzero(result.active_mask)))


def evaluate(selected, before, after, lo, hi, duration, identity, origin=ORIGIN, base=.6, support_radius=3.):
    ids0, pts0, rr0 = extract(selected, duration, origin, base)
    count, distinct, blocks = support(selected, pts0, lo, hi, support_radius)
    count2, _, _ = support(selected, pts0, lo, hi, 2.)
    _, bd, bt = support(before, pts0, lo-(hi-lo), lo, support_radius)
    _, ad, at = support(after, pts0, hi, hi+(hi-lo), support_radius)
    supported = (distinct >= 3) & (blocks >= 2)
    masks = [np.ones(len(pts0), dtype=bool), supported]
    m2 = np.zeros(len(pts0), dtype=bool)
    indices = np.flatnonzero(supported)
    for run in runs(ids0[supported], rr0[supported]):
        if len(run) >= 11:
            m2[indices[run]] = True
    masks.append(m2)
    rows, stored = [], {}
    for method, keep in zip(METHODS, masks):
        ids, points, rr = ids0[keep], pts0[keep], rr0[keep]
        longest = max((len(run) for run in runs(ids, rr)), default=0)
        # Twelve bins across the original upper semicircle; rightmost ray belongs to final bin.
        bins = len(np.unique(np.minimum(11, ids//10)))
        gate = bool(len(points) >= 30 and bins >= 6 and longest >= 21)
        avg = lambda x: float(np.mean(x[keep])) if keep.any() else np.nan
        med = lambda x: float(np.median(x[keep])) if keep.any() else np.nan
        diag = dict(**identity, method=method, duration_ms=duration, origin_x=origin[0], origin_y=origin[1],
                    support_radius_px=support_radius, base_threshold=base, threshold=base*duration/100, events_total=len(selected),
                    candidates=len(points), m0_candidates=len(pts0), dropped_fraction=1-len(points)/max(1, len(pts0)),
                    occupied_15deg_bins=bins, angular_span_deg=float(np.ptp(ids)*1.5) if len(ids) else 0.,
                    longest_radial_run_deg=max(0, longest-1)*1.5, engineering_gate=gate,
                    isolated_event_fraction=avg(count2 <= 1), isolated_event_fraction_support_neighborhood=avg(count <= 1), median_distinct_pixels=med(distinct),
                    median_active_time_blocks=med(blocks), support_gate_fraction=avg(supported),
                    before_support_fraction=avg((bd >= 3) & (bt >= 2)),
                    after_support_fraction=avg((ad >= 3) & (at >= 2)),
                    selected_origin_radius_median=med(rr0),
                    m0_origin_radius_median=float(np.median(rr0)) if len(rr0) else np.nan,
                    # Filtering never moves an individual point; composition can shift inward/outward.
                    median_radius_change_from_m0=med(rr0)-float(np.median(rr0)) if keep.any() else np.nan,
                    radial_boundary_hits=int(np.sum((rr == RADII[0]) | (rr == RADII[-1]))))
        for mode in ['free', 'fixed48']:
            fitted = fit(points, mode == 'fixed48')
            rows.append(dict(**diag, mode=mode, **fitted, retain_candidate=bool(gate and fitted['fit_valid'])))
        stored[method] = (ids, points, rr)
    return rows, stored


def transform(events, kind, segment):
    if kind == 'outer_erasure_r_le36':
        center = MASK_CENTERS[segment]
        return events[np.hypot(events['x'].astype(float)-center[0], events['y'].astype(float)-center[1]) <= 36]
    if kind == 'background_plus180':
        use = events[(events['x'] >= 650) & (events['x'] < 823)].copy()
        use['x'] -= 180
        return use
    return events


def semantic_checks():
    dtype = np.dtype([('x', 'u2'), ('y', 'u2'), ('p', 'i2'), ('t', 'i8')])
    hot = np.zeros(100, dtype=dtype); hot['x'] = 550; hot['y'] = 180; hot['t'] = np.arange(100)*200
    _, d, b = support(hot, np.array([[550., 180.]]), 0, 20000)
    assert not ((d >= 3) & (b >= 2)).any()
    cluster = np.array([(x, y, 1, 1000) for x in [549,550,551] for y in [179,180,181]], dtype=dtype)
    _, d2, b2 = support(cluster, np.array([[550., 180.]]), 0, 20000)
    assert not ((d2 >= 3) & (b2 >= 2)).any()
    arc = []
    for theta in ANGLES[20:101]:
        x, y = np.rint(ORIGIN+48*np.array([np.cos(theta), np.sin(theta)])).astype(int)
        for time in [2000, 7000]:
            for dx, dy in [(0,0),(1,0),(-1,0),(0,1),(0,-1)]:
                arc.append((x+dx, y+dy, 1, time))
    arc = np.array(arc, dtype=dtype)
    rows, _ = evaluate(arc, arc[:0], arc[:0], 0, 20000, 20, {'segment':'synthetic', 'midpoint_s':.01})
    assert next(r for r in rows if r['method']=='M2')['engineering_gate']
    dense = np.concatenate([arc, arc.copy()])
    dense['t'][len(arc):] += 10000
    dense_rows, _ = evaluate(dense, dense[:0], dense[:0], 0, 20000, 20, {'segment':'synthetic_dense', 'midpoint_s':.01})
    return dict(repeated_hotspot=dict(distinct_pixels=int(d[0]), active_blocks=int(b[0]), passes=False),
                single_time_cluster=dict(distinct_pixels=int(d2[0]), active_blocks=int(b2[0]), passes=False),
                continuous_multiblock_arc_M2_passes=True,
                high_density_four_block_arc={r['method']: {k:r[k] for k in ['candidates','engineering_gate','longest_radial_run_deg']} for r in dense_rows if r['mode']=='free'})


def main():
    previous_hashes = digest_directory(PREVIOUS)
    checks = semantic_checks()
    manifest = json.loads((DATA/'manifest.json').read_text())
    events = np.concatenate([np.load(DATA/f) for f in manifest['parts']])
    t0 = int(events['t'][0]); times = events['t'].astype(np.int64)
    rows, sensitivity, controls, reps = [], [], [], {}
    def window(lo, hi):
        a,b = np.searchsorted(times, [lo, hi]); return events[a:b]
    for segment, (start,end) in SEGMENTS.items():
        print(segment, flush=True)
        for mid in np.arange(round(start*100), round(end*100)+1)/100:
            for duration in [20,50,100]:
                lo,hi = t0+np.rint(np.array([mid-duration/2000,mid+duration/2000])*1e6).astype(np.int64)
                selected, before, after = window(lo,hi), window(2*lo-hi,lo), window(hi,2*hi-lo)
                identity = dict(segment=segment, midpoint_s=mid, start_s=(lo-t0)/1e6, end_s=(hi-t0)/1e6)
                rr, stored = evaluate(selected,before,after,lo,hi,duration,identity)
                rows.extend(rr)
                if abs(mid-(start+end)/2) < 1e-6:
                    reps[segment,duration] = (selected.copy(), stored)
                    rr,_ = evaluate(selected,before,after,lo,hi,duration,dict(**identity,perturbation='support_radius_2'),support_radius=2.)
                    sensitivity.extend(rr)
                    for base in [.3,1.]:
                        rr,_ = evaluate(selected,before,after,lo,hi,duration,dict(**identity,perturbation=f'threshold_{base}'),base=base)
                        sensitivity.extend(rr)
                    for dx,dy in [(-2,0),(2,0),(0,-2),(0,2)]:
                        rr,_ = evaluate(selected,before,after,lo,hi,duration,dict(**identity,perturbation=f'origin_{dx}_{dy}'),origin=ORIGIN+[dx,dy])
                        sensitivity.extend(rr)
                    if segment != 'translation' and duration in [20,100]:
                        for kind in ['outer_erasure_r_le36','background_plus180']:
                            transformed = [transform(x,kind,segment) for x in [selected,before,after]]
                            rr, ss = evaluate(*transformed,lo,hi,duration,dict(**identity,control=kind))
                            for row in rr:
                                row['control_evaluable'] = row['m0_candidates'] >= 30
                            controls.extend(rr)
                            reps[segment,duration,kind] = (transformed[0].copy(),ss)
    write_csv('observations.csv',rows);write_csv('sensitivity.csv',sensitivity);write_csv('controls.csv',controls)
    metadata = dict(source=str(DATA),source_t0_us=t0,events=len(events),previous_files_sha256=previous_hashes,
                    input_parts_sha256={f:hashlib.sha256((DATA/f).read_bytes()).hexdigest() for f in manifest['parts']},
                    manifest_sha256=hashlib.sha256((DATA/'manifest.json').read_bytes()).hexdigest(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    protocol_sha256=hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest(),
                    expected_rows=4*31*3*3*2,semantic_checks=checks)
    assert previous_hashes == digest_directory(PREVIOUS)
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    plot(reps,rows)
    print(f'Wrote observations={len(rows)}, sensitivity={len(sensitivity)}, controls={len(controls)}',flush=True)


def panel(ax, events, stored, title):
    use=crop(events)
    ax.scatter(use['x'],use['y'],s=1,c=np.where(use['p']>0,'#ce6622','#4477aa'),alpha=.24,rasterized=True)
    for method,color,size in [('M0','#777777',14),('M1','#f5a000',9),('M2','#007777',5)]:
        points=stored[method][1]
        ax.scatter(points[:,0],points[:,1],s=size,facecolors='none',edgecolors=color,linewidths=.7,label=method)
    ax.plot(*ORIGIN,'kx',ms=4)
    ax.set(xlim=(480,632),ylim=(243,151),aspect='equal',title=title)


def plot(reps, rows):
    fig,axes=plt.subplots(4,3,figsize=(15,12))
    for i,segment in enumerate(SEGMENTS):
        for j,d in enumerate([20,50,100]):
            events,stored=reps[segment,d]
            panel(axes[i,j],events,stored,f'{segment} mid={sum(SEGMENTS[segment])/2:.2f}s {d}ms')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('All raw events (ON orange/OFF blue), M0 gray/M1 amber/M2 teal; no fitted circle used for selection\nFiltering does not establish silhouette identity',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(OUT/'candidate_overlays.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(3,4,figsize=(16,9))
    for i,segment in enumerate(['early_turn','middle','late']):
        for j,(d,kind) in enumerate([(20,'outer_erasure_r_le36'),(100,'outer_erasure_r_le36'),(20,'background_plus180'),(100,'background_plus180')]):
            events,stored=reps[segment,d,kind]
            panel(axes[i,j],events,stored,f'{segment} {d}ms {kind}')
            if kind=='outer_erasure_r_le36':
                axes[i,j].add_patch(plt.Circle(MASK_CENTERS[segment],36,fill=False,ls=':',color='purple',lw=.8))
    fig.suptitle('Negative diagnostics: real-event artificial outer erasure / background ROI; mask itself can manufacture an edge',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(OUT/'negative_controls.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(4,3,figsize=(15,11))
    for i,segment in enumerate(SEGMENTS):
        for j,field in enumerate(['cx','cy','radius']):
            ax=axes[i,j]
            for d,color in [(20,'#4477aa'),(50,'#ee7733'),(100,'#228833')]:
                for method,ls in [('M0',':'),('M1','--'),('M2','-')]:
                    q=[r for r in rows if r['segment']==segment and r['duration_ms']==d and r['method']==method and r['mode']=='free']
                    ax.plot([r['midpoint_s'] for r in q],[r[field] if r['fit_valid'] else np.nan for r in q],ls,color=color,lw=1,label=f'{d}ms {method}')
            ax.set_title(f'{segment} {field} (numerically valid fits, gate may reject)');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=6,ncol=3)
    fig.tight_layout();fig.savefig(OUT/'fit_traces.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    main()
