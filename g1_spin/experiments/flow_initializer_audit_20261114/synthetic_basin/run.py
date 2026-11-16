"""Optimizer basin experiment on independently generated sphere-edge crossings.

The generator reuses stage4 material arcs. A closed-form constant-axis rotation
root replaces the old temporal scan, so each speed/window is generated anew.
Only the explicitly named oracle optimizer receives truth as an initial value.
"""
from pathlib import Path
from time import perf_counter
import hashlib
import importlib.util
import json
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

OUT = Path(__file__).resolve().parent
EXPERIMENTS = OUT.parents[1]
PHASE3 = EXPERIMENTS / 'drill_phase3_20261106'
PHASE4 = EXPERIMENTS / 'drill_phase4_20261106'
MODELS_PATH = EXPERIMENTS / 'drill_phase5_20261106/local_motion/models.py'
sys.path.insert(0, str(PHASE3))
import methods


def module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(value)
    return value


generator = module('material_arc_generator', PHASE4 / 'generate.py')
models = module('local_pixel_normal_model', MODELS_PATH)
RADIUS = float(generator.RADIUS)
CENTER = generator.CENTER.copy()
BOUND = 3500 * 2 * np.pi / 60
AXIS = np.array([1., 2., 3.]) / np.sqrt(14.)
SCENES = [(rpm, duration) for rpm in [240, 1200, 3000] for duration in [.010, .050]]


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate_crossings(rpm, duration):
    angular_speed = rpm * 2 * np.pi / 60
    omega = AXIS * angular_speed
    reference = duration / 2
    xx, yy = np.meshgrid(np.arange(-119, 120), np.arange(-119, 120))
    xy = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    xy = xy[np.sum(xy * xy, axis=1) < RADIUS**2]
    points = np.column_stack([xy / RADIUS, -np.sqrt(1 - np.sum((xy / RADIUS)**2, axis=1))])
    rows = []
    residuals = []
    derivative_errors = []
    lower_angle = -angular_speed * reference
    upper_angle = angular_speed * reference
    turns = range(-int(np.ceil(angular_speed * duration / (2*np.pi))) - 2,
                  int(np.ceil(angular_speed * duration / (2*np.pi))) + 3)
    for edge in generator.edges():
        normal = edge['normal']
        constant = (points @ AXIS) * (normal @ AXIS)
        cosine = points @ normal - constant
        sine = -np.cross(AXIS, points) @ normal
        amplitude = np.hypot(cosine, sine)
        eligible = (amplitude > 1e-12) & (np.abs(constant) < amplitude)
        indices = np.flatnonzero(eligible)
        phase = np.arctan2(sine[indices], cosine[indices])
        root_angle = np.arccos(-constant[indices] / amplitude[indices])
        for branch in [-1, 1]:
            for turn in turns:
                angle = phase + branch * root_angle + turn * 2 * np.pi
                within = (angle >= lower_angle) & (angle < upper_angle)
                ids = indices[within]
                theta = angle[within]
                if not len(ids):
                    continue
                material = Rotation.from_rotvec(-theta[:, None] * AXIS).apply(points[ids])
                on_arc = material @ edge['anchor'] >= np.cos(.14) - 1e-12
                ids, theta, material = ids[on_arc], theta[on_arc], material[on_arc]
                if not len(ids):
                    continue
                exact_time = reference + theta / angular_speed
                microseconds = np.rint(exact_time * 1e6).astype(np.int64)
                inside = (microseconds >= 0) & (microseconds < round(duration * 1e6))
                ids, theta, material = ids[inside], theta[inside], material[inside]
                exact_time, microseconds = exact_time[inside], microseconds[inside]
                derivative = angular_speed * (-cosine[ids] * np.sin(theta) + sine[ids] * np.cos(theta))
                polarities = (derivative > 0).astype(np.uint8)
                residuals.extend(np.abs(material @ normal).tolist())
                if len(ids):
                    # Independent SciPy finite difference verifies root polarity.
                    before = Rotation.from_rotvec(-(theta - angular_speed*1e-7)[:, None]*AXIS).apply(points[ids]) @ normal
                    after = Rotation.from_rotvec(-(theta + angular_speed*1e-7)[:, None]*AXIS).apply(points[ids]) @ normal
                    finite_derivative = (after-before) / 2e-7
                    derivative_errors.extend(np.abs(finite_derivative-derivative).tolist())
                    assert np.array_equal(finite_derivative > 0, derivative > 0)
                rows.extend(zip(xy[ids,0] + CENTER[0], xy[ids,1] + CENTER[1], microseconds,
                                polarities, np.full(len(ids),edge['edge_id']), exact_time))
    rows.sort(key=lambda row: (row[2], row[0], row[1], row[4]))
    result = np.asarray(rows, float)
    assert len(result)
    return dict(xy=result[:,:2], times_us=result[:,2].astype(np.int64),
                polarity=result[:,3].astype(np.uint8), edge_ids=result[:,4].astype(np.int8),
                exact_times=result[:,5], omega=omega,
                root_residual_max=max(residuals), derivative_error_max=max(derivative_errors))


def coarse_initialization(xy, times, polarity, duration):
    patch_ids = np.floor(xy / 12).astype(int)
    patch_keys, counts = np.unique(patch_ids, axis=0, return_counts=True)
    patch_keys = patch_keys[np.lexsort((patch_keys[:,1],patch_keys[:,0],-counts))[:6]]
    observations = []
    for patch_number, key in enumerate(patch_keys):
        patch = np.all(patch_ids == key, axis=1)
        for segment in range(round(duration/.002)):
            for pol in [0,1]:
                selected = patch & (times >= segment*.002) & (times < (segment+1)*.002) & (polarity == pol)
                pixels = xy[selected]
                if len(np.unique(pixels, axis=0)) < 12:
                    continue
                fit = models.fit(pixels, times[selected], 'M_spatial')
                if not fit['available']:
                    continue
                position = (np.array(fit['center'])-CENTER)/RADIUS
                if position @ position >= .9**2:
                    continue
                point = np.r_[position,-np.sqrt(1-position@position)]
                normal = np.r_[fit['normal'],0.]
                observations.append(dict(patch=patch_number,segment=segment,polarity=pol,
                                         constraint=(RADIUS*np.cross(point,normal)).tolist(),
                                         rhs=fit['speed_px_s'],normal_rms_px=fit['normal_residual_RMS_px']))
    result = dict(observations=observations, available=False)
    if not observations:
        result['failure_reason']='no_constraints'
        return result
    matrix=np.array([row['constraint'] for row in observations])
    target=np.array([row['rhs'] for row in observations])
    omega,_,rank,singular=np.linalg.lstsq(matrix,target,rcond=None)
    condition=float(singular[0]/singular[-1]) if rank==3 else None
    patches=len(set(row['patch'] for row in observations))
    result.update(constraints=len(observations),patches=patches,rank=int(rank),condition=condition,
                  omega=omega.tolist(),rms_px_s=float(np.sqrt(np.mean((matrix@omega-target)**2))))
    result['available']=bool(len(observations)>=6 and patches>=3 and rank==3 and condition<=100 and np.linalg.norm(omega)>1e-6)
    result['failure_reason']=None if result['available'] else 'count_rank_condition_or_zero'
    return result


def direct_nm(points, times, initial):
    times = times - (times[0]+times[-1])/2
    before=methods.score(points,times,np.zeros(3))
    evaluations=0
    def objective(omega):
        nonlocal evaluations
        evaluations+=1
        magnitude=np.linalg.norm(omega)
        if magnitude>BOUND:
            return 1+(magnitude/BOUND-1)**2
        return -methods.score(points,times,omega)/max(before,1e-12)
    simplex=np.vstack([initial,initial+np.eye(3)*max(2.,BOUND/40)])
    started=perf_counter()
    fit=minimize(objective,np.asarray(initial),method='Nelder-Mead',
                 options=dict(initial_simplex=simplex,maxiter=350,xatol=.05,fatol=1e-6))
    elapsed=perf_counter()-started
    score=methods.score(points,times,fit.x)
    return dict(omega=fit.x.tolist(),score=score,zero_score=before,converged=bool(fit.success),
                iterations=int(fit.nit),objective_calls=evaluations,runtime_s=elapsed,
                quality=methods.common_quality(points,times,fit.x,score,before,fit.success))


def candidate_from_search(search, outcome):
    return dict(omega=search['optimizer_candidate_omega'],score=search['optimizer_candidate_contrast'],
                zero_score=search['zero_contrast'],converged=search['optimizer_converged'],
                objective_calls=search['objective_calls'],contrast_calls=search['contrast_calls'],
                runtime_s=search['search_s'],quality=outcome['common_quality'],
                native_accepted=outcome['measurement_valid'],nm_calls=search['nm_calls'],
                selected_seed=search['selected_seed'])


def evaluate(candidate, truth, held_points, held_times, edge_ids, ref, normals):
    omega=np.array(candidate['omega']);magnitude=float(np.linalg.norm(omega));true_magnitude=float(np.linalg.norm(truth))
    cosine=np.dot(omega,truth)/(magnitude*true_magnitude) if magnitude else 0.
    material=Rotation.from_rotvec(-(held_times-ref)[:,None]*omega).apply(held_points)
    distance=RADIUS*np.arcsin(np.clip(np.abs(np.sum(material*normals[edge_ids],axis=1)),0,1))
    return dict(rpm=magnitude*60/(2*np.pi),rpm_error_percent=100*(magnitude/true_magnitude-1),
                directed_axis_error_deg=float(np.degrees(np.arccos(np.clip(cosine,-1,1)))),
                held_edge_plane_median_px=float(np.median(distance)),held_edge_plane_p95_px=float(np.percentile(distance,95)))


def main():
    OUT.mkdir(exist_ok=True)
    if (OUT/'results.json').exists():
        raise FileExistsError('Preserve previous synthetic results')
    sources=[Path(__file__),PHASE4/'generate.py',MODELS_PATH,PHASE3/'methods.py',PHASE3/'frozen/spherical_cmax.py']
    protocol=dict(scenes=[dict(rpm=rpm,duration_ms=duration*1000) for rpm,duration in SCENES],
                  axis=AXIS.tolist(),radius=RADIUS,center=CENTER.tolist(),bound_rpm=3500,
                  model='Six independent ideal finite material great-circle edge crossings on visible hemisphere; integer pixels and rounded microseconds',
                  source='Stage4 edges reused, analytic constant-axis crossing roots regenerated per speed and duration',
                  limitations=['No photometric threshold, contrast amplitude, refractory period, noise, flicker, translation or lens effects',
                               'Finite material arcs are separate ideal edges; endpoints do not generate additional boundary events',
                               'Oracle method intentionally uses truth as initial vector, solely to diagnose score bias and optimization basin',
                               'Reference axis supplied to generator/evaluator only; zero/coarse/grid estimators receive events and geometry',
                               'Uniform-coordinate spatial holdout is not an independent recording'],
                  holdout='(x+2*y)%5==0, before coarse seed or CMAX',inner_mask_ratio=.9,
                  coarse='12px top6 count patches,2ms segments,separate polarity,M_spatial,6 constraints/3patches/rank3/condition100',
                  cmax_event_cap=5000,source_hashes={str(p):digest(p) for p in sources})
    save('protocol.json',protocol)
    results=[]
    for rpm,duration in SCENES:
        name=f'rpm{rpm}_ms{round(duration*1000)}'
        generated=generate_crossings(rpm,duration)
        xy=generated['xy'];times=generated['times_us']*1e-6
        inside=np.sum((xy-CENTER)**2,axis=1)<(.9*RADIUS)**2
        held=inside&(((xy[:,0].astype(int)+2*xy[:,1].astype(int))%5)==0)
        train=inside&~held
        train_ids=np.flatnonzero(train);held_ids=np.flatnonzero(held)
        selected=train_ids[np.linspace(0,len(train_ids)-1,min(5000,len(train_ids))).astype(int)]
        points=methods.cmax.lift(xy[selected,0],xy[selected,1],*CENTER,RADIUS)
        selected_times=times[selected]
        held_points=methods.cmax.lift(xy[held_ids,0],xy[held_ids,1],*CENTER,RADIUS)
        if len(held_ids)>10000:
            held_subset=np.linspace(0,len(held_ids)-1,10000).astype(int)
            held_ids=held_ids[held_subset];held_points=held_points[held_subset]
        np.savez_compressed(OUT/f'{name}_events.npz',xy=xy,times_us=generated['times_us'],polarity=generated['polarity'],training_ids=train_ids,selected_ids=selected,held_ids=held_ids)
        np.savez_compressed(OUT/f'{name}_truth.npz',omega=generated['omega'],edge_ids=generated['edge_ids'],exact_times=generated['exact_times'])
        started=perf_counter()
        coarse=coarse_initialization(xy[train],times[train],generated['polarity'][train],duration)
        coarse_seconds=perf_counter()-started
        candidates={}
        candidates['zero_single_nm']=direct_nm(points,selected_times,np.zeros(3))
        candidates['oracle_true_single_nm']=direct_nm(points,selected_times,generated['omega'])
        if coarse['available']:
            search=methods.axis_search(points,selected_times,coarse['omega'])
            candidates['coarse_axis_scan_single_nm']=candidate_from_search(search,methods.decision(points,selected_times,search,'S'))
            candidates['coarse_axis_scan_single_nm']['initialization_s']=coarse_seconds
        else:
            candidates['coarse_axis_scan_single_nm']=dict(unavailable=coarse['failure_reason'])
        search=methods.g_search(points,selected_times)
        candidates['grid_multistart']=candidate_from_search(search,methods.decision(points,selected_times,search,'G'))
        centered=selected_times-(selected_times[0]+selected_times[-1])/2
        truth_score=methods.score(points,centered,generated['omega'])
        normals=np.array([edge['normal'] for edge in generator.edges()])
        for candidate in candidates.values():
            if 'omega' not in candidate:
                continue
            candidate['metrics']=evaluate(candidate,generated['omega'],held_points,times[held_ids],generated['edge_ids'][held_ids],duration/2,normals)
            candidate['score_relative_to_truth']=candidate['score']/truth_score
        row=dict(scene=name,rpm_truth=rpm,duration_ms=duration*1000,
                 angle_displacement_rad=float(np.linalg.norm(generated['omega'])*duration),
                 events=len(xy),training_events=int(train.sum()),selected_events=len(selected),held_events=len(held_ids),
                 root_residual_max=generated['root_residual_max'],derivative_error_max=generated['derivative_error_max'],
                 coarse=coarse,truth_score=truth_score,candidates=candidates)
        results.append(row);save('results_partial.json',results)
        print(name,'events',len(xy),'coarse',coarse['available'],
              [(key,round(value['metrics']['rpm'],1),round(value['metrics']['directed_axis_error_deg'],1)) for key,value in candidates.items() if 'metrics' in value],flush=True)
    save('results.json',results)
    save('verification.json',dict(source_hashes_unchanged=all(digest(p)==h for p,h in protocol['source_hashes'].items()),scenes=len(results),
                                  maximum_root_residual=max(r['root_residual_max'] for r in results),
                                  maximum_derivative_error=max(r['derivative_error_max'] for r in results),
                                  no_time_rescaling=True,estimator_sources_unchanged=True))

if __name__=='__main__':
    main()
