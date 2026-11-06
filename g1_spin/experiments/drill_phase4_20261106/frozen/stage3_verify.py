"""Independent 108-window reconstruction, frozen-native parity, and analytic A/S bridge."""
import csv
import json
import sys
from collections import Counter
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from common import OUT,ROOT,P1,P2,load,digest,array_hash,write_json,check_freeze
from flow import initialize,extract
from methods import cmax,axis_search,decision

def interp(ts,coords,queries):
    out=[]
    for t in queries:
        right=int(np.searchsorted(ts,t,'left'))
        if right<len(ts) and ts[right]==t:out.append(coords[right]);continue
        if right==0 or right==len(ts):raise ValueError('geometry_outside_coverage')
        if ts[right]-ts[right-1]>30000:raise ValueError('geometry_gap_over_30ms')
        a=(int(t)-int(ts[right-1]))/(int(ts[right])-int(ts[right-1]))
        out.append(coords[right-1]*(1-a)+coords[right]*a)
    return np.array(out)

def independent_score(points,times,omega):
    dt=times-(times[0]+times[-1])/2
    warped=Rotation.from_rotvec(-dt[:,None]*np.asarray(omega)).apply(points)
    lat=np.arcsin(np.clip(warped[:,2],-1,1));lon=np.arctan2(warped[:,1],warped[:,0])
    yy=(lat/np.pi+.5)*63;xx=((lon/(2*np.pi)+.5)*128)%128
    iy=np.floor(yy).astype(int);ix=np.floor(xx).astype(int);fy=yy-iy;fx=xx-ix
    hist=np.zeros((64,128))
    for dy in [0,1]:
        for dx in [0,1]:
            weight=(fy if dy else 1-fy)*(fx if dx else 1-fx)
            np.add.at(hist,(np.minimum(iy+dy,63),(ix+dx)%128),weight)
    iwe=gaussian_filter(hist,sigma=1,truncate=2,mode=('reflect','wrap'))
    return float(np.var(iwe/len(points)))

def independent_system(points,gradients,radius,method):
    matrix=[];rhs=[]
    for p,g in zip(points,gradients):
        x,y,z=p;a,b=g
        if method=='S':matrix.append(radius*np.cross(p,[a,b,0.]));rhs.append(1.)
        else:
            u,w=g/np.dot(g,g);v=np.array([u/radius,w/radius,-(x*u+y*w)/(radius*z)])
            projection=np.outer(v,v)/np.dot(v,v)
            cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0.]])
            matrix.extend(projection@(-cross));rhs.extend(projection@v)
    return np.array(matrix).reshape(-1,3),np.array(rhs)

def analytic_search_checks():
    import methods as algorithms
    xy=np.array([[.6,0],[0,.6],[-.5,.2],[.2,-.5],[.4,.4],[-.3,-.4]])
    points=np.column_stack([xy,-np.sqrt(1-np.sum(xy*xy,axis=1))]);truth=np.array([10.,20.,30.])
    n=np.array([[1,1],[1,-.3],[.2,1],[.7,-1],[1,.1],[.4,1.]])
    velocity=np.cross(truth,points);g=n/np.sum(n*velocity[:,:2],axis=1)[:,None]
    plane_points=np.tile(points,(2,1));gradients=np.tile(g,(2,1))
    times=np.repeat(np.linspace(0,.01,100),len(points));texture=np.tile(points,(100,1))
    observed=Rotation.from_rotvec(times[:,None]*truth).apply(texture)
    records=[]
    for method in ['A','S']:
        init=initialize(plane_points,gradients,1.,method);assert init['available']
        matrix,rhs=independent_system(plane_points,gradients,1.,method);omega,_,rank,_=np.linalg.lstsq(matrix,rhs,rcond=None)
        assert rank==3 and np.allclose(omega,init['raw_omega_rad_s'],atol=1e-10)
        original=cmax.estimate;original_nm=algorithms.minimize;nm_calls=0
        def forbidden_G(*args,**kwargs):raise AssertionError('Analytic A/S called G')
        def counted_NM(*args,**kwargs):
            nonlocal nm_calls
            nm_calls+=1;return original_nm(*args,**kwargs)
        cmax.estimate=forbidden_G
        algorithms.minimize=counted_NM
        try:search=axis_search(observed,times,init['raw_omega_rad_s'])
        finally:cmax.estimate=original;algorithms.minimize=original_nm
        assert nm_calls==search['nm_calls']==1 and search['magnitude_scan'][0]['rpm']==0
        assert np.array_equal(search['magnitude_scan'][0]['omega'],np.zeros(3))
        directed=np.array(init['raw_omega_rad_s']);directed/=np.linalg.norm(directed)
        for scan in search['magnitude_scan']:
            vector=np.array(scan['omega']);assert np.linalg.norm(np.cross(vector,directed))<1e-9
            assert np.dot(vector,directed)>=-1e-12
        errors=[]
        for scan in search['magnitude_scan']:
            independent=independent_score(observed,times,scan['omega']);errors.append(abs(scan['contrast']-independent))
            assert np.isclose(scan['contrast'],independent,atol=1e-15,rtol=1e-10)
        independent=independent_score(observed,times,search['optimizer_candidate_omega'])
        assert np.isclose(independent,search['optimizer_candidate_contrast'],atol=1e-15,rtol=1e-10)
        records.append(dict(method=method,initializer=init,independent_LS=omega.tolist(),search=search,
            independent_final_score=independent,max_scan_score_difference=max(errors),zero_included=True,
            directed_positive_magnitude_scan_only=True,exactly_one_NM=True,actual_wrapped_minimize_calls=nm_calls,no_G_called=True,
            decision=decision(observed,times,search,method),scope='Analytic coordinate trajectories, not a DVS simulation or additional real recording window.'))
    dummy=dict(optimizer_candidate_omega=[10,20,30],optimizer_candidate_contrast=.5,zero_contrast=1.,optimizer_converged=True)
    rejected=decision(observed,times,dummy,'A');assert rejected['zero_won'] and not rejected['measurement_valid']
    assert rejected['measurement_omega'] is None and dummy['optimizer_candidate_omega']==[10,20,30]
    write_json(OUT/'analytic_search_verification.json',dict(status='passed',records=records,zero_winner_test=rejected,
        caveat='All real A/S initializers failed; successful-search bridge is verified only analytically, not hidden in 324 real-data states.'))

def extraction_positive_fixture():
    i,j=np.meshgrid(np.arange(-6,7),np.arange(-6,7));xy=np.column_stack([640+i.ravel(),360+j.ravel()]).astype(float)
    times=10000+100*(i.ravel()+6)+150*(j.ravel()+6);order=np.argsort(times,kind='stable');xy=xy[order];times=times[order]
    data,metadata=extract(xy,(times-times[0])*1e-6,np.ones(len(times),np.uint8),np.array([640.,360.]),60.)
    assert metadata['accepted_planes']>=12
    expected=np.array([1e-4,1.5e-4]);error=float(np.max(abs(data['gradients']-expected)))
    assert error<1e-14 and float(np.max(data['rms_s']))<1e-12
    init={m:initialize(data['points'],data['gradients'],60.,m) for m in ['A','S']}
    assert not init['S']['available'] and init['S']['rank']<3
    np.savez_compressed(OUT/'extraction_positive_fixture.npz',input_xy=xy,input_time_us=times,**data)
    write_json(OUT/'extraction_positive_fixture.json',dict(status='passed',events=len(times),expected_gradient_s_per_pixel=expected.tolist(),
        accepted_planes=metadata['accepted_planes'],max_gradient_error=error,max_plane_rms_s=float(np.max(data['rms_s'])),
        extraction_metadata=metadata,initializers=init,
        caveat='Exact time-plane extraction branch only. Single image normal leaves scalar spin rank deficient as expected; this is not DVS data or evidence of real flow coverage.'))

def native_plane_checks(audit,details):
    sys.path.insert(0,str(ROOT/'g1_spin/src'))
    from g1_spin.events import EventBatch
    from g1_spin.optical_flow import estimate_flow
    checks=[];seed_audit=[]
    for name in ['run1_01495000','run2_03425000','run4_07495000']:
        w=next(w for w in audit if w['window']==name);d=np.load(OUT/w['input_file'])
        event=EventBatch(d['xy'][:,0],d['xy'][:,1],d['time_us'],d['polarity'])
        native=estimate_flow(event,*d['center'],w['radius_px'])
        saved=next(r for r in details if r['window']==name and r['method']=='A')
        assert native['failure_reason']=='insufficient_stable_local_planes' and native['quality']['accepted_planes']==saved['flow_metadata']['accepted_planes']==0
        checks.append(dict(window=name,native_result=native,saved_accepted_planes=0))
        if name!='run1_01495000':continue
        # Independently audit every seed's least-squares residual on the same selected events.
        xy=d['xy'];times=(d['time_us']-int(d['time_us'][0]))*1e-6;tree=cKDTree(xy)
        for seed in np.linspace(0,len(xy)-1,min(300,len(xy)),dtype=int):
            ids=np.asarray(tree.query_ball_point(xy[seed],4),int)
            ids=ids[(abs(times[ids]-times[seed])<=.002)&(d['polarity'][ids]==d['polarity'][seed])]
            row=dict(seed=int(seed),neighbors=len(ids),status='too_few_neighbors',spatial_ratio=None,std_s=None,rms_s=None,rms_over_std=None,threshold_s=None)
            if len(ids)>=10:
                positions=xy[ids];position=positions.mean(axis=0);sp=positions-position;singular=np.linalg.svd(sp,compute_uv=False)
                ratio=float(singular[1]/singular[0]) if singular[0]>1e-6 else 0.;row['spatial_ratio']=ratio
                dt=times[ids]-times[seed];std=float(np.std(dt));design=np.column_stack([sp,np.ones(len(ids))]);coef,_,rank,_=np.linalg.lstsq(design,dt,rcond=None)
                rms=float(np.sqrt(np.mean((design@coef-dt)**2)));threshold=max(20e-6,.15*std)
                row.update(std_s=std,rms_s=rms,rms_over_std=rms/std if std else None,threshold_s=threshold,
                    status='spatial_plane_degenerate' if ratio<.08 else ('insufficient_local_time_span' if std<1e-6 else ('plane_residual_or_rank' if rank<3 or rms>threshold else 'passes_plane_filters')))
            seed_audit.append(row)
        assert not any(r['status']=='passes_plane_filters' for r in seed_audit)
    write_json(OUT/'native_plane_verification.json',dict(status='passed',representative_checks=checks,independent_seed_audit=seed_audit,
        seed_failures=dict(Counter(r['status'] for r in seed_audit)),
        interpretation='Zero accepted planes agrees with the unmodified production optical-flow interface on three selected input batches; no thresholds were changed. This does not isolate downsampling from scene/model effects.'))

def main():
    check_freeze();audit=json.loads((OUT/'window_audit.json').read_text());details=[json.loads(l) for l in (OUT/'results.jsonl').read_text().splitlines()]
    assert len(audit)==108 and len(details)==324 and len({(r['window'],r['method']) for r in details})==324
    checks=[];max_xyz=0.;replay=None;objective_checks=0
    for run in [1,2,3,4,5]:
        events,origin=load(run)
        for w in [w for w in audit if w['run']==run]:
            rr=[r for r in details if r['window']==w['window']];assert {r['method'] for r in rr}=={'G','A','S'}
            assert len({r['input_hash'] for r in rr})==1 and len({r['sphere_hash'] for r in rr})==1
            if w['geometry_status']=='unavailable':
                assert all(r['failure_reason']=='geometry_unavailable' and r['search'] is None for r in rr)
                checks.append(dict(window=w['window'],status='preserved_geometry_failure'));continue
            d=np.load(OUT/w['input_file']);assert digest(OUT/w['input_file'])==w['input_file_sha256'];ids=d['ids'];raw=events[ids]
            assert np.array_equal(raw['t'].astype(np.int64)-origin,d['time_us'])
            assert np.array_equal(np.column_stack([raw['x'],raw['y']]).astype(float),d['xy']);assert np.array_equal(raw['p'],d['polarity'])
            assert np.all((raw['t']>=origin+w['start_us'])&(raw['t']<origin+w['end_us']))
            if w['source_file']:
                source=np.load(w['source_file']);assert np.array_equal(source['ids'],ids) and np.array_equal(source['common_ids'],d['common_ids'])
                assert np.allclose(source['B1'],d['center'][None,:],atol=0,rtol=0)
            else:
                gs=[g for g in json.loads((OUT/'geometry.json').read_text()) if g['run']==run and g['block_s']==w['block_s'] and g['geometry_valid']]
                ts=np.array([round(g['time_s']*1e6) for g in gs]);cs=np.array([[g['cx'],g['cy']] for g in gs]);ref=(w['start_us']+w['end_us'])//2
                center=interp(ts,cs,[ref])[0];assert np.allclose(center,d['center'],atol=1e-10)
                whole=np.flatnonzero((events['t']>=origin+w['start_us'])&(events['t']<origin+w['end_us']))
                xy=np.column_stack([events[whole]['x'],events[whole]['y']]).astype(float)
                sel=np.flatnonzero(np.sum((xy-center)**2,axis=1)<=(.8*w['radius_px'])**2)
                assert np.array_equal(whole[sel],d['common_ids'])
                if len(sel)>5000:sel=sel[np.linspace(0,len(sel)-1,5000).astype(int)]
                assert np.array_equal(whole[sel],ids)
            points=np.column_stack([(d['xy'][:,0]-d['center'][0])/w['radius_px'],(d['xy'][:,1]-d['center'][1])/w['radius_px'],np.zeros(len(ids))])
            points[:,2]=-np.sqrt(1-np.sum(points[:,:2]**2,axis=1));seconds=(d['time_us']-w['start_us'])*1e-6
            sphere=np.load(OUT/'spheres'/f'{w["window"]}.npz');err=float(np.max(abs(points-sphere['points']))) if len(ids) else 0.;max_xyz=max(max_xyz,err);assert err<1e-12
            assert np.all(points[:,2]<=0) and np.allclose(np.linalg.norm(points,axis=1),1)
            assert np.array_equal(sphere['seconds'],seconds) and array_hash(sphere['points'],seconds)==rr[0]['sphere_hash']
            a,s=[next(r for r in rr if r['method']==m) for m in ['A','S']]
            assert a['plane_hash']==s['plane_hash'] and digest(OUT/a['plane_file'])==a['plane_hash']
            assert a['flow_extraction_s']==s['flow_extraction_s']
            for r in rr:
                if not r['measurement_valid']:assert all(r[k] is None for k in ['measurement_omega','measurement_axis','measurement_rpm'])
                assert np.isclose(r['pipeline_s'],r['lift_s']+r['flow_extraction_s']+r['ls_s']+r['search_s'])
                if r['method']!='G' and not r['init_available']:assert r['search'] is None and r['terminal_phase'] in ['input','initialization']
            g=next(r for r in rr if r['method']=='G')
            if g['search'] and g['search']['optimizer_candidate_omega'] is not None:
                independent=independent_score(points,seconds,g['search']['optimizer_candidate_omega'])
                assert np.isclose(independent,g['search']['optimizer_candidate_contrast'],atol=1e-15,rtol=1e-9);objective_checks+=1
            if replay is None and len(ids)>=100:
                result=cmax.estimate(sphere['points'],seconds,max_rpm=3500,initial=None);assert result==g['native_result'];replay=w['window']
            checks.append(dict(window=w['window'],status='input_and_three_routes_verified',events=len(ids)))
        del events
    native_plane_checks(audit,details);analytic_search_checks();extraction_positive_fixture()
    protected=json.loads((OUT/'source_hashes.json').read_text())
    for p,h in protected.items():assert digest(ROOT/p)==h,p
    verify=dict(status='passed',planned_windows=108,result_states=324,available_inputs=sum(w['geometry_status']=='available' for w in audit),
        preserved_geometry_failures=sum(w['geometry_status']=='unavailable' for w in audit),max_xyz_reconstruction_error=max_xyz,
        three_method_input_and_sphere_hashes_equal=True,A_S_shared_planes_equal=True,real_A_S_failed_paths_never_called_search_or_G=True,
        unwrapped_G_exact_replay=replay,independent_real_G_objective_checks=objective_checks,
        A_S_LS_and_objective_checks='analytic_search_verification.json (all real plane sets empty; no claim of real A/S search success)',
        native_plane_crosscheck='native_plane_verification.json',protected_files_unchanged=len(protected),checks=checks,
        caveat='Computation/immutability checks, not physical truth. Four recurrence attempts were frozen before CMax execution.')
    write_json(OUT/'verification.json',verify);print(json.dumps({k:v for k,v in verify.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
