"""Ideal material great-circle edge crossings. No estimator imports or feedback."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation

OUT=Path(__file__).resolve().parent
CENTER=np.array([640.,360.]);RADIUS=120.;REFERENCE=.005;DURATION=.010
ANCHORS=np.array([[.6,0],[0,.6],[-.5,.2],[.2,-.5],[.4,.4],[-.3,-.4]])
NORMALS=np.array([[1.,1.],[1.,-.3],[.2,1.],[.7,-1.],[1.,.1],[.4,1.]])
SCENES={'base':([10.,20.,30.],[0.,0.]),'reverse':([-10.,-20.,-30.],[0.,0.]),
        'second':([-20.,30.,10.],[0.,0.]),'zero':([0.,0.,0.],[0.,0.]),
        'translation':([0.,0.,0.],[400.,-300.]),'mixed':([10.,20.,30.],[400.,-300.])}

def sphere(xy):
    xy=np.asarray(xy,float)
    return np.column_stack([xy,-np.sqrt(1-np.sum(xy*xy,axis=1))])

def edges():
    points=sphere(ANCHORS);result=[]
    for i,(point,n2) in enumerate(zip(points,NORMALS)):
        J=np.vstack([np.eye(2),-point[:2]/point[2]])
        n3=J@np.linalg.solve(J.T@J,n2);n3/=np.linalg.norm(n3)
        assert abs(n3@point)<1e-12
        tangent=np.cross(n3,point);tangent/=np.linalg.norm(tangent)
        result.append(dict(edge_id=i,anchor=point,normal=n3,tangent=tangent,image_normal=n2,J=J))
    return result

def generate(name,omega,velocity):
    omega=np.array(omega,float);velocity=np.array(velocity,float);rows=[];root_errors=[];round_errors=[];checks=[]
    scan=np.linspace(0,DURATION,401);rot=Rotation.from_rotvec(-(scan-REFERENCE)[:,None]*omega)
    for edge in edges():
        anchor=edge['anchor'];normal=edge['normal'];tangent=edge['tangent']
        arc=np.cos(np.linspace(-.14,.14,81))[:,None]*anchor+np.sin(np.linspace(-.14,.14,81))[:,None]*tangent
        forward=Rotation.from_rotvec((scan-REFERENCE)[:,None]*omega)
        projected=np.concatenate([forward.apply(q)[:,:2]*RADIUS+CENTER+(scan-REFERENCE)[:,None]*velocity for q in arc])
        low=np.floor(projected.min(axis=0)-2).astype(int);high=np.ceil(projected.max(axis=0)+2).astype(int)
        count=0
        for x in range(low[0],high[0]+1):
            for y in range(low[1],high[1]+1):
                pixel=np.array([x,y],float)
                xy=(pixel-CENTER-(scan-REFERENCE)[:,None]*velocity)/RADIUS
                valid=np.sum(xy*xy,axis=1)<1
                if not valid.all():continue
                material=rot.apply(sphere(xy));signed=material@normal
                intervals=np.flatnonzero(signed[:-1]*signed[1:]<0)
                # Handle an exact scan root by requiring opposite nonzero signs around it.
                exact=np.flatnonzero((signed[1:-1]==0)&(signed[:-2]*signed[2:]<0))+1
                intervals=np.unique(np.r_[intervals,exact-1])
                def signed_at(t):
                    sxy=(pixel-CENTER-(t-REFERENCE)*velocity)/RADIUS
                    if np.dot(sxy,sxy)>=1:return np.nan
                    point=np.r_[sxy,-np.sqrt(1-np.dot(sxy,sxy))]
                    return float(Rotation.from_rotvec(-omega*(t-REFERENCE)).apply(point)@normal)
                last=None
                for i in intervals:
                    a=scan[i];b=scan[i+1]
                    if signed[i+1]==0:b=scan[min(i+2,len(scan)-1)]
                    t=brentq(signed_at,a,b,xtol=1e-14,rtol=1e-14)
                    if last is not None and abs(t-last)<1e-10:continue
                    last=t
                    sxy=(pixel-CENTER-(t-REFERENCE)*velocity)/RADIUS
                    point=np.r_[sxy,-np.sqrt(1-np.dot(sxy,sxy))]
                    q=Rotation.from_rotvec(-omega*(t-REFERENCE)).apply(point)
                    patch_angle=float(np.arccos(np.clip(q@anchor,-1,1)))
                    if patch_angle>.14+1e-12 or point[2]>=0:continue
                    # Same material bright side for all edges: polarity is sign of crossing derivative.
                    before=signed_at(max(0.,t-1e-7));after=signed_at(min(DURATION,t+1e-7))
                    if before==after:continue
                    polarity=int(after>before);time_us=int(np.rint(t*1e6))+10000
                    if time_us>=20000:continue
                    root_errors.append(abs(signed_at(t)));round_errors.append(abs((time_us-10000)*1e-6-t))
                    rows.append((x,y,time_us,polarity,edge['edge_id'],t,patch_angle));count+=1
        checks.append(dict(edge_id=edge['edge_id'],events=count,pixel_box=[low.tolist(),high.tolist()]))
    rows.sort(key=lambda r:(r[2],r[0],r[1],r[4]))
    x=np.array([r[0] for r in rows],np.int16);y=np.array([r[1] for r in rows],np.int16);times=np.array([r[2] for r in rows],np.int64);pol=np.array([r[3] for r in rows],np.uint8)
    edge_ids=np.array([r[4] for r in rows],np.int8);exact=np.array([r[5] for r in rows]);angles=np.array([r[6] for r in rows])
    if len(rows):
        # Distinct ideal edges may coincide; keep IDs distinct as independent contrast crossings.
        assert np.all(np.diff(times)>=0) and np.isin(pol,[0,1]).all()
    events=OUT/'inputs'/f'{name}.npz';truth=OUT/'truth'/f'{name}.npz'
    np.savez_compressed(events,x=x,y=y,t_us=times,polarity=pol)
    np.savez_compressed(truth,edge_ids=edge_ids,exact_time_s=exact,patch_angle=angles,omega=omega,velocity=velocity,
        anchors=np.array([e['anchor'] for e in edges()]),normals=np.array([e['normal'] for e in edges()]))
    (OUT/'inputs'/f'{name}_geometry.json').write_text(json.dumps(dict(center_px=CENTER.tolist(),radius_px=RADIUS,
        reference_us=15000,window_start_us=10000,window_end_us=20000,center_velocity_px_s=velocity.tolist()),indent=2))
    return dict(scene=name,events=len(rows),root_residual_max=max(root_errors,default=0.),integer_time_rounding_max_s=max(round_errors,default=0.),
        edges=checks,truth_file=str(truth.relative_to(OUT)),events_file=str(events.relative_to(OUT)),
        note='Geometric material edge crossings only; edge IDs and omega are evaluator truth, never runner inputs.')

def main():
    if (OUT/'generation.json').exists():raise FileExistsError('Preserve generated scenes')
    (OUT/'inputs').mkdir(exist_ok=True);(OUT/'truth').mkdir(exist_ok=True)
    config=dict(radius_px=RADIUS,center_px=CENTER.tolist(),duration_s=DURATION,reference_s=REFERENCE,scan_step_us=25,
        patch_angle_rad=.14,anchors=ANCHORS.tolist(),image_normals=NORMALS.tolist(),scenes=SCENES,
        pixel_box='Projected material arc across fixed scan times, padded2px; not chosen using any estimated omega')
    (OUT/'generation_config.json').write_text(json.dumps(config,indent=2))
    results=[]
    for name,(omega,velocity) in SCENES.items():
        result=generate(name,omega,velocity);results.append(result);print(name,result['events'],result['root_residual_max'],flush=True)
    (OUT/'generation.json').write_text(json.dumps(results,indent=2))

if __name__=='__main__':main()
