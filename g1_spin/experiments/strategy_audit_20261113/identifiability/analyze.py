"""Optimistic full-flow geometry conditioning, never a spin measurement."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(name,data):(OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def lift(pixels,centers,radius):
    xy=(pixels-centers)/radius;q=np.sum(xy*xy,axis=1)
    if np.any(q>1+1e-10):raise ValueError('Outside sphere')
    return np.column_stack([xy,-np.sqrt(np.maximum(0,1-q))])

def jacobian(points):
    n=len(points);x,y,z=points.T
    j=np.zeros((n,2,5));j[:,0,0]=1;j[:,1,1]=1;j[:,0,3]=z;j[:,0,4]=-y;j[:,1,2]=-z;j[:,1,4]=x
    return j.reshape(-1,5)

def svd_stats(j):
    singular=np.linalg.svd(j,compute_uv=False)
    tol=max(j.shape)*np.finfo(float).eps*singular[0]
    rank=int(np.count_nonzero(singular>tol))
    return dict(singular=singular.tolist(),rank=rank,condition=float(singular[0]/singular[-1]) if rank==j.shape[1] else None)

def metric(points):
    n=len(points);j=jacobian(points);jn=j/np.sqrt(n);t=j[:,:2];rot=j[:,2:]
    # Orthogonal projection onto the complement of translational flow.
    q=np.linalg.qr(t,mode='reduced')[0];remaining=rot-q@(q.T@rot)
    explicit=rot.reshape(n,2,3);explicit=explicit-explicit.mean(axis=0,keepdims=True)
    assert np.allclose(remaining,explicit.reshape(-1,3),rtol=0,atol=1e-12)
    cosine=np.diag(t.T@rot[:,[1,0]])/np.linalg.norm(t,axis=0)/np.linalg.norm(rot[:,[1,0]],axis=0)
    expected=np.array([points[:,2].mean(),-points[:,2].mean()])/np.sqrt(np.mean(points[:,2]**2))
    assert np.allclose(cosine,expected,rtol=0,atol=1e-12)
    energy=np.sum(remaining**2,axis=0)/np.sum(rot**2,axis=0)
    analytic=np.var(points[:,2])/np.mean(points[:,2]**2)
    assert np.allclose(energy[:2],analytic,rtol=0,atol=1e-12)
    known=svd_stats(rot/np.sqrt(n));unknown=svd_stats(remaining/np.sqrt(n));joint=svd_stats(jn)
    return dict(points=n,z_mean=float(points[:,2].mean()),z_std=float(points[:,2].std()),rho_max=float(np.max(np.linalg.norm(points[:,:2],axis=1))),
        tx_wy_cos=float(cosine[0]),ty_wx_cos=float(cosine[1]),known_translation=known,joint=joint,translation_projected_rotation=unknown,
        rotation_column_energy_retained=energy.tolist(),minimum_singular_ratio=unknown['singular'][-1]/known['singular'][-1])

def geometry_cells(pixels,points,dt):
    cells=np.floor(pixels/4).astype(np.int64);unique,inverse=np.unique(cells,axis=0,return_inverse=True);indices=[]
    for i in range(len(unique)):
        group=np.flatnonzero(inverse==i);middle=np.median(pixels[group],axis=0);indices.append(group[np.argmin(np.linalg.norm(pixels[group]-middle,axis=1))])
    ix=np.array(indices);return points[ix],dt[ix],ix

def ambiguity(points,dt,radius):
    dt=dt-(dt.min()+dt.max())/2
    # Compare added spin against constant center velocity at all occupied locations.
    records=[]
    for axis in [0,1,2]:
        w=np.eye(3)[axis];rot_v=radius*np.cross(w,points)[:,:2]
        weights=dt*dt;velocity=-np.sum(rot_v*weights[:,None],axis=0)/weights.sum()
        linear=(rot_v+velocity)*dt[:,None]
        exact=radius*(Rotation.from_rotvec(dt[:,None]*w).apply(points)-points)[:,:2]+velocity*dt[:,None]
        records.append(dict(delta_omega_rad_s=w.tolist(),delta_rpm_magnitude=60/(2*np.pi),compensating_center_velocity_px_s=velocity.tolist(),
            linear_pixel_rms=float(np.sqrt(np.mean(np.sum(linear*linear,axis=1)))),exact_pixel_rms=float(np.sqrt(np.mean(np.sum(exact*exact,axis=1)))),
            exact_pixel_max=float(np.max(np.linalg.norm(exact,axis=1)))))
    return records

def verify_jacobian():
    points=np.array([[.2,.1,-np.sqrt(.95)],[-.4,.3,-np.sqrt(.75)],[.7,-.2,-np.sqrt(.47)]])
    exact=jacobian(points);eps=1e-6;cols=[]
    for k in range(5):
        param=np.eye(5)[k]*eps
        positive=Rotation.from_rotvec(param[2:]).apply(points)[:,:2]+param[:2]
        negative=Rotation.from_rotvec(-param[2:]).apply(points)[:,:2]-param[:2]
        cols.append(((positive-negative)/(2*eps)).ravel())
    finite=np.column_stack(cols);error=float(np.max(abs(finite-exact)));assert error<1e-8
    # Constant-z ring: transverse spin exactly confounded with translation.
    a=np.linspace(0,2*np.pi,48,endpoint=False);ring=np.column_stack([.5*np.cos(a),.5*np.sin(a),np.full(48,-np.sqrt(.75))]);r=metric(ring)
    assert r['joint']['rank']==3 and np.max(r['rotation_column_energy_retained'][:2])<1e-25
    return dict(finite_difference_max_abs_error=error,constant_z_ring_joint_rank=r['joint']['rank'],constant_z_ring_rotation_energy=r['rotation_column_energy_retained'])

def main():
    if (OUT/'results.json').exists():raise FileExistsError('Preserve evidence')
    rows=[];files=[OUT/'protocol.md',OUT/'analyze.py',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py'];examples=[]
    for path in sorted((ROOT/'g1_spin/experiments/disk_recheck_20261113/inputs').glob('*.npz')):
        data=np.load(path);raw=data['events'];pixels=np.column_stack([raw['x'],raw['y']]).astype(float);centers=data['linear50'];r=48.;t=(raw['t'].astype(np.int64)-int(data['origin_us']))*1e-6
        points=lift(pixels,centers,r);pc,tc,ix=geometry_cells(pixels,points,t)
        rows.append(dict(dataset='disk',window=path.stem,radius=r,event_weighted=metric(points),cell_weighted=metric(pc),ambiguity=ambiguity(points,t,r),input_sha256=digest(path)))
        files.append(path)
    audit=json.loads((ROOT/'g1_spin/experiments/drill_candidate_validation_20261110/windows.json').read_text())
    for job in audit:
        if job['run']==5:continue
        name=f"run{job['run']}_{job['end_us']:08d}";path=ROOT/f'g1_spin/experiments/drill_candidate_validation_20261110/inputs/{name}_candidate.npz'
        if not path.exists():
            rows.append(dict(dataset='drill',window=name,run=job['run'],available=False));continue
        d=np.load(path);e=d['events'];r=float(d['radius']);pixels=np.column_stack([e['x'],e['y']]).astype(float);t=(e['t'].astype(np.int64)-int(d['origin_us']))*1e-6;points=lift(pixels,d['centers'],r);pc,tc,ix=geometry_cells(pixels,points,t)
        rows.append(dict(dataset='drill',window=name,run=job['run'],available=True,radius=r,event_weighted=metric(points),cell_weighted=metric(pc),ambiguity=ambiguity(points,t,r),input_sha256=digest(path)));files.append(path)
    synthetic=[]
    for rho in [.2,.4,.6,.8,.95]:
        grid=np.arange(-rho,rho+.0001,.01);xx,yy=np.meshgrid(grid,grid);keep=xx*xx+yy*yy<=rho*rho;xy=np.column_stack([xx[keep],yy[keep]]);points=np.column_stack([xy,-np.sqrt(1-np.sum(xy*xy,axis=1))])
        for cropped in [False,True]:
            p=points[points[:,1]<=.333] if cropped else points
            synthetic.append(dict(kind='uniform_projected_disk',rho=rho,occluded_lower_y_above_0333=cropped,metrics=metric(p)))
    checks=verify_jacobian();save('results.json',rows);save('uniform_geometry.json',synthetic);save('verification.json',checks)
    save('sources.json',{str(p.relative_to(ROOT)):digest(p) for p in files})
    # Comparable dimensionless conditions, not sensor-error bounds.
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for dataset,color in [('disk','#188c98'),('drill','#c27b2a')]:
        rr=[r for r in rows if r['dataset']==dataset and 'event_weighted' in r]
        for weighting,marker in [('event_weighted','o'),('cell_weighted','x')]:
            axes[0].scatter([abs(r[weighting]['tx_wy_cos']) for r in rr],[r[weighting]['joint']['condition'] for r in rr],s=18,color=color,marker=marker,label=f'{dataset} {weighting}')
            axes[1].scatter([r[weighting]['known_translation']['condition'] for r in rr],[r[weighting]['translation_projected_rotation']['condition'] for r in rr],s=18,color=color,marker=marker)
        axes[2].scatter([r['event_weighted']['points'] for r in rr],[r['event_weighted']['rotation_column_energy_retained'][0]*100 for r in rr],s=18,color=color,label=dataset)
    axes[0].set(xlabel='|cos(translation x, omega y)|',ylabel='Joint normalized geometry condition',title='Strong translation / transverse-spin coupling');axes[0].legend(fontsize=6)
    axes[1].set(xlabel='Known translation: rotation condition',ylabel='After removing unknown translation: condition',title='Optimistic full 2D flow assumption')
    axes[2].set(xlabel='Number of event locations',ylabel='Omega x/y squared column norm retained (%)',title='Event count does not resolve geometry ambiguity');axes[2].legend()
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(OUT/'conditioning.png',dpi=150);plt.close(fig)
    print(json.dumps(checks))

if __name__=='__main__':main()
