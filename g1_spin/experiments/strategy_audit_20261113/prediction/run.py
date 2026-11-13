import importlib.util,json,hashlib
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
s=importlib.util.spec_from_file_location('frozen_cmax',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py');cmax=importlib.util.module_from_spec(s);s.loader.exec_module(cmax)
DISK=ROOT/'g1_spin/experiments/disk_recheck_20261113/inputs';DRILL=ROOT/'g1_spin/experiments/drill_phase3_20261106/windows'
FILES=[('disk_motion',DISK/'high_motion_00450000_50ms.npz',425000,50000),('disk_early',DISK/'early_turn_01750000_50ms.npz',1725000,50000),('disk_middle',DISK/'middle_06350000_50ms.npz',6325000,50000),('disk_late',DISK/'late_11350000_50ms.npz',11325000,50000)]+[(f'drill_{run}',DRILL/f'run{run}_{start:08d}.npz',start,10000) for run,start in [(1,1495000),(2,2095000),(3,1495000),(4,895000),(5,4095000)]]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def load_file(path,start):
 d=np.load(path)
 if 'events' in d:
  e=d['events'];time=(e['t'].astype(np.int64)-int(d['origin_us'])-start)*1e-6;xy=np.column_stack([e['x'],e['y']]).astype(float);pol=e['p'];center=d['linear50'];radius=48.
 else:
  time=(d['time_us']-start)*1e-6;xy=d['xy'];pol=d['polarity'];center=np.tile(d['center'],(len(xy),1));radius=float(d['radius_px'])
 return xy,time,pol,center,radius

def synthetic(axis):
 rng=np.random.default_rng(1);land=rng.normal(size=(24,3));land/=np.linalg.norm(land,axis=1,keepdims=True);times=np.arange(420)/8400
 omega=np.array(axis,float);omega=omega/np.linalg.norm(omega)*247*2*np.pi/60
 pts=np.einsum('tij,kj->tki',Rotation.from_rotvec(times[:,None]*omega).as_matrix(),land);visible=pts[:,:,2]<-.12
 t=np.broadcast_to(times[:,None],visible.shape)[visible];p=np.broadcast_to((np.arange(24)%2)[None,:],visible.shape)[visible];xy=pts[:,:,:2][visible]*60+[128,128];cs=np.tile([128,128],(len(xy),1));keep=np.linalg.norm(xy-cs,axis=1)<=48
 return (xy[keep],t[keep],p[keep],cs[keep],60.),omega

def iwe(points):
 h,w=64,128;yy=(np.arcsin(np.clip(points[:,2],-1,1))/np.pi+.5)*(h-1);xx=((np.arctan2(points[:,1],points[:,0])/(2*np.pi)+.5)*w)%w
 iy=np.floor(yy).astype(int);ix=np.floor(xx).astype(int);fy=yy-iy;fx=xx-ix;im=np.zeros((h,w))
 for dy,wy in [(0,1-fy),(1,fy)]:
  for dx,wx in [(0,1-fx),(1,fx)]:np.add.at(im,(np.minimum(iy+dy,h-1),(ix+dx)%w),wy*wx)
 return gaussian_filter(im,1,truncate=2,mode=('reflect','wrap'))

def score(points,t,pol,centers,radius,omega,duration):
 cut=.6*duration;train=t<cut;hold=~train
 wp=Rotation.from_rotvec(-(t-cut)[:,None]*omega).apply(points);a=iwe(wp[train]);b=iwe(wp[hold]);den=np.linalg.norm(a)*np.linalg.norm(b);cos=float(np.sum(a*b)/den) if den else None
 source=train&(t>=cut-min(.005,cut/3));predictions=[]
 for low,high in zip(np.linspace(cut,duration,5)[:-1],np.linspace(cut,duration,5)[1:]):
  target=(t>=low)&(t<high);mid=(low+high)/2
  if source.sum()<10 or target.sum()<10:predictions.append(dict(evaluable=False));continue
  moved=Rotation.from_rotvec((mid-t[source])[:,None]*omega).apply(points[source]);pred=moved[:,:2]*radius+np.array([np.interp(mid,t,centers[:,i]) for i in range(2)]);visible=moved[:,2]<0;truth=points[target,:2]*radius+centers[target]
  forward=[];back=[]
  for polarity in [0,1]:
   sm=pol[source]==polarity;hm=pol[target]==polarity;px=pred[sm];seen=visible[sm];tx=truth[hm]
   f=np.full(len(px),3.);g=np.full(len(tx),3.)
   if seen.any() and len(tx):f[seen]=np.minimum(cKDTree(tx).query(px[seen])[0],3);g=np.minimum(cKDTree(px[seen]).query(tx)[0],3)
   forward.extend(f.tolist());back.extend(g.tolist())
  predictions.append(dict(evaluable=True,forward_mean=float(np.mean(forward)),reverse_mean=float(np.mean(back)),balanced_mean=float((np.mean(forward)+np.mean(back))/2),source_count=int(source.sum()),target_count=int(target.sum())))
 valid=[r for r in predictions if r['evaluable']]
 return dict(iwe_holdout_cosine=cos,future_balanced_mean= float(np.mean([r['balanced_mean'] for r in valid])) if valid else None,evaluable_slices=len(valid),slices=predictions)

def main():
 if (OUT/'results.json').exists():raise FileExistsError('Preserve experiment')
 paths=[p for _,p,_,_ in FILES]+[OUT/'run.py',OUT/'protocol.md',ROOT/'g1_spin/src/g1_spin/spherical_cmax.py'];hashes={str(p.relative_to(ROOT)):digest(p) for p in paths};save('freeze.json',hashes)
 cases=[(name,load_file(path,start),duration/1e6,None) for name,path,start,duration in FILES]
 for name,axis in [('synthetic_y',[0,1,0]),('synthetic_tilt',[.6,.7,.4])]:
  data,truth=synthetic(axis);cases.append((name,data,.05,truth))
 results=[]
 for name,(xy,t,pol,cs,radius),duration,truth in cases:
  p=(xy-cs)/radius;points=np.column_stack([p,-np.sqrt(1-np.sum(p*p,axis=1))]);train=t<duration*.6;assert train.any() and (~train).any()
  fits=[];candidates={};start=perf_counter()
  for kind in ['real','shuffle101','shuffle102','shuffle103','reversed_training']:
   data=points[train].copy()
   if kind.startswith('shuffle'):data=data[np.random.default_rng(int(kind[7:])).permutation(len(data))]
   elif kind=='reversed_training':data=data[::-1]
   r=cmax.estimate(data,t[train],max_rpm=3500);candidate=r['omega_rad_s'] if r['valid'] else r['quality'].get('rejected_candidate_omega_rad_s')
   fits.append(dict(kind=kind,estimate=r));candidates[kind]=np.array(candidate) if candidate is not None else None
  candidates['zero']=np.zeros(3);candidates['opposite']= -candidates['real'] if candidates['real'] is not None else None
  evaluation={key:score(points,t,pol,cs,radius,omega,duration) if omega is not None else None for key,omega in candidates.items() if key!='reversed_training'}
  reverse_error=None
  # Reversing only points at nonuniform event times is NOT exact physical time reversal; report diagnosis only.
  if candidates['real'] is not None and candidates['reversed_training'] is not None:reverse_error=float(np.linalg.norm(candidates['real']+candidates['reversed_training']))
  result=dict(case=name,duration_s=duration,train_events=int(train.sum()),holdout_events=int((~train).sum()),fits=fits,evaluation=evaluation,reversal_sum_omega_norm=reverse_error,elapsed_s=perf_counter()-start,
       geometry='previously frozen full-window hypothesis; conditional not end-to-end causal',truth_omega=truth.tolist() if truth is not None else None)
  results.append(result);save('partial_results.json',results);print(name,'valid',fits[0]['estimate']['valid'],'rpm',fits[0]['estimate']['spin_rpm'],'eval', {k:(round(v['iwe_holdout_cosine'],3),round(v['future_balanced_mean'],3)) if v and v['future_balanced_mean'] is not None else None for k,v in evaluation.items()},flush=True)
 save('results.json',results);assert all(digest(ROOT/p)==h for p,h in hashes.items());save('run_check.json',dict(cases=len(cases),fit_calls=len(cases)*5,frozen_sources_unchanged=True))

if __name__=='__main__':main()
