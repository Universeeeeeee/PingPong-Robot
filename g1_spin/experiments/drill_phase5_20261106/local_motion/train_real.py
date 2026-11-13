"""T-only model/bootstrap training; cannot read evaluation files or truth."""
import hashlib
import json
from pathlib import Path
import numpy as np
from models import representation,fit

OUT=Path(__file__).resolve().parent
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False))

def bootstrap(xy,times,kind,identity,baseline):
    if not baseline['available']:return dict(draws=32,valid_draws=0,samples=[],normal_spread95_deg=None,speed_deviation95_px_s=None,speed_CI95=None)
    pixels,inverse=np.unique(xy,axis=0,return_inverse=True);seed=int.from_bytes(hashlib.sha256(identity.encode()).digest()[:8],'little');rng=np.random.default_rng(seed)
    rows=[];angles=[];speeds=[];normal=np.array(baseline['normal']);speed=baseline['speed_px_s']
    for draw in range(32):
        sampled=rng.integers(0,len(pixels),len(pixels));indices=np.concatenate([np.flatnonzero(inverse==i) for i in sampled]);model=fit(xy[indices],times[indices],kind)
        row=dict(draw=draw,drawn_pixel_indices=sampled.tolist(),available=model['available'],failure_reason=model['failure_reason'],aligned_speed_px_s=None,normal_angle_deg=None)
        if model['available']:
            n=np.array(model['normal']);v=model['speed_px_s']
            if n@normal<0:n=-n;v=-v
            angle=float(np.degrees(np.arccos(np.clip(n@normal,-1,1))));angles.append(angle);speeds.append(v);row.update(aligned_speed_px_s=v,normal_angle_deg=angle)
        rows.append(row)
    return dict(draws=32,seed=seed,valid_draws=len(speeds),samples=rows,
        normal_spread95_deg=float(np.percentile(angles,95)) if angles else None,
        speed_deviation95_px_s=float(np.percentile(abs(np.array(speeds)-speed),95)) if speeds else None,
        speed_CI95=np.percentile(speeds,[2.5,97.5]).tolist() if speeds else None,
        note='Cluster draws fixed by case/patch/segment/polarity/representation, identical for both model kinds; no H/V or residual resampling.')

def main():
    target=OUT/'real';path=target/'development_training.json'
    if path.exists():raise FileExistsError('Preserve T-only training output')
    config=json.loads((OUT/'pre_real_freeze.json').read_text())
    for rel,h in config['files'].items():assert digest(OUT/rel)==h,rel
    manifest=json.loads((target/'development_inputs.json').read_text());rows=[]
    for case in manifest:
        for p in case['patches']:
            input_path=OUT/p['T_file'];assert digest(input_path)==p['T_sha256'];d=np.load(input_path);xy=d['xy'];times=d['time_s'];pol=d['polarity'];ids=d['ids']
            assert np.all(times<.006) and np.all((xy[:,0]+2*xy[:,1])%5!=0)
            for segment,(low,high) in enumerate(config['segments_s']):
                for polarity in config['polarities']:
                    mask=(times>=low-1e-12)&(times<high-1e-12)&(pol==polarity);training=np.flatnonzero(mask)
                    for rep in config['representations']:
                        chosen=training[representation(xy[training],times[training],ids[training],rep)]
                        identity=f'{case["case"]}/patch{p["patch"]}/segment{segment}/pol{polarity}/{rep}'
                        for kind in config['models']:
                            model=fit(xy[chosen],times[chosen],kind);boot=bootstrap(xy[chosen],times[chosen],kind,identity,model)
                            initial_support=model['available'] and model['unique_pixels']>=12
                            rows.append(dict(case=case['case'],patch=p['patch'],box=p['box'],segment=segment,polarity=polarity,representation=rep,model_kind=kind,
                                start_s=low,end_s=high,T_file=p['T_file'],T_sha256=p['T_sha256'],all_segment_T_ids=ids[training].tolist(),support_ids=ids[chosen].tolist(),
                                model=model,bootstrap=boot,T_support_min12=bool(initial_support),
                                bootstrap_direction_stable=None if boot['normal_spread95_deg'] is None else boot['normal_spread95_deg']<=15,
                                qualification='Pending final supplemental classification; raw T models and all bootstrap draws frozen first'))
        print('T-only saved in memory',case['case'],'patches',len(case['patches']),flush=True)
    write(path,rows);write(target/'development_training_freeze.json',dict(stage='All T-only fits/bootstrap before H/future evaluation',records=len(rows),
        files={'pre_real_freeze.json':digest(OUT/'pre_real_freeze.json'),'train_real.py':digest(OUT/'train_real.py'),'models.py':digest(OUT/'models.py'),
            'real/development_training.json':digest(path),'real/development_inputs.json':digest(target/'development_inputs.json')},
        real_evaluation_performed=False,trainer_read_scope='Only T NPZ; no evaluation_only files, H/V IDs, residuals or labels opened'))
    print('FROZEN',len(rows),'T-only model/bootstrap records')

if __name__=='__main__':main()
