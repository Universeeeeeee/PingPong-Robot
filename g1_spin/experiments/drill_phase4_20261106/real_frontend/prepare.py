"""Freeze complete common-ID support and geometry hypotheses before real estimates."""
import csv
import json
import numpy as np
from helpers import OUT,PHASE4,ROOT,P3,digest,write_json,table,load_archive,hash_arrays

WINDOWS=[('development',1,1495000),('development',2,2095000),('development',5,4095000),('development',5,1525000),
    ('confirmation',1,1425000),('confirmation',3,1495000),('confirmation',4,2495000),('confirmation',5,5095000)]

def interpolate(ts,xy,t):
    right=int(np.searchsorted(ts,t,'left'))
    if right<len(ts) and ts[right]==t:return xy[right]
    if right==0 or right==len(ts):raise ValueError('outside_geometry')
    if ts[right]-ts[right-1]>30000:raise ValueError('geometry_gap_over30ms')
    a=(t-ts[right-1])/(ts[right]-ts[right-1]);return xy[right-1]*(1-a)+xy[right]*a

def main():
    if (OUT/'input_manifest.json').exists():raise FileExistsError('Preserve frozen complete-support inputs')
    (OUT/'inputs').mkdir(exist_ok=True);audit=json.loads((P3/'window_audit.json').read_text());rows=[];protected={}
    for run in [1,2,3,4,5]:
        selected=[w for w in WINDOWS if w[1]==run]
        if not selected:continue
        e,origin=load_archive(run)
        for partition,_,begin in selected:
            name=f'run{run}_{begin:08d}';w=next(w for w in audit if w['window']==name)
            source=P3/w['input_file'];d=np.load(source);physical=d['common_ids'];raw=e[physical]
            xy=np.column_stack([raw['x'],raw['y']]).astype(float);times=(raw['t'].astype(np.int64)-origin-begin)*1e-6
            assert len(physical)==w['common_count'] and np.all((times>=0)&(times<.01))
            path=OUT/'inputs'/f'{name}.npz';np.savez_compressed(path,xy=xy,times=times,polarity=raw['p'],ids=physical,center=d['center'],radius=w['radius_px'],raw_time_us=raw['t'],origin_us=origin,start_us=begin)
            lane=[dict(name='B1',velocity=[0.,0.],available=True,source='stage3 fixed B1',reference_time_s=.005)]
            reason=None;velocity=None;startcenter=None;endcenter=None;geometry_path=None
            if w['condition']=='static_prior':reason='manual_static_prior_has_no_observed_secant'
            else:
                if run in [1,5]:
                    geometry_path=P3.parent/'drill_phase1_20261106'/('rotation_geometry.csv' if run==1 else 'control_geometry.csv')
                    gs=[g for g in csv.DictReader(geometry_path.open()) if g['geometry_valid']=='True']
                else:
                    geometry_path=P3.parent/'drill_phase2_20261106/geometry.json'
                    gs=[g for g in json.loads(geometry_path.read_text()) if g['run']==run and g['geometry_valid']]
                ts=np.array([round(float(g['time_s'])*1e6) for g in gs]);centers=np.array([[float(g['cx']),float(g['cy'])] for g in gs]);ix=np.argsort(ts);ts=ts[ix];centers=centers[ix]
                try:startcenter=interpolate(ts,centers,begin);endcenter=interpolate(ts,centers,begin+10000);velocity=((endcenter-startcenter)/.01).tolist()
                except ValueError as exc:reason=str(exc)
            lane.append(dict(name='secant',available=velocity is not None,velocity=velocity,failure_reason=reason,
                source=None if geometry_path is None else str(geometry_path.relative_to(ROOT)),
                geometry_source_sha256=None if geometry_path is None else digest(geometry_path),
                start_center=None if startcenter is None else startcenter.tolist(),end_center=None if endcenter is None else endcenter.tolist(),reference_time_s=.005))
            row=dict(case=name,run=run,start_us=begin,partition=partition,condition=w['condition'],source='complete_stage3_common_ids',
                source_window_input=str(source.relative_to(ROOT)),source_input_sha256=digest(source),source_mask=w['source_mask'],
                physical_events=len(raw),prior_selected_events=w['event_count'],input_file=str(path.relative_to(OUT)),input_sha256=digest(path),lanes=lane,
                center=d['center'].tolist(),radius_px=w['radius_px'],input_hash=hash_arrays(physical,xy,times,raw['p'],d['center']))
            rows.append(row);protected[str(source.relative_to(ROOT))]=digest(source)
            if geometry_path is not None:protected[str(geometry_path.relative_to(ROOT))]=digest(geometry_path)
    write_json(OUT/'input_manifest.json',rows);write_json(OUT/'input_source_hashes.json',protected);table(OUT/'input_manifest.csv',rows)
    print([(r['case'],r['physical_events'],r['lanes'][1]['velocity']) for r in rows])

if __name__=='__main__':main()
