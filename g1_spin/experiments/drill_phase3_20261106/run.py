"""One fixed pass: 108 inputs times G/A/S; no fallback and no re-selection."""
import json
from time import perf_counter
import numpy as np
from common import OUT,check_freeze,digest,array_hash,write_json,table
from flow import extract,initialize
from methods import cmax,precheck,g_search,axis_search,decision

def empty_measurement():
    return dict(candidate_omega=None,candidate_axis=None,candidate_rpm=None,candidate_contrast=None,
        measurement_valid=False,measurement_omega=None,measurement_axis=None,measurement_rpm=None,
        zero_won=None,common_quality=dict(passed=False,failure_reason=None,gain=None))

def main():
    check_freeze();audit=json.loads((OUT/'window_audit.json').read_text());path=OUT/'results.jsonl'
    if path.exists():raise FileExistsError('Preserve all prior state records')
    (OUT/'planes').mkdir(exist_ok=True);(OUT/'spheres').mkdir(exist_ok=True);rows=[]
    with path.open('w') as out:
        for w in audit:
            data=None;points=None;seconds=None;lift_s=0.;flow_s=0.;planes=None;plane_meta=None;plane_file=None;plane_hash=None;sphere_hash=None
            if w['geometry_status']=='available':
                data=np.load(OUT/w['input_file']);assert digest(OUT/w['input_file'])==w['input_file_sha256']
                start=perf_counter();relative=data['xy']-data['center'];points=cmax.lift(relative[:,0],relative[:,1],0,0,w['radius_px'])
                seconds=(data['time_us']-w['start_us'])*1e-6;lift_s=perf_counter()-start
                sphere_hash=array_hash(points,seconds)
                np.savez_compressed(OUT/'spheres'/f'{w["window"]}.npz',points=points,seconds=seconds)
            for method in ['G','A','S']:
                row={k:w[k] for k in ['window','run','condition','partition','geometry_source','geometry_status','start_us','end_us','radius_px','event_count','input_hash']}
                row.update(method=method,sphere_hash=sphere_hash,plane_file=None,plane_hash=None,flow_extraction_s=0.,lift_s=lift_s,
                    ls_s=0.,search_s=0.,pipeline_s=lift_s,flow_metadata=None,initializer=None,init_available=None,
                    native_status='not_executed',native_result=None,multistart_check='not_executed',search=None,
                    terminal_phase='input',failure_reason=None,hypothesis_only=True,physics_validated=False)
                row.update(empty_measurement())
                if points is None:
                    row['failure_reason']='geometry_unavailable';row['common_quality']['failure_reason']='geometry_unavailable'
                elif method=='G':
                    search=g_search(points,seconds);outcome=decision(points,seconds,search,'G')
                    row.update(outcome,search=search,search_s=search['search_s'],pipeline_s=lift_s+search['search_s'],terminal_phase='quality',
                        native_status='accepted' if search['native']['valid'] else search['native']['failure_reason'],native_result=search['native'],
                        multistart_check=search['multistart_check'])
                    row['failure_reason']=None if outcome['measurement_valid'] else (outcome['common_quality']['failure_reason'] or search['native']['failure_reason'])
                else:
                    if planes is None:
                        start=perf_counter();planes,plane_meta=extract(data['xy'],seconds,data['polarity'],data['center'],w['radius_px']);flow_s=perf_counter()-start
                        plane_file=f'planes/{w["window"]}.npz';np.savez_compressed(OUT/plane_file,**planes);plane_hash=digest(OUT/plane_file)
                        write_json(OUT/'planes'/f'{w["window"]}.json',dict(**plane_meta,extraction_s=flow_s,input_hash=w['input_hash'],sphere_hash=sphere_hash,plane_hash=plane_hash))
                    start=perf_counter();init=initialize(planes['points'],planes['gradients'],w['radius_px'],method);ls_s=perf_counter()-start
                    row.update(plane_file=plane_file,plane_hash=plane_hash,flow_extraction_s=flow_s,ls_s=ls_s,flow_metadata=plane_meta,
                        initializer=init,init_available=init['available'],pipeline_s=lift_s+flow_s+ls_s,
                        native_status=init['native_scalar_status'] if method=='S' else 'literal_text_initializer_'+('available' if init['available'] else 'failed'),
                        multistart_check='not_performed_single_start_not_a_pass')
                    input_failure=precheck(points,seconds)
                    if input_failure or not init['available']:
                        row.update(terminal_phase='input' if input_failure else 'initialization',failure_reason=input_failure or init['failure_reason'])
                        row['common_quality']['failure_reason']=row['failure_reason']
                    else:
                        search=axis_search(points,seconds,init['raw_omega_rad_s']);outcome=decision(points,seconds,search,method)
                        row.update(outcome,search=search,search_s=search['search_s'],pipeline_s=lift_s+flow_s+ls_s+search['search_s'],terminal_phase='quality')
                        row['failure_reason']=outcome['common_quality']['failure_reason']
                if not row['measurement_valid']:
                    assert all(row[k] is None for k in ['measurement_omega','measurement_axis','measurement_rpm'])
                out.write(json.dumps(row,allow_nan=False)+'\n');out.flush();rows.append(row)
            last=rows[-3:];print(w['window'],w['geometry_status'],w['event_count'],
                '/'.join(f'{r["method"]}:{r["failure_reason"] or "accepted"}' for r in last),flush=True)
    assert len(rows)==324
    keep=[k for k in rows[0] if k not in ['search','native_result','initializer','flow_metadata']]
    table(OUT/'results.csv',[{k:r[k] for k in keep} for r in rows])
    check_freeze();print('DONE 108 windows / 324 states; no additional experiments',flush=True)

if __name__=='__main__':main()
