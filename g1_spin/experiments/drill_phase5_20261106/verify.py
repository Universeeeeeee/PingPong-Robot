"""Minimal independent raw-ID, selection and variance-decomposition verification."""
import json
import numpy as np
from diagnose import OUT,ROOT,FRONT,CASES,CONTROL,load,digest,write

def independent_lower(xy,time,pol):
    groups={}
    for pixel,t,p in zip(xy,time,pol):groups.setdefault((int(pixel[0]),int(pixel[1]),int(p)),[]).append(float(t))
    sse=0.;between=0.;mean=float(np.mean(time)) if len(time) else 0.
    for values in groups.values():
        values=np.array(values);m=float(values.mean());sse+=float(np.sum((values-m)**2));between+=len(values)*(m-mean)**2
    return sse/len(time) if len(time) else 0.,between/len(time) if len(time) else 0.

def main():
    selected=json.loads((OUT/'patch_selection.json').read_text());metrics=json.loads((OUT/'patch_metrics.json').read_text());checks=[];maxdiff=0.
    manifest=json.loads((FRONT/'input_manifest.json').read_text())
    for run in [1,2,5]:
        e,origin,_,_=load(run)
        for case in [c for c in CASES+[CONTROL] if c.startswith(f'run{run}_')]:
            source=next(r for r in manifest if r['case']==case);data=np.load(FRONT/source['input_file']);selection=next(r for r in selected if r['case']==case)
            front=data['xy'][data['times']<.006];counts={}
            for x,y in front:counts[(int(x)//12,int(y)//12)]=counts.get((int(x)//12,int(y)//12),0)+1
            expected=[]
            for cell,count in sorted(counts.items(),key=lambda item:(-item[1],item[0][0],item[0][1])):
                if count<20:continue
                center=np.array(cell)*12+6
                if any(np.linalg.norm(center-np.array(r)*12-6)<16 for r in expected):continue
                expected.append(cell)
                if len(expected)==6:break
            assert expected==[tuple(r['cell']) for r in selection['patches']]
            for m in [m for m in metrics if m['case']==case]:
                patch=m['patch'];box=patch['box'];saved=np.load(OUT/m['events_file']);start=origin+source['start_us']
                # Independent full-file boolean selector verifies physical IDs and complete10ms pixel support.
                mask=(e['t']>=start)&(e['t']<start+10000)&(e['x']>=box[0])&(e['x']<box[2])&(e['y']>=box[1])&(e['y']<box[3])
                ids=np.flatnonzero(mask);assert np.array_equal(saved['ids'],ids)
                assert np.array_equal(saved['xy'],np.column_stack([e[ids]['x'],e[ids]['y']])) and np.array_equal(saved['polarity'],e[ids]['p'])
                assert np.allclose(saved['time_s'],(e[ids]['t'].astype(np.int64)-start)*1e-6,atol=0,rtol=0)
                before=saved['time_s']<.006;within,between=independent_lower(saved['xy'][before],saved['time_s'][before],saved['polarity'][before]);record=m['training_only_decomposition']
                diff=abs(within-record['within_variance_s2']);maxdiff=max(maxdiff,diff);assert diff<1e-18
                assert abs(within+between-record['total_variance_s2'])<1e-18
                checks.append(dict(case=case,patch=patch['patch'],physical_events=len(ids),selection_first6ms_only=True,full10ms_raw_IDs_match=True))
        del e
    old=json.loads((OUT/'old_F0_neighbor_lower_bound_summary.json').read_text())
    expected_counts={'run1_01495000':(256,269),'run2_02095000':(173,186),'run5_04095000':(71,95)}
    for r in old:assert (r['single_value_lower_bound_exceeds_old_gate'],r['eligible_neighbors'])==expected_counts[r['case']]
    for rel,h in json.loads((OUT/'source_hashes.json').read_text()).items():assert digest(ROOT/rel)==h,rel
    assert json.loads((OUT/'measurement_calibration.json').read_text())['status']=='passed'
    result=dict(status='passed',completed_date='2026-11-07',patches_verified=len(checks),selection_counts={r['case']:r['selected_count'] for r in selected},
        max_independent_within_variance_difference_s2=maxdiff,all_raw_IDs_pixels_times_polarities_match=True,selection_uses_first6ms_counts_only=True,
        calibration_passed=True,source_hashes_unchanged=len(json.loads((OUT/'source_hashes.json').read_text())),old_F0_reviewer_counts_reproduced=True,
        checks=checks,scope='Raw/math checkpoint only; no CMax, new model threshold or108 rerun. Confirmation windows not analyzed.')
    write(OUT/'verification.json',result);print({k:v for k,v in result.items() if k!='checks'})

if __name__=='__main__':main()
