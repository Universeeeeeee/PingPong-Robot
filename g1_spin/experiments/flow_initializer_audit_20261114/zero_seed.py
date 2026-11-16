"""Single Nelder-Mead from zero, same objective/simplex/quality as flow route."""
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy.optimize import minimize

OUT = Path(__file__).resolve().parent
PHASE3 = OUT.parent / 'drill_phase3_20261106'
sys.path.insert(0, str(PHASE3))
from methods import BOUND, cmax, score, decision


def main():
    original = [json.loads(line) for line in (PHASE3/'results.jsonl').open()]
    rows = []
    for old in [row for row in original if row['method']=='G']:
        row = dict(case=old['window'], condition=old['condition'],
                   old_G_valid=old['measurement_valid'],old_G_rpm=old['measurement_rpm'],
                   old_G_axis=old['measurement_axis'])
        path = PHASE3/'windows'/f"{old['window']}.npz"
        if old['geometry_status']!='available' or not path.exists():
            row['status']='geometry_unavailable'
            rows.append(row)
            continue
        data=np.load(path)
        if len(data['xy'])<100:
            row['status']='insufficient_events'
            rows.append(row)
            continue
        started=perf_counter()
        xy=data['xy']
        points=cmax.lift(xy[:,0],xy[:,1],*data['center'],float(data['radius_px']))
        times=(data['time_us']-data['time_us'][0])*1e-6
        dt=times-(times[0]+times[-1])/2
        before=score(points,dt,np.zeros(3))
        calls=0
        def objective(omega):
            nonlocal calls
            calls+=1
            speed=np.linalg.norm(omega)
            if speed>BOUND:
                return 1+(speed/BOUND-1)**2
            return -score(points,dt,omega)/max(before,1e-12)
        simplex=np.vstack([np.zeros(3),np.eye(3)*max(2.,BOUND/40)])
        fit=minimize(objective,np.zeros(3),method='Nelder-Mead',options=dict(
            initial_simplex=simplex,maxiter=350,xatol=.05,fatol=1e-6))
        after=score(points,dt,fit.x)
        search=dict(optimizer_candidate_omega=fit.x.tolist(),zero_contrast=before,
                    optimizer_candidate_contrast=after,optimizer_converged=bool(fit.success),
                    nm_calls=1,objective_calls=calls)
        outcome=decision(points,times,search,'S')
        row.update(status='accepted' if outcome['measurement_valid'] else outcome['common_quality']['failure_reason'],
                   search=search,final=outcome,compute_s=perf_counter()-started)
        rows.append(row)
    (OUT/'zero_seed_results.json').write_text(json.dumps(dict(
        protocol='One NM from exactly zero, no direction/magnitude scan. Same simplex, bound, stopping conditions and common quality as flow-based single-start. Includes sphere lift in timing.',
        rows=rows),indent=2,allow_nan=False)+'\n')
    accepted=[r for r in rows if r['status']=='accepted']
    print('total',len(rows),'accepted',len(accepted),'oldG',sum(r['old_G_valid'] for r in rows),flush=True)
    print('differences',[(r['case'],r['status'],r['old_G_valid']) for r in rows if (r['status']=='accepted')!=r['old_G_valid']],flush=True)

if __name__=='__main__':
    main()
