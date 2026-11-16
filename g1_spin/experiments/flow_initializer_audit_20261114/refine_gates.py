"""Send all eligible relaxed-plane candidates through single-start CMAX."""
import json
from pathlib import Path
import sys

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
PHASE3 = ROOT / 'g1_spin/experiments/drill_phase3_20261106'
sys.path.insert(0, str(PHASE3))
from methods import axis_search, decision, g_search, cmax


def summarize(search):
    return {key: value for key, value in search.items()
            if key not in ['native', 'objective_history', 'magnitude_scan', 'local_solutions']}


def main():
    diagnostics = json.loads((OUT/'gate_diagnostics.json').read_text())['rows']
    chosen = [r for r in diagnostics if r['half_span_ms']==2 and r['residual_ratio_limit'] in [.5,.8]]
    baselines = {}
    outputs = []
    for row in chosen:
        candidate = row['unchecked_candidate']
        gate = bool(row['accepted']>=12 and candidate and candidate['rank']==3 and
                    candidate['condition'] is not None and candidate['condition']<=100)
        record = dict(case=row['case'], representation=row['representation'],
                      residual_ratio_limit=row['residual_ratio_limit'], planes=row['accepted'],
                      production_solver_valid=row['production_solver']['valid'],
                      coarse_candidate=candidate, seed_available=gate,
                      qualification='12 planes, rank 3, condition <=100; final measurement residual gate withheld for coarse-seed diagnostic')
        if gate:
            data = np.load(PHASE3/'windows'/f"{row['case']}.npz")
            xy=data['xy']
            points = cmax.lift(xy[:,0],xy[:,1],*data['center'],float(data['radius_px']))
            times=(data['time_us']-data['time_us'][0])*1e-6
            if row['case'] not in baselines:
                grid = g_search(points,times)
                baselines[row['case']] = dict(search=summarize(grid),outcome=decision(points,times,grid,'G'))
            fit=axis_search(points,times,candidate['omega_rad_s'])
            record.update(search=summarize(fit),outcome=decision(points,times,fit,'S'))
        outputs.append(record)
        (OUT/'gate_refinement.json').write_text(json.dumps(dict(rows=outputs,baselines=baselines),indent=2,allow_nan=False)+'\n')
        print(row['case'],row['representation'],row['residual_ratio_limit'],
              record.get('outcome',{}).get('candidate_rpm'),flush=True)

if __name__ == '__main__':
    main()
