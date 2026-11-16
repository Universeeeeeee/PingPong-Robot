"""Equal-search-budget generic directed seeds for the six initial audit cases."""
import json
from pathlib import Path
import sys

import numpy as np

OUT = Path(__file__).resolve().parent
PHASE3 = OUT.parent / 'drill_phase3_20261106'
sys.path.insert(0, str(PHASE3))
from methods import axis_search, decision, cmax


def main():
    cases = ['run1_01495000','run2_02095000','run3_01495000',
             'run4_02495000','run5_04095000','run5_01525000']
    directions = [('x+', [1,0,0]), ('x-', [-1,0,0]), ('y+', [0,1,0]),
                  ('y-', [0,-1,0]), ('z+', [0,0,1]), ('z-', [0,0,-1])]
    random = np.random.default_rng(20261114).normal(size=(4,3))
    directions += [(f'random{i}', vector.tolist()) for i,vector in enumerate(random)]
    records = []
    for case in cases:
        d = np.load(PHASE3/'windows'/f'{case}.npz')
        xy = d['xy']
        if len(xy)<100:
            records.append(dict(case=case, failure_reason='insufficient_events'))
            continue
        p = cmax.lift(xy[:,0],xy[:,1],*d['center'],float(d['radius_px']))
        t = (d['time_us']-d['time_us'][0])*1e-6
        for name, axis in directions:
            search = axis_search(p,t,axis)
            outcome = decision(p,t,search,'S')
            records.append(dict(case=case, seed=name, axis=axis,
                                selected_seed=search['selected_seed'],
                                search_s=search['search_s'],nm_calls=search['nm_calls'],
                                objective_calls=search['objective_calls'], outcome=outcome))
        print(case,flush=True)
    (OUT/'seed_controls.json').write_text(json.dumps(dict(
        protocol='Six coordinate directions and four fixed random directions, evaluated separately with identical directed magnitude scan and one NM; no best-of selection claimed as single-start method',
        rows=records),indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
