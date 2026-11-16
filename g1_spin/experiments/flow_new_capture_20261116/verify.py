"""Independent data integrity, normal-system and historical-grid checks."""
import hashlib
import json
from pathlib import Path

import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
CAPTURE=ROOT/'g1_spin/dual_capture/analysis_20261113'
REMOTE_HASHES={17:'b6fd367ae7476a789499770e30b866cffad6f5b34edcbd3367099a3bff6437de',
19:'9a2a442dfca8ffdba874bad67db73155eff5750946af3445e8a37df4fe7016e7',
31:'b6791a2480b23cf1eca8c7444aaacb7379a0e724c68b6bd91f6c3f028f5d71b0',
42:'123449b9ca6829e49a0fe1f502d4f46e86754b6b8963e52f83f96c32e4986d6e'}


def main():
    rows=json.loads((OUT/'results.json').read_text())
    protocol=json.loads((OUT/'protocol.json').read_text())
    assert len(rows)==20
    assert len({(r['time_s'],r['window_ms']) for r in rows})==20
    for k,digest in REMOTE_HASHES.items():
        assert hashlib.sha256((OUT/'inputs'/f't{k}.npz').read_bytes()).hexdigest()==digest
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h for path,h in protocol['sources'].items())
    differences=[]
    for row in rows:
        corrected=row['corrected_coarse']
        matrix=np.array(corrected['matrix']);rhs=np.array(corrected['rhs_px_s'])
        omega=np.linalg.lstsq(matrix,rhs,rcond=None)[0]
        assert np.max(abs(omega-corrected['omega_rad_s']))<1e-10
        for name in ['zero','flow','flow_uncorrected']:
            assert row[name] is None or row[name]['nm_calls']==1
        if row['window_ms'] in [10,50,250]:
            filename={10:'spin_10ms.json',50:'spin_short50ms.json',250:'spin_preliminary.json'}[row['window_ms']]
            old=json.loads((CAPTURE/filename).read_text())
            old=next(r for r in old if r['time']==row['time_s'] and r.get('bound_rpm',1200)==1200)
            difference=abs(row['grid']['outcome']['rpm']-old['free']['spin_rpm'])
            assert difference<.01,(row['time_s'],row['window_ms'],difference)
            differences.append(difference)
    result=dict(remote_input_hashes_match=True,source_hashes_match=True,
                independently_reconstructed_global_constraints=True,rows=20,
                historical_grid_max_rpm_difference=max(differences),historical_grid_checks=len(differences),
                one_NM_per_single_route=True,
                timestamp_boundary_correction='First float-boundary run preserved as results_float_boundaries.json; final run uses integer-us half-open intervals')
    (OUT/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
