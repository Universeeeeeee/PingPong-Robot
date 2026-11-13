"""Compare saved recurrence-registration shifts with separately fitted arc centers.

No peak reselection or model tuning. The arc centers remain uncertain measurements,
so discrepancy diagnoses non-independent evidence rather than proving either wrong.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[3]
DISK=ROOT/'g1_spin/experiments/disk_recheck_20261113'
QC=ROOT/'g1_spin/outputs/real_222403/qc/circle_fits.csv'
records=json.loads((DISK/'recurrence.json').read_text())+json.loads((DISK/'high_motion_recurrence/recurrence.json').read_text())
qc=np.genfromtxt(QC,delimiter=',',names=True)
# outline.py explicitly accumulates [relative_s, relative_s+.1); compare at midpoint.
t=qc['relative_s']+.05
rows=[]
for item in records:
    if item['recording']!='222403' or item['background']:continue
    peak=item['highest_peak'];idx=peak['index'];checks=[]
    for source in item['sources']:
        start=source['source_s'];end=start+peak['lag_s']
        assert t[0]<=start<=end<=t[-1]
        center_delta=np.array([np.interp(end,t,qc[k])-np.interp(start,t,qc[k]) for k in ['cx','cy']])
        registration=np.array(source['shifts'][idx],float)
        # Explicit component-wise sum cross-checks norm calculation.
        residual=registration-center_delta
        norm=float(np.linalg.norm(residual));assert np.isclose(norm,np.sqrt(sum(v*v for v in residual)))
        checks.append(dict(source_s=start,target_s=end,registration_shift_px=registration.tolist(),
                           arc_center_delta_px=center_delta.tolist(),discrepancy_px=norm,
                           registered_similarity=source['scores'][idx]))
    rows.append(dict(recording='222403',anchor_s=item['anchor_s'],lag_s=peak['lag_s'],
                     median_discrepancy_px=float(np.median([r['discrepancy_px'] for r in checks])),
                     details=checks))
paths=[QC,DISK/'recurrence.json',DISK/'high_motion_recurrence/recurrence.json',
       ROOT/'g1_spin/outputs/real_222403/qc/outline.py']
result=dict(results=rows,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    interpretation='Registration can absorb phase mismatch. Uncertain 100ms arc centers are not truth; discrepancy does not establish true RPM or center error.')
(OUT/'recurrence_translation_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k!='details'} for r in rows],indent=2))
