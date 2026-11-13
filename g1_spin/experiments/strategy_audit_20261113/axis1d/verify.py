import json,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[3]
r=json.loads((p/'results.json').read_text());assert len(r)==87 and all(x['valid'] for x in r)
for f,h in json.loads((p/'freeze.json').read_text()).items():assert hashlib.sha256((root/f).read_bytes()).hexdigest()==h
# Every scalar omega must lie on its declared axis and reproduce signed RPM exactly.
for x in r:
 a=np.asarray(x['axis_camera']);w=np.asarray(x['omega_rad_s']);assert np.linalg.norm(np.cross(a,w))<1e-8;assert abs(np.linalg.norm(w)*60/(2*np.pi)-x['spin_rpm'])<1e-8
(p/'verification.json').write_text(json.dumps(dict(rows=len(r),valid=sum(x['valid'] for x in r),frozen_sources_unchanged=True),indent=2)+'\n')
print('verified',len(r))
