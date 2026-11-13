import json,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[3];r=json.loads((p/'direct_results.json').read_text())
assert len(r)==53
for x in r:
 a=np.array(x['direct_axis']['axis_camera']);w=np.array(x['direct_axis']['omega_rad_s']);assert abs(np.linalg.norm(np.cross(a,w)))<1e-8;assert abs(np.linalg.norm(w)*60/(2*np.pi)-x['direct_axis']['spin_rpm'])<1e-8
for f,h in json.loads((p/'direct_freeze.json').read_text()).items():assert hashlib.sha256((root/f).read_bytes()).hexdigest()==h
(p/'direct_verification.json').write_text(json.dumps(dict(rows=len(r),run2_status='unavailable_no_independent_shaft_projection',sources_unchanged=True),indent=2)+'\n')
print('verified',len(r))
