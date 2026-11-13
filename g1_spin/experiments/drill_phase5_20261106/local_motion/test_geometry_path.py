"""Future path guard only: out-of-sphere T observation invalidates the whole system."""
import json
import numpy as np
from qualify_training import OUT,build_system,write,digest

def main():
    def observation(i,center):
        return dict(observation_id=i,patch=i,T_qualification={'training_qualified':True},
            model=dict(center=center,normal=[1.,0.],speed_px_s=100.))
    rows=[observation(0,[0.,0.]),observation(1,[120.,0.]),observation(2,[10.,5.])]
    system=build_system(rows,np.zeros(2),60.)
    assert system['failure_reason']=='geometry_invalid' and system['spherical_geometry_invalid']==[1]
    assert system['observation_ids']==[0,1,2] and not system['rotation']['available'] and system['leave_patch_out']==[]
    old=json.loads((OUT/'real/development_T_global_systems.json').read_text());assert all(not x['spherical_geometry_invalid'] for x in old)
    result=dict(status='passed',test='one invalid qualified sphere point rejects whole global system without silently skipping row or LOPO mask mismatch',
        development_invalid_systems=0,development_values_rerun=False,old_dev_freezes_unchanged=True,
        old_source='frozen_preconfirmation/qualify_training_dev_v1.py',old_source_sha256=digest(OUT/'frozen_preconfirmation/qualify_training_dev_v1.py'),
        new_source_sha256=digest(OUT/'qualify_training.py'),future_only_change='No dev system had invalid geometry; only invalid-input branch changes for later confirmation.')
    write(OUT/'geometry_path_amendment.json',result);print(result)

if __name__=='__main__':main()
