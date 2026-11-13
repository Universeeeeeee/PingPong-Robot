"""Fixed confirmation pass after frozen no-winner decision; no parameter changes."""
import json
import numpy as np
from helpers import OUT,methods,digest,write_json,table,hash_arrays
from pipeline import process_case,split
from run import bundle_real,save_detail,compact

def main():
    selection=json.loads((OUT/'selection_development.json').read_text());assert selection['selected_frontend'] is None
    assert json.loads((OUT/'split_diagnostic/verification.json').read_text())['status']=='passed'
    for rel,h in json.loads((OUT/'parameter_freeze.json').read_text())['files'].items():assert digest(OUT/rel)==h,rel
    directory=OUT/'confirmation';directory.mkdir(exist_ok=True)
    if (directory/'summary.json').exists():raise FileExistsError('Preserve confirmation results')
    rows=[];Grefs=[]
    for record in json.loads((OUT/'input_manifest.json').read_text()):
        if record['partition']!='confirmation':continue
        bundle=bundle_real(record);T=split(bundle['xy'],bundle['times'])['T']
        for lane in record['lanes']:
            if not lane['available']:
                for frontend in ['F0','F1','F2']:rows.append(dict(case=record['case'],frontend=frontend,lane=lane['name'],status='geometry_lane_unavailable',reason=lane['failure_reason']))
                continue
            velocity=np.array(lane['velocity']);relative=bundle['xy'][T]-bundle['center']-(bundle['times'][T]-.005)[:,None]*velocity
            points=methods.cmax.lift(relative[:,0],relative[:,1],0,0,bundle['radius']);search=methods.g_search(points,bundle['times'][T]);decision=methods.decision(points,bundle['times'][T],search,'G')
            Grefs.append(dict(case=record['case'],lane=lane['name'],T_ids=bundle['ids'][T].tolist(),input_hash=hash_arrays(bundle['ids'][T],points,bundle['times'][T]),search=search,decision=decision))
            for frontend in ['F0','F1','F2']:
                result,planes,raw=process_case(bundle,lane,frontend,include_LOO=True);name=f'{record["case"]}_{lane["name"]}_{frontend}'
                path=directory/f'{name}.json.gz';save_detail(path,result);np.savez_compressed(directory/f'{name}_planes.npz',**planes);np.savez_compressed(directory/f'{name}_raw_planes.npz',**raw)
                entry=compact(result,path);rows.append(entry);print('confirmation',name,'planes',entry['accepted_planes'],'S',entry['routes']['S_comp']['failure_reason'],flush=True)
    write_json(directory/'summary.json',rows);table(directory/'summary.csv',rows);save_detail(directory/'G_reference.json.gz',Grefs)

if __name__=='__main__':main()
