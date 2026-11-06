"""Read-only archive/stage2 geometry helpers; stage3 output is always explicit."""
import importlib.util
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
P1=OUT.parent/'drill_phase1_20261106';P2=OUT.parent/'drill_phase2_20261106'
spec=importlib.util.spec_from_file_location('stage2_geometry_readonly',P2/'run.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
load=base.load;select=base.select;picture=base.picture;centers=base.centers;fit=base.fit
digest=base.digest;write_json=base.write_json;table=base.table

def array_hash(*arrays):
    import hashlib
    h=hashlib.sha256()
    for a in arrays:
        a=np.ascontiguousarray(a);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()

def protected():
    paths=list(P1.rglob('*'))+list(P2.rglob('*'))+list((ROOT/'g1_spin/src').rglob('*.py'))
    arch=ROOT/'g1_spin/data/archive_20261105'
    paths+=list((arch/'raw').glob('*.raw'))+[arch/'recording_labels.json',arch/'verification.json']
    for name in base.NAMES.values():
        folder=arch/'decoded'/name;m=json.loads((folder/'manifest.json').read_text())
        paths+=[folder/'manifest.json']+[folder/p for p in m['parts']]
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths)) if p.is_file() and '__pycache__' not in str(p)}

def check_freeze():
    for p,h in json.loads((OUT/'freeze.json').read_text())['files'].items():assert digest(OUT/p)==h,p
