"""Explicit paths; old artifacts are read-only."""
import importlib.util
import hashlib
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent;PHASE4=OUT.parent;ROOT=OUT.parents[3];P3=PHASE4.parent/'drill_phase3_20261106'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def write_json(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False))
def table(path,rows):
    import csv
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows([{k:r.get(k) for k in keys} for r in rows])

flow=module('real_stage3_flow',PHASE4/'frozen/stage3_flow.py')
# Frozen methods imports phase4 common through a controlled search path, never a real-front-end common.
import sys
if str(PHASE4) not in sys.path:sys.path.insert(0,str(PHASE4))
methods=module('real_stage3_methods',PHASE4/'frozen/stage3_methods.py')

def hash_arrays(*arrays):
    h=hashlib.sha256()
    for value in arrays:
        a=np.ascontiguousarray(value);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()

def load_archive(run):
    names={1:'DSP03F100025_20261105_201058',2:'DSP03F100025_20261105_201123',3:'DSP03F100025_20261105_201143',4:'DSP03F100025_20261105_203829',5:'DSP03F100025_20261105_215030'}
    p=ROOT/'g1_spin/data/archive_20261105/decoded'/names[run];manifest=json.loads((p/'manifest.json').read_text())
    e=np.concatenate([np.load(p/name) for name in manifest['parts']]);return e,int(e['t'][0])
