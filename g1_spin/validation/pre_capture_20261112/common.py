import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent

def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def verify_lock():
    for name,sha in json.loads((ROOT/'freeze.json').read_text())['files'].items():
        if digest(ROOT/name)!=sha:raise ValueError('Frozen code changed: '+name)

def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def smooth(ts,xy):
    result=xy.copy()
    for i in range(1,len(ts)-1):
        if np.all(np.diff(ts[i-1:i+2])<=30000):result[i]=np.median(xy[i-1:i+2],axis=0)
    return result

def integers(values,name):
    a=np.asarray(values)
    if not np.issubdtype(a.dtype,np.integer):raise ValueError(name+' must use integer microseconds')
    if len(a) and (np.any(a<0) or np.any(a>np.iinfo(np.int64).max)):raise ValueError('Invalid '+name)
    return a.astype(np.int64)

def read_geometry(path):
    import csv
    rows=list(csv.DictReader(Path(path).open()))
    times=integers([int(r['t_us']) for r in rows],'geometry times')
    if len(times)>1 and np.any(np.diff(times)<=0):raise ValueError('Geometry times must increase strictly')
    if any(r['valid'] not in ['0','1'] for r in rows):raise ValueError('Geometry valid must be 0 or 1')
    good=[i for i,r in enumerate(rows) if r['valid']=='1']
    xy=np.array([[float(rows[i]['cx_px']),float(rows[i]['cy_px'])] for i in good]).reshape(-1,2)
    if not np.isfinite(xy).all():raise ValueError('Nonfinite geometry')
    return times[good],xy

class Events:
    def __init__(self,path,origin_us):
        path=Path(path)
        if not isinstance(origin_us,(int,np.integer)):raise ValueError('Event origin must be integer microseconds')
        self.origin=int(origin_us)
        if self.origin<0:raise ValueError('Negative timestamp origin')
        self.paths=[path];self.parts=[]
        if path.suffix=='.json':
            manifest=json.loads(path.read_text())
            if manifest.get('timestamp_unit')!='us':raise ValueError('Event timestamps must be us')
            self.paths += [path.parent/p for p in manifest['parts']]
            self.parts=[np.load(p,mmap_mode='r',allow_pickle=False) for p in self.paths[1:]]
        else:
            with np.load(path,allow_pickle=False) as data:
                t=integers(data['t_us'],'events');n=len(t)
                a=np.empty(n,dtype=[('x','f8'),('y','f8'),('t','i8'),('p','u1')])
                for key,field in [('x','x'),('y','y'),('t','t_us'),('p','polarity')]:
                    if len(data[field])!=n:raise ValueError('Mismatched events')
                    if field=='polarity' and not np.isin(data[field],[0,1]).all():raise ValueError('Invalid polarity')
                    a[key]=data[field]
                self.parts=[a]
        last=None
        for p in self.parts:
            if len(p):
                integers(p['t'],'events')
                if np.any(np.diff(p['t'].astype(np.int64))<0):raise ValueError('Unsorted events')
                if last is not None and p['t'][0]<last:raise ValueError('Cross-part timestamp regression')
                last=int(p['t'][-1])
        self.offsets=np.cumsum([0]+[len(p) for p in self.parts])
    def window(self,begin,end):
        chunks=[];ids=[]
        for i,p in enumerate(self.parts):
            lo,hi=np.searchsorted(p['t'],[self.origin+int(begin),self.origin+int(end)])
            if hi>lo:chunks.append(p[lo:hi]);ids.append(np.arange(self.offsets[i]+lo,self.offsets[i]+hi,dtype=np.int64))
        if chunks:
            raw=np.concatenate(chunks);physical=np.concatenate(ids)
        else:
            raw=np.empty(0,dtype=[('x','f8'),('y','f8'),('t','i8'),('p','u1')]);physical=np.empty(0,dtype=np.int64)
        if not np.isfinite(raw['x']).all() or not np.isfinite(raw['y']).all() or not np.isin(raw['p'],[0,1]).all():raise ValueError('Invalid event values')
        return raw,physical
