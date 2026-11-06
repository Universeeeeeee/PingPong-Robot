"""Inspect archived event recordings before choosing geometry or spin parameters."""
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
DATA=ROOT/'g1_spin/data/archive_20261105/decoded'
RUNS={'rotation':'DSP03F100025_20261105_201058', 'control':'DSP03F100025_20261105_215030'}

def load(key):
    folder=DATA/RUNS[key]
    m=json.loads((folder/'manifest.json').read_text())
    e=np.concatenate([np.load(folder/f,allow_pickle=False) for f in m['parts']])
    return e, int(e['t'][0])

def select(e,t0,start,ms):
    a,b=np.searchsorted(e['t'],[t0+round(start*1e6),t0+round(start*1e6)+round(ms*1000)])
    return e[a:b]

def picture(e,box=(0,0,1280,720)):
    x0,y0,x1,y1=box
    q=e[(e['x']>=x0)&(e['x']<x1)&(e['y']>=y0)&(e['y']<y1)]
    im=np.full((y1-y0,x1-x0,3),18,np.uint8)
    for p,color in [(0,[75,175,255]),(1,[255,175,65])]:
        a=q[q['p']==p]; im[a['y']-y0,a['x']-x0]=color
    return im

def main():
    frozen=OUT/'frozen';frozen.mkdir(exist_ok=True)
    files={}
    for rel in ['g1_spin/src/g1_spin/spherical_cmax.py','g1_spin/src/g1_spin/optical_flow.py',
                'g1_spin/src/g1_spin/events.py','g1_spin/src/g1_spin/io.py',
                'g1_spin/src/g1_spin/cli.py', 'g1_spin/experiments/translation_compensation_20261102/translation.py']:
        p=ROOT/rel;dest=frozen/p.name
        if dest.exists():assert dest.read_bytes()==p.read_bytes()
        else:shutil.copy2(p,dest)
        files[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
    (OUT/'baseline_hashes.json').write_text(json.dumps(files,indent=2))
    for key in RUNS:
        e,t0=load(key);duration=(int(e['t'][-1])-t0)*1e-6
        samples=np.linspace(.1,duration-.1,24)
        if key=='control':samples=np.r_[.1,.5,1.,1.5,2.,2.5,2.8,3.,3.2,3.4,3.6,3.8,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10.,11.,11.5]
        for ms in ([5,50] if key=='control' else [5]):
            fig,axes=plt.subplots(6,4,figsize=(18,16))
            for ax,s in zip(axes.flat,samples):
                q=select(e,t0,s,ms);ax.imshow(picture(q));ax.set_title(f'{s:.3f}s | {ms}ms | {len(q)} ev',fontsize=9);ax.set_axis_off()
            fig.suptitle(f'{key}: raw events before model fitting; relative to first CD event')
            fig.tight_layout();fig.savefig(OUT/f'{key}_survey_{ms}ms.jpg',dpi=120);plt.close(fig)
        counts=np.bincount(((e['t']-t0)//10000).astype(int))
        fig,ax=plt.subplots(figsize=(12,3));ax.plot(np.arange(len(counts))*.01,counts/1e4)
        ax.set(xlabel='Seconds from first CD event',ylabel='Mevents/s',title=key)
        fig.tight_layout();fig.savefig(OUT/f'{key}_rate.png',dpi=140);plt.close(fig)
        print(key,len(e),duration,flush=True)

if __name__=='__main__':main()
