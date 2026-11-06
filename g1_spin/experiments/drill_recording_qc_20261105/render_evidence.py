"""Compact raw-event evidence figure from preserved 100 ms samples."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
boxes=[(540,260,740,450),(520,395,730,600),(470,330,680,530),(490,260,700,460)]
labels=['1 | vertical-axis recording','2 | camera-facing recording',
        '3 | oblique-axis recording','4 | reverse, before axis changes']
fig,axes=plt.subplots(2,2,figsize=(9,9))
for ax,folder,box,label in zip(axes.flat,sorted((ROOT/'results').glob('DSP*')),boxes,labels):
    data=np.load(folder/'sample_2s_100ms.npz',allow_pickle=False)
    e=data['events']; t0=int(data['origin_us']); e=e[(e['t']>=t0+2_000_000)&(e['t']<t0+2_002_000)]
    x0,y0,x1,y1=box;e=e[(e['x']>=x0)&(e['x']<x1)&(e['y']>=y0)&(e['y']<y1)]
    im=np.full((y1-y0,x1-x0,3),18,np.uint8)
    for p,color in [(0,[75,175,255]),(1,[255,175,65])]:
        q=e[e['p']==p];im[q['y']-y0,q['x']-x0]=color
    ax.imshow(im,extent=[x0,x1,y1,y0],interpolation='nearest')
    ax.set_title(label,fontsize=10);ax.set_xlabel('x (pixels)');ax.set_ylabel('y (pixels)')
fig.suptitle('Raw event evidence | 2.000-2.002 s in each recording\nAxes are user-provided conditions, not estimated results',fontsize=11)
fig.tight_layout();fig.savefig(ROOT/'raw_evidence.png',dpi=150)
