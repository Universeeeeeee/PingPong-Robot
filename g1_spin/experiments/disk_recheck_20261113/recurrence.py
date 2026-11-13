"""Translation-registered image recurrence; independent of spherical omega fitting."""
import json,hashlib
from pathlib import Path
from functools import lru_cache
import numpy as np,cv2
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
JOBS={'222403':[1.75,6.35,11.35],'221756':[2.,14.,24.]}

def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
    if (OUT/'recurrence.json').exists():raise FileExistsError('Preserve recurrence')
    allrows=[];fig,axes=plt.subplots(2,3,figsize=(15,8));row_index=0
    for name,anchors in JOBS.items():
        folder=ROOT/f'g1_spin/data/decoded/local_analysis_{name}';meta=json.loads((folder/'manifest.json').read_text());e=np.concatenate([np.load(folder/p) for p in meta['parts']]);origin=int(e['t'][0]);t=e['t'].astype(np.int64)-origin
        @lru_cache(maxsize=1200)
        def picture(mid_us):
            a,b=np.searchsorted(t,[mid_us-5000,mid_us+5000]);q=e[a:b];image=np.zeros((720,1280,2),np.float32)
            image[q['y'],q['x'],q['p']]=1
            return cv2.GaussianBlur(image,(5,5),1)
        if name=='222403':
            qc=np.genfromtxt(ROOT/'g1_spin/outputs/real_222403/qc/circle_fits.csv',delimiter=',',names=True);valid=np.isfinite(qc['cx'])&np.isfinite(qc['cy']);qt=qc['relative_s'][valid]+.05
        for col,anchor in enumerate(anchors):
            sources=np.array([anchor+x for x in [-.10,-.06,-.02,.02,.06,.10]]);lags=np.arange(.15,1.60001,.01);curves=[];detail=[]
            for background in [False,True]:
                per_source=[]
                for source in sources:
                    cx=np.interp(source,qt,qc['cx'][valid]) if name=='222403' else 552
                    cy=np.interp(source,qt,qc['cy'][valid]) if name=='222403' else 215
                    cx-=180 if background else 0;cy-=10 # up in image to avoid the occluded lower region
                    x,y=round(cx)-20,round(cy)-20
                    template=picture(round(source*1e6))[y:y+40,x:x+40].copy()
                    support=int(np.count_nonzero(template.sum(axis=2)>.1));scores=[];shifts=[]
                    for lag in lags:
                        target=picture(round((source+lag)*1e6))[y-8:y+48,x-25:x+65]
                        if support<20 or np.count_nonzero(target.sum(axis=2)>.1)<20:
                            scores.append(None);shifts.append(None);continue
                        match=cv2.matchTemplate(target,template,cv2.TM_CCORR_NORMED);_,score,_,loc=cv2.minMaxLoc(match)
                        scores.append(float(score));shifts.append([int(loc[0]-25),int(loc[1]-8)])
                    per_source.append(dict(source_s=float(source),template_xy=[x,y],support=support,scores=scores,shifts=shifts))
                matrix=np.array([[np.nan if x is None else x for x in r['scores']] for r in per_source]);count=np.isfinite(matrix).sum(axis=0)
                curve=np.array([float(np.nanmedian(matrix[:,i])) if count[i]>=4 else np.nan for i in range(len(lags))])
                clean=np.nan_to_num(curve,nan=-1);peaks,_=find_peaks(clean,prominence=.03,distance=8)
                candidates=[dict(lag_s=float(lags[k]),score=float(curve[k]),rpm_if_one_turn=float(60/lags[k]),support=int(count[k]),index=int(k)) for k in peaks if count[k]>=4]
                # Preserve all peaks in lag order; top by similarity is diagnostic only.
                item=dict(recording=name,anchor_s=anchor,background=background,lags_s=lags.tolist(),curve=[None if not np.isfinite(x) else float(x) for x in curve],peaks=candidates,
                    highest_peak=max(candidates,key=lambda x:x['score']) if candidates else None,sources=per_source,
                    note='Local repeated image pattern under 2D translation, not independently measured full revolutions. Frequency pertains to source/target interval, not necessarily source instant.')
                allrows.append(item);ax=axes[row_index,col];ax.plot(lags,curve,label='background' if background else 'ball');ax.set(title=f'{name} source~{anchor:g}s',xlabel='Forward lag (s)',ylabel='Translation-registered similarity');ax.legend();ax.grid(alpha=.2)
                if not background and candidates:
                    # Top 3 independent recurrence peaks; no CMAX used to choose.
                    chosen=sorted(candidates,key=lambda x:x['score'],reverse=True)[:3];srcrow=per_source[2];sx,sy=srcrow['template_xy'];st=srcrow['source_s'];tiles=[]
                    for lag in [0]+[c['lag_s'] for c in chosen]:
                        img=picture(round((st+lag)*1e6));crop=img[sy-8:sy+48,sx-25:sx+65]
                        bgr=np.zeros((*crop.shape[:2],3),np.float32);bgr[:,:,0]=crop[:,:,0];bgr[:,:,2]=crop[:,:,1];bgr[:,:,1]=crop.sum(axis=2)*.4
                        bgr=np.clip(bgr*255,0,255).astype(np.uint8);bgr=cv2.resize(bgr,(360,224));cv2.putText(bgr,f't={st+lag:.3f}s lag={lag:.2f}',(8,20),0,.5,(255,255,255),1);tiles.append(bgr)
                    cv2.imwrite(str(OUT/f'recurrence_{name}_{anchor:.2f}.png'),np.hstack(tiles))
                print(name,anchor,'background',background,'peaks',[(round(x['lag_s'],2),round(x['score'],2)) for x in candidates],flush=True)
        row_index+=1;picture.cache_clear();del e
    fig.tight_layout();fig.savefig(OUT/'recurrence_curves.png',dpi=150);plt.close(fig);save('recurrence.json',allrows)

if __name__=='__main__':main()
