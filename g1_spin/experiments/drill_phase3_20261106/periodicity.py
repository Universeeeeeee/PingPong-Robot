"""Four frozen raw texture recurrence attempts; never reads a CMax result."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import OUT,load,select,picture,digest,write_json,table,array_hash

ATTEMPTS=[(1,1.,2.),(1,4.5,5.5),(4,1.,2.),(4,2.5,3.5)]
CENTERS={1:[642.,357.],4:[598.,355.]}

def main():
    if (OUT/'periodicity.json').exists():raise FileExistsError('Preserve independently selected peaks')
    rows=[];fig,axes=plt.subplots(4,2,figsize=(13,14));(OUT/'periodicity').mkdir(exist_ok=True)
    for (run,begin,end),(rawax,ax) in zip(ATTEMPTS,axes):
        e,t0=load(run);c=np.array(CENTERS[run]);r=60.;vectors=[];counts=[];hashes=[]
        for t in range(round(begin*1e6),round(end*1e6),2000):
            q,a=select(e,t0,t,2000);xy=np.column_stack([q['x'],q['y']]).astype(float)-c
            keep=np.sum(xy*xy,axis=1)<=(.65*r)**2;xy=xy[keep];pol=q['p'][keep]
            bins=np.floor((xy+40)/4).astype(int);frame=np.zeros((2,20,20))
            np.add.at(frame,(pol,bins[:,1],bins[:,0]),1)
            vector=frame.ravel();norm=np.linalg.norm(vector);vectors.append(vector/norm if norm>0 else vector);counts.append(int(keep.sum()))
            hashes.append(array_hash(np.flatnonzero(keep)+a))
        vectors=np.array(vectors);valid=np.array(counts)>0;lags=np.arange(10,251);scores=[];pairs=[]
        for lag in lags:
            supported=valid[:-lag]&valid[lag:];pairs.append(int(supported.sum()))
            scores.append(float(np.mean(np.sum(vectors[:-lag][supported]*vectors[lag:][supported],axis=1))) if supported.sum()>=30 else None)
        peaks=[]
        for i in range(1,len(lags)-1):
            if all(scores[j] is not None for j in [i-1,i,i+1]) and scores[i]>scores[i-1] and scores[i]>=scores[i+1]:
                peaks.append(dict(period_s=float(lags[i]*.002),frequency_hz=float(1/(lags[i]*.002)),rpm_if_one_repeat_per_rotation=float(60/(lags[i]*.002)),similarity=scores[i],pairs=pairs[i]))
        peaks.sort(key=lambda p:p['similarity'],reverse=True)
        chosen=peaks[0] if peaks else None;harmonics=[]
        for i,a in enumerate(peaks):
            for j,b in enumerate(peaks):
                if j>i:
                    ratio=max(a['period_s'],b['period_s'])/min(a['period_s'],b['period_s'])
                    if abs(ratio-2)<=.08:harmonics.append(dict(peak_indices=[i,j],period_ratio=ratio))
        name=f'run{run}_{begin:g}_{end:g}';path=OUT/'periodicity'/f'{name}.npz'
        np.savez_compressed(path,vectors=vectors,counts=np.array(counts),valid=valid,lag_s=lags*.002,scores=np.array([np.nan if s is None else s for s in scores]),pairs=np.array(pairs))
        row=dict(attempt=name,run=run,start_s=begin,end_s=end,status='raw_recurrence_peaks_not_rpm_truth' if chosen else 'unidentifiable',
            fixed_roi_center=c.tolist(),radius_px=r,frames=len(vectors),valid_frames=int(valid.sum()),
            texture_input_hash=array_hash(vectors,np.array(counts)),frame_event_id_hashes=hashes,all_peaks=peaks,selected_by_similarity_only=chosen,
            half_double_peak_pairs=harmonics,input_file=str(path.relative_to(OUT)),input_file_sha256=digest(path),
            caveat='Fixed ROI does not compensate vibration; recurrence may be a subperiod, multiple period or activity pattern. Peaks are not independent RPM truth.')
        rows.append(row)
        q,_=select(e,t0,round((begin+.05)*1e6),2000);box=np.r_[np.floor(c-65).astype(int),np.ceil(c+65).astype(int)]
        rawax.imshow(picture(q,box),extent=[box[0],box[2],box[3],box[1]]);rawax.add_patch(plt.Circle(c,.65*r,fill=False,color='lime',lw=.7));rawax.set_title(f'{name} | one raw2ms example | fixed ROI')
        ax.plot(lags*.002,[np.nan if s is None else s for s in scores]);ax.scatter([p['period_s'] for p in peaks],[p['similarity'] for p in peaks],s=12)
        if chosen:ax.axvline(chosen['period_s'],color='orange',ls='--')
        ax.set(title='All recurrence peaks; selected by similarity before CMax',xlabel='Lag (s)',ylabel='Mean frame cosine similarity');ax.grid(alpha=.2)
        print('periodicity',name,'peaks',len(peaks),'selected',None if chosen is None else chosen['period_s'],flush=True)
    write_json(OUT/'periodicity.json',rows);table(OUT/'periodicity.csv',[{k:v for k,v in r.items() if k not in ['frame_event_id_hashes','all_peaks','half_double_peak_pairs']} for r in rows])
    fig.tight_layout();fig.savefig(OUT/'periodicity_raw_peaks.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
