"""B0 raw continuous previews. No candidate labels or tracking in this script."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DATA = ROOT/'g1_spin/data/decoded/local_analysis_222403'
SEGMENTS = {'early_turn': (1.70, 1.80), 'middle': (6.30, 6.40), 'late': (11.30, 11.40)}
ROI = (485, 630, 155, 250)


def load_events():
    manifest = json.loads((DATA/'manifest.json').read_text())
    events = np.concatenate([np.load(DATA/f) for f in manifest['parts']])
    return events, manifest


def crop(events):
    a,b,c,d = ROI
    return events[(events['x']>=a)&(events['x']<=b)&(events['y']>=c)&(events['y']<=d)]


def draw_raw(ax, events, point_size=1.):
    on, off = events[events['p']>0], events[events['p']<=0]
    ax.scatter(on['x'],on['y'],s=point_size,color='#dc7424',marker='s',linewidths=0,rasterized=True)
    ax.scatter(off['x'],off['y'],s=point_size,color='#377eb8',marker='s',linewidths=0,rasterized=True)
    ids = lambda e: e['y'].astype(np.int64)*1280+e['x'].astype(np.int64)
    dual = np.intersect1d(np.unique(ids(on)),np.unique(ids(off)))
    ax.plot(dual%1280,dual//1280,linestyle='none',marker='o',markersize=1.4,fillstyle='left',
            markerfacecolor='#dc7424',markerfacecoloralt='#377eb8',markeredgewidth=0,rasterized=True)
    ax.set(xlim=(ROI[0]-.5,ROI[1]+.5),ylim=(ROI[3]+.5,ROI[2]-.5),aspect='equal')


def main():
    events,manifest = load_events()
    times=events['t'].astype(np.int64);t0=int(times[0]);rows=[]
    for segment,(start,end) in SEGMENTS.items():
        # Full sensor timestamps are verified before ROI cropping, independent of searchsorted.
        segment_events=events[(times>=t0+round(start*1e6))&(times<t0+round(end*1e6))]
        for duration in [2,5]:
            n=100//duration;columns=10 if duration==2 else 5;nr=n//columns
            fig,axes=plt.subplots(nr,columns,figsize=(20,12 if duration==2 else 11))
            for k,ax in enumerate(axes.flat):
                lo=t0+round(start*1e6)+k*duration*1000;hi=lo+duration*1000
                a,b=np.searchsorted(times,[lo,hi]);selected=events[a:b]
                direct=segment_events[(segment_events['t']>=lo)&(segment_events['t']<hi)]
                assert np.array_equal(selected,direct)
                raw=crop(selected)
                unique=len(np.unique(raw['y'].astype(np.int64)*1280+raw['x']))
                rows.append(dict(segment=segment,duration_ms=duration,frame=k,start_s=(lo-t0)/1e6,
                                 end_s=(hi-t0)/1e6,midpoint_s=((lo+hi)/2-t0)/1e6,
                                 start_absolute_us=lo,end_absolute_us=hi,events_sensor=len(selected),
                                 events_roi=len(raw),distinct_pixels_roi=unique,
                                 on_events_roi=int(np.sum(raw['p']>0)),off_events_roi=int(np.sum(raw['p']<=0)),
                                 independent_absolute_mask_matches=True))
                draw_raw(ax,raw,1.1 if duration==2 else 1.5)
                ax.set_title(f'{k:02d} [{(lo-t0)/1e6:.3f},{(hi-t0)/1e6:.3f})s\nN={len(raw)} / pixels={unique}',fontsize=7 if duration==2 else 8)
                ax.tick_params(labelsize=5 if duration==2 else 6)
            fig.suptitle(f'{segment}: consecutive, non-overlapping {duration}ms windows, {start:.2f}–{end:.2f}s\n'
                         'Raw ON orange / OFF blue; same-pixel dual polarity split color; no candidates or fits',fontsize=12)
            fig.subplots_adjust(left=.025,right=.995,bottom=.04,top=.90,hspace=.58 if duration==2 else .5,wspace=.22)
            fig.savefig(OUT/f'continuous_{segment}_{duration}ms.png',dpi=180)
            plt.close(fig)
    with (OUT/'window_counts.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    metadata=dict(source=str(DATA),source_t0_us=t0,events_total=len(events),segments=SEGMENTS,
                  roi_inclusive=ROI,durations_ms=[2,5],nonoverlapping=True,expected_windows=210,
                  actual_windows=len(rows),all_independent_absolute_masks_match=True,
                  manifest_sha256=hashlib.sha256((DATA/'manifest.json').read_bytes()).hexdigest(),
                  input_parts_sha256={f:hashlib.sha256((DATA/f).read_bytes()).hexdigest() for f in manifest['parts']},
                  continuous_preview_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    assert len(rows)==210
    (OUT/'preview_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Wrote six continuous raw contact sheets and 210 non-overlapping window counts; all masks verified.')


if __name__=='__main__':main()
