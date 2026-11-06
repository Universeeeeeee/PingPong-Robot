"""Raw short-window visibility and diagnostic circle candidates, not axis accuracy."""
import csv
import json
from pathlib import Path

import cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path.home()/'g1_spin/phase1_demo/outputs/drill_recording_qc_20261105'
BOXES = [(490,230,800,510), (470,340,800,650), (430,270,780,590), (440,210,810,530)]
SEEDS = [(640,355), (620,495), (580,425), (580,355)]

def picture(e, box):
    x0,y0,x1,y1=box
    e=e[(e['x']>=x0)&(e['x']<x1)&(e['y']>=y0)&(e['y']<y1)]
    im=np.full((y1-y0,x1-x0,3),18,np.uint8)
    for polarity,color in [(0,[75,175,255]),(1,[255,175,65])]:
        q=e[e['p']==polarity]; im[q['y']-y0,q['x']-x0]=color
    return im,e

all_rows=[]
for number, folder in enumerate(sorted(ROOT.glob('DSP*')),1):
    meta=json.loads((folder/'manifest.json').read_text())
    events=np.concatenate([np.load(folder/p) for p in meta['parts']])
    t0=int(events['t'][0]); duration=(int(events['t'][-1])-t0)*1e-6
    box=BOXES[number-1]; x0,y0,x1,y1=box
    def select(s,ms):
        a,b=np.searchsorted(events['t'],[t0+round(s*1e6), t0+round(s*1e6)+round(ms*1000)])
        return events[a:b]
    fig,axes=plt.subplots(3,4,figsize=(12,9))
    for j,s in enumerate([1.5,2.,3.,3.6]):
        for i,ms in enumerate([1,2,5]):
            im,q=picture(select(s,ms),box); ax=axes[i,j]
            ax.imshow(im,extent=[x0,x1,y1,y0],interpolation='nearest')
            ax.set_title(f'{s:.3f}s | {ms}ms | ROI {len(q)} events',fontsize=9)
    fig.suptitle(f'Run {number}: raw polarity, fixed crop, no filtering or fitted geometry')
    fig.tight_layout(); fig.savefig(folder/'detail_windows.png',dpi=150); plt.close(fig)
    fig,axes=plt.subplots(2,6,figsize=(16,6))
    for i,ax in enumerate(axes.flat):
        s=2.+i*.002; im,q=picture(select(s,2),box)
        ax.imshow(im,extent=[x0,x1,y1,y0],interpolation='nearest')
        ax.set_title(f'{s:.3f}-{s+.002:.3f}s',fontsize=9)
    fig.suptitle(f'Run {number}: consecutive 2ms raw event windows')
    fig.tight_layout();fig.savefig(folder/'consecutive_2ms.png',dpi=150);plt.close(fig)
    # 50fps overview samples 5ms of every 20ms, so it is explicitly a sampled preview.
    writer=cv2.VideoWriter(str(folder/'sampled_preview_50fps.mp4'),cv2.VideoWriter_fourcc(*'mp4v'),50,(1280,720))
    assert writer.isOpened()
    frame_count=0
    for s in np.arange(0,duration,.02):
        im,q=picture(select(s,5),(0,0,1280,720))
        bgr=cv2.cvtColor(im,cv2.COLOR_RGB2BGR)
        cv2.putText(bgr,f'Run {number} {s:.3f}s | 5ms sampled every 20ms',(15,30),0,.7,(255,255,255),1)
        writer.write(bgr);frame_count+=1
    writer.release()
    cap=cv2.VideoCapture(str(folder/'sampled_preview_50fps.mp4')); decoded=0
    while True:
        ok,frame=cap.read()
        if not ok:break
        assert frame.shape[:2]==(720,1280);decoded+=1
    cap.release(); assert decoded==frame_count
    (folder/'preview_validation.json').write_text(json.dumps(dict(frames_written=frame_count,
        frames_decoded=decoded, fps=50, sampled_window_ms=5, sample_step_ms=20,
        note='Overview samples 25% of event time; raw data remain complete.')))
    # Hough is used only as a diagnostic candidate, never as independent center truth.
    rows=[]; stable_end=3.8 if number==4 else duration-.1
    stable_start=1.2 if number==3 else .4
    for s in np.arange(stable_start,stable_end,.05):
        raw=select(s,5);im,q=picture(raw,box)
        gray=np.zeros((y1-y0,x1-x0),np.uint8);gray[q['y']-y0,q['x']-x0]=255
        gray=cv2.GaussianBlur(gray,(5,5),1)
        circles=cv2.HoughCircles(gray,cv2.HOUGH_GRADIENT,dp=1,minDist=70,
                                param1=80,param2=20,minRadius=45,maxRadius=75)
        row=dict(run=number,start_s=float(s),roi_events=len(q),full_events=len(raw),
                 cx=None,cy=None,radius=None,arc_bins=None)
        if circles is not None:
            circles=circles[0]; centers=circles[:,:2]+[x0,y0]
            seed=np.array(SEEDS[number-1]); good=np.linalg.norm(centers-seed,axis=1)<80
            choices=[]
            for c in circles[good]:
                cx,cy,r=c+[x0,y0,0]
                dx=q['x'].astype(float)-cx;dy=q['y'].astype(float)-cy
                near=abs(np.hypot(dx,dy)-r)<2
                bins=np.unique(np.floor((np.arctan2(dy[near],dx[near])+np.pi)/(2*np.pi)*24).astype(int)%24)
                choices.append((len(bins),cx,cy,r))
            if choices:
                bins,cx,cy,r=max(choices)
                row.update(cx=float(cx),cy=float(cy),radius=float(r),arc_bins=bins)
        rows.append(row)
    with (folder/'circle_candidates.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    # Return representative raw events for optional independent local checks.
    raw=select(2.,100);_,crop=picture(raw,box)
    np.savez_compressed(folder/'sample_2s_100ms.npz',events=crop,origin_us=t0)
    valid=[r for r in rows if r['cx'] is not None and r['arc_bins']>=12]
    summary=dict(run=number,window_count=len(rows),circle_candidates_half_bins=len(valid),
                 geometry_status='unvalidated Hough candidates; includes fitting error and motion',
                 center_x_range=[min(r['cx'] for r in valid),max(r['cx'] for r in valid)] if valid else None,
                 center_y_range=[min(r['cy'] for r in valid),max(r['cy'] for r in valid)] if valid else None,
                 diameter_median_px=float(2*np.median([r['radius'] for r in valid])) if valid else None,
                 median_roi_events_5ms=float(np.median([r['roi_events'] for r in rows])),
                 median_roi_fraction=float(np.median([r['roi_events']/r['full_events'] for r in rows])))
    all_rows.append(summary); print(json.dumps(summary),flush=True)
    fig,axes=plt.subplots(2,4,figsize=(13,6))
    for ax,r in zip(axes.flat,[valid[i] for i in np.linspace(0,len(valid)-1,min(8,len(valid)),dtype=int)]):
        im,q=picture(select(r['start_s'],5),box)
        ax.imshow(im,extent=[x0,x1,y1,y0]);ax.add_patch(plt.Circle((r['cx'],r['cy']),r['radius'],fill=False,color='lime',lw=.7))
        ax.set_title(f'{r["start_s"]:.2f}s | r={r["radius"]:.1f}px',fontsize=9)
    fig.suptitle(f'Run {number}: diagnostic Hough circle candidates, NOT ground truth')
    fig.tight_layout();fig.savefig(folder/'circle_review.png',dpi=140);plt.close(fig)
    del events
(ROOT/'detail_summary.json').write_text(json.dumps(all_rows,indent=2))
