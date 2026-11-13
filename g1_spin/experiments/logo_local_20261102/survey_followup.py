"""Exploratory follow-up of three irregular structures in the uniform survey."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
DATA = OUT.parents[2] / 'g1_spin/data/decoded/local_analysis_222403'

def main():
    manifest = json.loads((DATA/'manifest.json').read_text())
    events = np.concatenate([np.load(DATA/p) for p in manifest['parts']])
    times = events['t'].astype(np.int64); t0 = int(times[0]); rows=[]
    fig, axes = plt.subplots(3, 3, figsize=(13, 12))
    for i, center in enumerate([2.75, 5.25, 11.75]):
        for j, offset in enumerate([-.005, 0, .005]):
            mid = center+offset; lo=t0+round((mid-.0025)*1e6); hi=t0+round((mid+.0025)*1e6)
            a,b=np.searchsorted(times,[lo,hi]); e=events[a:b]
            assert np.array_equal(e,events[(times>=lo)&(times<hi)])
            e=e[(e['x']>=485)&(e['x']<=630)&(e['y']>=155)&(e['y']<=250)]
            for pol,color in [(0,'#377eb8'),(1,'#d47527')]:
                s=e[e['p']==pol];axes[i,j].scatter(s['x'],s['y'],s=4,color=color,marker='s',linewidths=0)
            axes[i,j].set(xlim=(485,630),ylim=(250,155),aspect='equal',title=f'{mid:.3f}s +/-2.5ms; N={len(e)}')
            rows.append(dict(survey_center_s=center,midpoint_s=mid,start_us=lo,end_us=hi,roi_events=len(e)))
    fig.suptitle('Post-survey exploratory inspection of irregular event structures\nNo confirmed logo identity; three consecutive disjoint 5ms blocks per row',fontsize=12)
    fig.subplots_adjust(top=.9,bottom=.05,left=.06,right=.98,hspace=.4,wspace=.25)
    fig.savefig(OUT/'irregular_survey_details.png',dpi=180);plt.close(fig)
    with (OUT/'survey_followup.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':main()
