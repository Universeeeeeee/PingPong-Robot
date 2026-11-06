"""Separately predeclared near branch; existing L/A1/A2 observations are immutable inputs."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import exploratory
import review
import validate_review
OUT=review.OUT
CONFIG=json.loads((OUT/'near_protocol.json').read_text())


def main():
    old_path=OUT/'exploratory_observations.csv';before=hashlib.sha256(old_path.read_bytes()).hexdigest()
    frames,_=exploratory.windows();rows=[]
    for duration in [5,2]:
        images=[review.b0.image(f['events']) for f in frames[duration]]
        for feature,seed in CONFIG['seeds'].items():
            rr=exploratory.forward_general(frames[duration],images,seed,duration,feature)
            for r in rr:
                if r['position_valid']:
                    e=frames[duration][int(r['frame'])]['events'];center=[r['anchor_x'],r['anchor_y']]
                    r.update(review.b0.geometry(e,center));r.update(review.morphology(e,center))
            rows.extend(rr)
            print(duration,feature,[(r['frame'],r['status']) for r in rr],flush=True)
    review.write_csv('near_observations.csv',rows)
    assert before==hashlib.sha256(old_path.read_bytes()).hexdigest()
    old_rows=validate_review.read_csv('exploratory_observations.csv');joint=old_rows+rows
    fig,axes=plt.subplots(2,5,figsize=(18,9));colors={'L':'#b03090','N1':'#008b45','N2':'#7755cc'}
    for i,duration in enumerate([5,2]):
        chosen=[10,12,14,16,19] if duration==5 else [25,30,35,40,49]
        for ax,k in zip(axes[i],chosen):
            f=frames[duration][k];review.b0.raw.draw_raw(ax,review.b0.raw.crop(f['events']),2);missing=[]
            for feature in ['L','N1','N2']:
                r=next(r for r in joint if r['duration_ms']==duration and r['feature']==feature and r['frame']==k)
                if r['position_valid']:
                    center=np.array([r['anchor_x'],r['anchor_y']]);ax.plot(*center,'x',color=colors[feature]);ax.text(center[0]+1,center[1]-2,feature,color=colors[feature],fontsize=8)
                    ax.add_patch(plt.Rectangle(center-13,26,26,fill=False,color=colors[feature],lw=.8))
                else:missing.append(feature+' lost')
            ax.set_title(f'NEAR {duration}ms f{k} [{f["start_s"]:.3f},{f["end_s"]:.3f})s',fontsize=8);ax.set_xlabel('; '.join(missing),fontsize=8,color='#aa3333')
    fig.suptitle('Post-survey near branch: independent patch tracks, overlapping support possible\nL trajectory is the unchanged prior output; no auxiliary backfilling',fontsize=12)
    fig.subplots_adjust(left=.035,right=.995,bottom=.06,top=.9,wspace=.18,hspace=.45)
    fig.savefig(OUT/'near_tracking_review.png',dpi=180);plt.close(fig)
    meta=dict(protocol_sha256=hashlib.sha256((OUT/'near_protocol.json').read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              reused_exploratory_script_sha256=hashlib.sha256((OUT/'exploratory.py').read_bytes()).hexdigest(),existing_exploratory_observations_sha256=before,
              original_three_observations_unchanged=True,rows=len(rows),seeds=CONFIG['seeds'])
    (OUT/'near_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')


if __name__=='__main__':main()
