"""Frozen one-candidate contour experiment on previously verified raw crops."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from contour import extract

OUT=Path(__file__).resolve().parent
SOURCE=OUT.parent/'drill_center_evidence_20261110'
s=importlib.util.spec_from_file_location('evidence',SOURCE/'run.py');base=importlib.util.module_from_spec(s);s.loader.exec_module(base)
old=base.old


def main():
    if (OUT/'measurements.json').exists():raise FileExistsError('Preserve output')
    rows=json.loads((SOURCE/'measurements.json').read_text());anchors=json.loads((SOURCE/'anchors.json').read_text())
    protected=[OUT/'run.py',OUT/'contour.py',OUT/'protocol.md',SOURCE/'run.py',SOURCE/'measurements.json',SOURCE/'anchors.json']+[SOURCE/x['input_file'] for x in rows]
    hashes={str(p.relative_to(base.ROOT)):old.digest(p) for p in protected};old.save(OUT/'freeze.json',hashes)
    results=[]
    for a in anchors:
        anchor=np.array(a['center']);radius=a['radius'];direction=np.array(a['direction']);perp=np.array([-direction[1],direction[0]])
        rr=[r for r in rows if r['case']==a['case']]
        fig,axes=plt.subplots(2,9,figsize=(20,6))
        for source in rr:
            raw=np.load(SOURCE/source['input_file'])['events']
            candidate=extract(raw,anchor,radius,direction);points=np.array(candidate['points']).reshape(-1,2);side=(points-anchor)@perp
            fits={key:base.fit_circle(points[mask],anchor,radius,minimum) for key,mask,minimum in [
                 ('all',np.ones(len(points),bool),12),('side_a',side<0,6),('side_b',side>=0,6)]}
            item=dict(case=a['case'],run=a['run'],end_us=a['end_us'],mid_us=source['mid_us'],duration_us=source['duration_us'],
                      input_file=source['input_file'],input_sha256=source['input_sha256'],candidate=candidate,fits=fits,old_fits=source['fits'],shaft=source['shaft'])
            results.append(item)
            row=0 if source['duration_us']==1000 else 1;col=(source['mid_us']-(a['end_us']-50000))//10000;ax=axes[row,col]
            box=np.r_[np.floor(anchor-90).astype(int),np.ceil(anchor+90).astype(int)]
            ax.imshow(base.p2.picture(raw,box),extent=[box[0],box[2],box[3],box[1]])
            if len(points):ax.scatter(points[:,0],points[:,1],s=8,color='red')
            for key,color in [('side_a','lime'),('side_b','magenta')]:
                if fits[key]['valid']:ax.add_patch(plt.Circle(fits[key]['center'],radius,fill=False,color=color,lw=.7))
            ax.set_title(f'{source["mid_us"]/1e6:.3f}s/{source["duration_us"]/1000:g}ms\nA:{fits["side_a"]["bins"]} B:{fits["side_b"]["bins"]}',fontsize=8)
            ax.set_xticks([]);ax.set_yticks([])
        fig.suptitle(a['case']+' | Red: actual selected outer-support pixels; green/magenta: half-circle fits')
        fig.tight_layout();fig.savefig(OUT/(a['case']+'_review.png'),dpi=130);plt.close(fig)
        print(a['case'],flush=True)
    old.save(OUT/'measurements.json',results)
    assert len(results)==144 and all(old.digest(base.ROOT/p)==h for p,h in hashes.items())
    print('144 exposures complete; sources unchanged')

if __name__=='__main__':main()
