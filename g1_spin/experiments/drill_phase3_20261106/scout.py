"""One-pass raw geometry scout, before choosing block priors or spin fitting."""
import importlib.util
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stage2_readonly',OUT.parent/'drill_phase2_20261106/run.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
BLOCKS={4:[4.4,5.4,6.4,7.4,8.4],5:[8.,9.5,11.]}

def main():
    for run,blocks in BLOCKS.items():
        e,t0=base.load(run);fig,axes=plt.subplots(len(blocks),2,figsize=(17,4.3*len(blocks)))
        for i,b in enumerate(blocks):
            for j,duration in enumerate([1000,5000]):
                q,_=base.select(e,t0,round((b+.1)*1e6)-duration//2,duration)
                ax=axes[i,j];ax.imshow(base.picture(q,(0,0,1280,720)),extent=[0,1280,720,0])
                ax.set_title(f'run{run} block {b:g} | nominal {b+.1:g}s | raw {duration/1000:g}ms')
        fig.tight_layout();fig.savefig(OUT/f'run{run}_raw_scout.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
