"""Uniform raw-event survey for human appearance review, without semantic detection."""
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
    manifest = json.loads((DATA / 'manifest.json').read_text())
    events = np.concatenate([np.load(DATA / part) for part in manifest['parts']])
    times = events['t'].astype(np.int64)
    t0 = int(times[0])
    midpoints = np.arange(.25, 17.26, .5)
    rows = []
    # Two pages keep the local structures large enough for visual inspection.
    for page, mids in enumerate([midpoints[:20], midpoints[20:]], 1):
        fig, axes = plt.subplots(4 if page == 1 else 3, 5, figsize=(18, 12 if page == 1 else 9))
        for ax, mid in zip(axes.flat, mids):
            lo = t0 + round((mid - .0025) * 1e6)
            hi = t0 + round((mid + .0025) * 1e6)
            left, right = np.searchsorted(times, [lo, hi])
            selected = events[left:right]
            assert np.array_equal(selected, events[(times >= lo) & (times < hi)])
            selected = selected[(selected['x'] >= 485) & (selected['x'] <= 630) &
                                (selected['y'] >= 155) & (selected['y'] <= 250)]
            on = selected[selected['p'] > 0]
            off = selected[selected['p'] <= 0]
            for part, color in [(on, '#d47527'), (off, '#377eb8')]:
                ax.scatter(part['x'], part['y'], s=1.5, marker='s', color=color, linewidths=0)
            a = np.unique(on['y'].astype(int) * 1280 + on['x'])
            b = np.unique(off['y'].astype(int) * 1280 + off['x'])
            dual = np.intersect1d(a, b)
            ax.plot(dual % 1280, dual // 1280, linestyle='none', marker='o', markersize=2,
                    fillstyle='left', markerfacecolor='#d47527', markerfacecoloralt='#377eb8', markeredgewidth=0)
            ax.set(xlim=(485, 630), ylim=(250, 155), aspect='equal',
                   title=f'{mid:.2f}s +/- 2.5ms | {len(selected)} events')
            ax.title.set_fontsize(9)
            ax.tick_params(labelsize=7)
            rows.append(dict(midpoint_s=float(mid), start_s=float(mid-.0025), end_s=float(mid+.0025),
                             start_us=lo, end_us=hi, roi_events=len(selected),
                             distinct_pixels=len(np.union1d(a, b)), independent_time_check=True))
        fig.suptitle('Uniform 0.5s sampling across the recording; each panel is a 5ms raw-event window\n'
                     'ON orange / OFF blue; no logo identity, fitted boundary or feature selection assumed', fontsize=12)
        fig.subplots_adjust(left=.04, right=.99, bottom=.04, top=.89, hspace=.45, wspace=.2)
        fig.savefig(OUT / f'whole_recording_survey_{page}.png', dpi=170)
        plt.close(fig)
    with (OUT / 'whole_recording_survey.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'{len(rows)} uniformly selected raw windows; all absolute-time checks passed')


if __name__ == '__main__':
    main()
