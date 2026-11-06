"""Raw short-window previews only: no fits, candidates, associations or Marker labels."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DATA = ROOT / 'g1_spin/data/decoded/local_analysis_222403'
SAMPLES = [('early_turn', 1.75), ('middle', 6.35), ('late', 11.35)]
DURATIONS = [1, 2, 5, 10, 20]
# Inclusive integer pixel limits; all panels use these exact limits and equal aspect.
ROI = (485, 630, 155, 250)


def main():
    manifest = json.loads((DATA/'manifest.json').read_text())
    events = np.concatenate([np.load(DATA/f) for f in manifest['parts']])
    times = events['t'].astype(np.int64)
    t0 = int(times[0])
    rows = []
    fig, axes = plt.subplots(3, 5, figsize=(20, 12))
    for i, (segment, midpoint) in enumerate(SAMPLES):
        for j, duration in enumerate(DURATIONS):
            start = t0 + round((midpoint-duration/2000)*1e6)
            end = t0 + round((midpoint+duration/2000)*1e6)
            a, b = np.searchsorted(times, [start, end])
            selected = events[a:b]
            direct = events[(times >= start) & (times < end)]
            assert np.array_equal(selected, direct)
            xmin, xmax, ymin, ymax = ROI
            raw = selected[(selected['x'] >= xmin) & (selected['x'] <= xmax) &
                           (selected['y'] >= ymin) & (selected['y'] <= ymax)]
            direct_roi = direct[(direct['x'] >= xmin) & (direct['x'] <= xmax) &
                                (direct['y'] >= ymin) & (direct['y'] <= ymax)]
            assert np.array_equal(raw, direct_roi)
            on = raw[raw['p'] > 0]
            off = raw[raw['p'] <= 0]
            pixel_id = lambda x: x['y'].astype(np.int64)*1280+x['x'].astype(np.int64)
            on_ids, off_ids = np.unique(pixel_id(on)), np.unique(pixel_id(off))
            dual = np.intersect1d(on_ids, off_ids)
            distinct = len(np.union1d(on_ids, off_ids))
            rows.append(dict(segment=segment, midpoint_s=midpoint, duration_ms=duration,
                             start_relative_s=(start-t0)/1e6, end_relative_s=(end-t0)/1e6,
                             start_absolute_us=start, end_absolute_us=end, source_t0_us=t0,
                             interval='[start,end)', roi_x_min=xmin, roi_x_max_inclusive=xmax,
                             roi_y_min=ymin, roi_y_max_inclusive=ymax,
                             events_full_sensor=len(selected), events_roi=len(raw),
                             distinct_pixels_roi=distinct, on_events_roi=len(on), off_events_roi=len(off),
                             on_distinct_pixels_roi=len(on_ids), off_distinct_pixels_roi=len(off_ids),
                             dual_polarity_pixels_roi=len(dual), independent_absolute_mask_matches=True))
            ax = axes[i, j]
            # Every event is plotted. Repeated pixels overlap without count-dependent point size.
            ax.scatter(on['x'], on['y'], s=1.5, color='#dc7424', marker='s', linewidths=0, rasterized=True)
            ax.scatter(off['x'], off['y'], s=1.5, color='#377eb8', marker='s', linewidths=0, rasterized=True)
            # Same-pixel ON/OFF occupancy is explicitly visible in both colors at the actual pixel.
            ax.plot(dual % 1280, dual // 1280, linestyle='none', marker='o', markersize=2,
                    fillstyle='left', markerfacecolor='#dc7424', markerfacecoloralt='#377eb8',
                    markeredgewidth=0, rasterized=True)
            ax.set(xlim=(xmin-.5, xmax+.5), ylim=(ymax+.5, ymin-.5), aspect='equal',
                   title=f'{segment} {duration}ms [{(start-t0)/1e6:.4f},{(end-t0)/1e6:.4f})s\n'
                         f'{len(raw)} events / {distinct} distinct pixels')
            ax.title.set_fontsize(9)
            ax.tick_params(labelsize=7)
            if j == 0:
                ax.set_ylabel('y (px)')
            if i == 2:
                ax.set_xlabel('x (px)')
    legend = [Line2D([], [], color='#dc7424', marker='s', linestyle='none', markersize=4, label='ON'),
              Line2D([], [], color='#377eb8', marker='s', linestyle='none', markersize=4, label='OFF'),
              Line2D([], [], marker='o', linestyle='none', markersize=5, fillstyle='left',
                     markerfacecolor='#dc7424', markerfacecoloralt='#377eb8', markeredgewidth=0,
                     label='Both polarities at same pixel (split color)')]
    fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(.5, .94), ncol=3, fontsize=9)
    fig.suptitle('Raw short windows at identical midpoints and ROI; all panels share pixel scale\n'
                 'No fitted circles, candidates, automatic Marker labels or associations; repeated events overlap', fontsize=12)
    fig.subplots_adjust(left=.04, right=.995, bottom=.055, top=.86, wspace=.14, hspace=.5)
    fig.savefig(OUT/'short_window_events.png', dpi=180)
    plt.close(fig)
    with (OUT/'short_window_counts.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'Wrote 15 panels and 15 CSV rows; independent absolute-mask checks all pass; t0={t0}us')


if __name__ == '__main__':
    main()
