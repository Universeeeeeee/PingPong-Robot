"""Render saved results; no estimator calls or changes to frozen protocol."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
models = json.loads((OUT/'models.json').read_text())
rows = list(csv.DictReader((OUT/'results.csv').open()))
variants = [r['variant'] for r in rows[:9]]
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for window in dict.fromkeys(r['window'] for r in rows):
    rr = [r for r in rows if r['window']==window]
    axes[0].plot(range(9), [float(r['directed_angle_to_A_fixed_deg']) for r in rr], '.-', label=window)
    axes[1].plot(range(9), [float(r['candidate_rpm']) for r in rr], '.-', label=window)
for ax in axes:
    ax.set_xticks(range(9), variants, rotation=60, ha='right'); ax.legend(fontsize=6)
axes[0].set_ylabel('Directed candidate axis difference to A (deg)')
axes[1].set_ylabel('Candidate RPM')
fig.suptitle('Unverified center hypotheses: computational validity is not physical accuracy')
fig.tight_layout(); fig.savefig(OUT/'axis_sensitivity.png', dpi=150); plt.close(fig)

fig, axes = plt.subplots(6, 2, figsize=(10, 21))
for axes_row, window in zip(axes, dict.fromkeys(r['window'] for r in rows)):
    row = next(r for r in rows if r['window']==window)
    ev = np.load(OUT/f'{window}_events.npz')
    model = models[row['segment']]['50']; ref = float(row['reference_s'])
    xy = np.column_stack((ev['x'], ev['y']))
    dt = ev['relative_time_us']/1e6-ref
    aligned = xy-dt[:, None]*np.array(model['velocity_px_s'])
    for ax, coords, title in zip(axes_row, [xy, aligned], ['Raw common events', 'B nominal center alignment']):
        ax.scatter(*coords.T, c=dt*1e3, s=.6, cmap='viridis')
        ax.set(title=f'{window}\n{title}', aspect='equal'); ax.invert_yaxis()
fig.suptitle('Unverified B hypothesis; same selected IDs in each row', y=.997)
fig.tight_layout(rect=(0, 0, 1, .985)); fig.savefig(OUT/'real_event_alignment.png', dpi=140); plt.close(fig)

fig = plt.figure(figsize=(12, 7))
for i, window in enumerate(dict.fromkeys(r['window'] for r in rows)):
    ax = fig.add_subplot(2, 3, i+1, projection='3d')
    rr = [r for r in rows if r['window']==window]
    for j, row in enumerate(rr):
        w = np.array([float(row[k]) for k in ['wx_rad_s', 'wy_rad_s', 'wz_rad_s']])
        u = w/np.linalg.norm(w)
        ax.quiver(0, 0, 0, *u, color=plt.cm.tab10(j), linewidth=1)
        ax.text(*u, str(j), fontsize=7)
    ax.set(xlim=(-1, 1), ylim=(-1, 1), zlim=(-1, 1),
           xlabel='x right', ylabel='y down', zlabel='z away', title=window)
fig.suptitle('Directed optimizer axes; unverified geometry, no true-axis comparison')
fig.text(.5, .01, '0:A  1:B50  2:C100  3:cx+2  4:cx-2  5:cy+2  6:cy-2  7:v×0.75  8:v×1.25', ha='center')
fig.tight_layout(rect=(0, .04, 1, .96)); fig.savefig(OUT/'axis_directions.png', dpi=150); plt.close(fig)
