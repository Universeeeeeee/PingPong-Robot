"""Read-only RAW audit using the existing official SDK adapter; derived files only."""
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path.home() / 'g1_spin/scripts'))
from event_io import source, EventStats

ROOT = Path.home() / 'g1_spin'
OUT = ROOT / 'phase1_demo/outputs/drill_recording_qc_20261105'
OUT.mkdir(exist_ok=True, parents=True)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def raw_image(e, width=1280, height=720):
    im = np.full((height, width, 3), 18, dtype=np.uint8)
    for polarity, color in [(0, [75, 175, 255]), (1, [255, 175, 65])]:
        q = e[e['p'] == polarity]
        im[q['y'], q['x']] = color
    return im

summaries = []
for number, path in enumerate(sorted((ROOT/'data/raw').glob('*20261105_*.raw')), 1):
    out = OUT / path.stem
    out.mkdir(exist_ok=True)
    before = sha(path)
    info, batches = source(path, batch_size=1_000_000)
    stats = EventStats(info['width'], info['height'])
    parts = []
    for i, events in enumerate(batches):
        name = f'part_{i:04d}.npy'
        np.save(out / name, events, allow_pickle=False)
        stats.add(events)
        parts.append(name)
    after = sha(path)
    assert before == after
    info.update(parts=parts, statistics=stats.report(), sha256=before,
                raw_unchanged=True, source_size_bytes=path.stat().st_size,
                recording_order=number)
    (out/'manifest.json').write_text(json.dumps(info, indent=2))
    print(json.dumps({k:v for k,v in info.items() if k!='parts'}), flush=True)
    if stats.invalid or stats.regressions:
        raise RuntimeError('Invalid events or clock regression: inspect before plotting')
    events = np.concatenate([np.load(out / p) for p in parts])
    t0 = int(events['t'][0]); duration = (int(events['t'][-1])-t0)*1e-6
    sample_s = np.linspace(.15, duration-.15, 16)
    if number == 4:
        sample_s = np.r_[.2, .8, 1.4, 2., 2.6, 3.2, 3.8, 4.2,
                         np.linspace(4.8, duration-.15, 8)]
    for window_ms in [5]:
        fig, axes = plt.subplots(4, 4, figsize=(20, 12))
        for ax, s in zip(axes.flat, sample_s):
            a,b = np.searchsorted(events['t'], [t0+round(s*1e6), t0+round(s*1e6)+window_ms*1000])
            ax.imshow(raw_image(events[a:b], info['width'], info['height']))
            ax.set_title(f'{s:.3f}s | {window_ms} ms | {b-a} events', fontsize=10)
            ax.set_axis_off()
        fig.suptitle(f'Run {number} {path.stem} | raw polarity | relative to first decoded event')
        fig.tight_layout(); fig.savefig(out/'overview_5ms.jpg', dpi=130); plt.close(fig)
    counts = np.bincount(((events['t']-t0)//10000).astype(int))
    np.save(out/'rate_10ms.npy', counts)
    fig,ax=plt.subplots(figsize=(12,3)); ax.plot(np.arange(len(counts))*.01,counts/10000)
    ax.set(xlabel='Seconds from first CD event',ylabel='Mevents/s',title=path.stem)
    fig.tight_layout(); fig.savefig(out/'event_rate.png',dpi=120); plt.close(fig)
    summaries.append({k:v for k,v in info.items() if k!='parts'})
    del events
(OUT/'summary.json').write_text(json.dumps(summaries,indent=2))
