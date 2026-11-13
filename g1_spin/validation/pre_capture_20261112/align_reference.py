"""Map already calibrated event-pixel reference observations into event time."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from common import verify_lock


def align(source,mapping,event_origin_us,out):
    m=json.loads(Path(mapping).read_text())
    if not m['residual_pass']:raise ValueError('Synchronization residual check failed')
    rows=list(csv.DictReader(Path(source).open()));low,high=m['reference_valid_interval_us'];mapped=[];last=None
    for r in rows:
        t=int(r['reference_us'])
        if not low<=t<=high:raise ValueError('Reference outside calibrated pulse span')
        # Avoid adding a large float epoch before rounding; anchors remain integers.
        relative=(int(m['event_anchor_us'])-int(event_origin_us))+round(m['offset_us']+m['slope']*(t-int(m['reference_anchor_us'])))
        if relative<0 or (last is not None and relative<=last):raise ValueError('Invalid mapped timestamps')
        if r['valid'] not in ['0','1']:raise ValueError('valid must be 0 or 1')
        if r['valid']=='1' and not np.isfinite([float(r['cx_event_px']),float(r['cy_event_px'])]).all():raise ValueError('Invalid center')
        mapped.append([relative,r['cx_event_px'],r['cy_event_px'],r['valid']]);last=relative
    with Path(out).open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['t_us','cx_px','cy_px','valid']);w.writerows(mapped)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('csv');p.add_argument('--mapping',required=True);p.add_argument('--event-origin-us',type=int,required=True);p.add_argument('--out',required=True);a=p.parse_args();verify_lock();align(a.csv,a.mapping,a.event_origin_us,a.out)
