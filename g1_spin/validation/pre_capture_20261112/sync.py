"""Affine clock mapping from explicitly matched pulses, with temporal holdout."""
import argparse,csv
from pathlib import Path
import numpy as np
from common import save,digest,verify_lock


def fit_sync(path,budget_us):
    rows=list(csv.DictReader(Path(path).open()))
    if len(rows)<6:raise ValueError('Need >=6 matched pulses for fit and holdout')
    if len({r['pulse_id'] for r in rows})!=len(rows):raise ValueError('Duplicate pulse ids')
    a=np.array([int(r['reference_us']) for r in rows],np.int64);b=np.array([int(r['event_us']) for r in rows],np.int64)
    if np.any(a<0) or np.any(b<0) or np.any(np.diff(a)<=0) or np.any(np.diff(b)<=0):raise ValueError('Pulse times must be nonnegative and strictly increasing')
    if not np.isfinite(budget_us) or budget_us<=0:raise ValueError('Invalid synchronization budget')
    # Integer subtraction before conversion protects large absolute timestamps.
    x=(a-a[0]).astype(float);y=(b-b[0]).astype(float);n=max(4,len(rows)*2//3);n=min(n,len(rows)-2)
    train=x[:n];center=train.mean();dx=train-center
    slope=float(dx@(y[:n]-y[:n].mean())/(dx@dx));offset=float(y[:n].mean()-slope*center)
    if slope<=0:raise ValueError('Invalid clock slope')
    residual=y-(slope*x+offset);hold=residual[n:]
    result=dict(mapping='event_us = event_anchor_us + offset_us + slope*(reference_us-reference_anchor_us)',
        reference_anchor_us=int(a[0]),event_anchor_us=int(b[0]),slope=slope,offset_us=offset,drift_ppm=(slope-1)*1e6,
        train_pulses=n,holdout_pulses=len(rows)-n,train_max_abs_residual_us=float(max(abs(residual[:n]))),
        holdout_max_abs_residual_us=float(max(abs(hold))),holdout_rms_us=float(np.sqrt(np.mean(hold**2))),
        budget_us=budget_us,residual_pass=bool(max(abs(residual))<=budget_us),residuals_us=residual.tolist(),
        reference_valid_interval_us=[int(a[0]),int(a[-1])],source_sha256=digest(path),
        limitation='Matched-pulse residual only. Does not certify exposure latency, rolling shutter, marker labeling, or cross-camera geometry. Do not extrapolate beyond pulse span.')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('csv');p.add_argument('--budget-us',type=float,required=True);p.add_argument('--out',required=True);args=p.parse_args();verify_lock()
    if Path(args.out).exists():raise FileExistsError(args.out)
    save(args.out,fit_sync(args.csv,args.budget_us))
