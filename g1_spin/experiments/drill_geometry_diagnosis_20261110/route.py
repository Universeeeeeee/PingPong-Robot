"""Geometry-only offline window selection; no spin scores or estimated RPM."""
import numpy as np


def covered(times_us, begin_us, end_us, max_gap_us=30000):
    ts=np.asarray(times_us)
    if ts.ndim!=1 or len(ts)<2 or np.any(np.diff(ts)<=0):
        raise ValueError('Geometry timestamps must be strictly increasing')
    if begin_us>=end_us:raise ValueError('Empty or inverted interval')
    if begin_us<ts[0] or end_us>ts[-1]:return False
    overlapping=(ts[:-1]<end_us)&(ts[1:]>begin_us)
    return not np.any(np.diff(ts)[overlapping]>max_gap_us)


def select(times_us,end_us):
    if covered(times_us,end_us-20000,end_us):return '20ms_median3'
    if covered(times_us,end_us-10000,end_us):return '10ms_fixed_fallback'
    return 'geometry_unavailable'
