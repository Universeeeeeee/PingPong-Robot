"""Offline frontend choice depends only on geometry coverage, never spin output."""
import numpy as np


def covers(ts, begin, end):
    ts=np.asarray(ts)
    if len(ts)<2 or begin<ts[0] or end>ts[-1]:return False
    if np.any(np.diff(ts)<=0):raise ValueError('Nonmonotonic geometry')
    return not np.any((np.diff(ts)>30000)&(ts[:-1]<end)&(ts[1:]>begin))


def choose(old_times,outer_times,end,static=False):
    old='old20' if static or covers(old_times,end-20000,end) else 'old10' if covers(old_times,end-10000,end) else 'unavailable'
    new='outer20' if covers(outer_times,end-20000,end) else old
    return old,new
