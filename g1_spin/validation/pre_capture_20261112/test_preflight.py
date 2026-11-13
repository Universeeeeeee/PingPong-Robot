import csv
import numpy as np
import pytest
from common import Events,read_geometry,smooth
from sync import fit_sync
from frozen.policy import covers,choose


def pulses(path,x,y):
    with path.open('w') as f:
        w=csv.writer(f);w.writerow(['pulse_id','reference_us','event_us']);w.writerows(zip(range(len(x)),x,y))


def test_sync_offset_drift_large_origin_holdout(tmp_path):
    x=10**15+np.arange(9,dtype=np.int64)*1000000;y=10**15+50000+np.arange(9,dtype=np.int64)*1000100
    p=tmp_path/'pulse.csv';pulses(p,x,y);r=fit_sync(p,1)
    assert abs(r['slope']-1.0001)<1e-12 and r['residual_pass'] and r['holdout_pulses']==3


def test_sync_holdout_detects_late_mismatch(tmp_path):
    x=np.arange(9)*1000000;y=x+100;y[-2:]+=10000
    p=tmp_path/'p.csv';pulses(p,x,y);r=fit_sync(p,1000);assert not r['residual_pass']


def test_bad_sync_rejected(tmp_path):
    p=tmp_path/'p.csv';pulses(p,[0,1,2],[0,1,2])
    with pytest.raises(ValueError):fit_sync(p,100)
    pulses(p,[0,1,1,3,4,5],[0,1,2,3,4,5])
    with pytest.raises(ValueError):fit_sync(p,100)


def test_event_boundaries_and_regression(tmp_path):
    p=tmp_path/'e.npz';np.savez(p,x=[1,2,3,4],y=[1,2,3,4],t_us=[100,110,110,120],polarity=[0,1,0,1])
    e=Events(p,100);raw,ids=e.window(10,20);assert ids.tolist()==[1,2] and len(raw)==2
    np.savez(p,x=[1,2],y=[1,2],t_us=[100,99],polarity=[0,1])
    with pytest.raises(ValueError):Events(p,100)


def test_invalid_reference_and_gap(tmp_path):
    p=tmp_path/'g.csv';p.write_text('t_us,cx_px,cy_px,valid\n0,1,1,1\n10000,,,0\n40000,1,1,1\n50000,1,1,1\n')
    ts,xy=read_geometry(p);assert not covers(ts,20000,40000);assert choose(ts,[],50000)==('old10','old10')
    p.write_text('t_us,cx_px,cy_px,valid\n0,nan,1,1\n')
    with pytest.raises(ValueError):read_geometry(p)


def test_smoother_does_not_cross_gap():
    ts=np.array([0,10000,20000,60000]);xy=np.array([[0.,0.],[5.,5.],[0.,0.],[6.,6.]])
    assert np.array_equal(smooth(ts,xy),[[0,0],[0,0],[0,0],[6,6]])


def test_empty_sensitivity_support_remains_invalid():
    from frozen import cmax
    points=cmax.lift(np.array([]),np.array([]),0,0,60)
    result=cmax.estimate(points,np.array([]),max_rpm=3500)
    assert not result['valid'] and result['spin_rpm'] is None and result['failure_reason']=='insufficient_events'


def test_alignment_large_epoch_and_no_extrapolation(tmp_path):
    from align_reference import align
    import json
    p=tmp_path/'mapping.json';p.write_text(json.dumps(dict(residual_pass=True,reference_valid_interval_us=[10**15,10**15+1000],reference_anchor_us=10**15,event_anchor_us=10**15+100,offset_us=.2,slope=1.0)))
    source=tmp_path/'ref.csv';source.write_text('reference_us,cx_event_px,cy_event_px,valid\n1000000000000010,4,5,1\n1000000000000020,,,0\n')
    out=tmp_path/'g.csv';align(source,p,10**15,out);ts,xy=read_geometry(out);assert ts.tolist()==[110]
    source.write_text('reference_us,cx_event_px,cy_event_px,valid\n1000000000001001,4,5,1\n')
    with pytest.raises(ValueError):align(source,p,10**15,tmp_path/'bad.csv')
