import numpy as np
from contour import extract


def events(points):
    a=np.zeros(len(points),dtype=[('x','f8'),('y','f8')]);a['x'],a['y']=np.array(points).T;return a


def circle(center,radius):
    t=np.linspace(-np.pi,np.pi,720,endpoint=False)
    return center+radius*np.column_stack([np.cos(t),np.sin(t)])


def test_outer_edge_survives_dense_inner_texture_and_isolated_outliers():
    c=np.array([100.,100.]);outer=circle(c,60);inner=circle(c,53)
    raw=events(np.vstack([outer,np.repeat(inner,4,axis=0),[[171.,100.],[100.,171.]]]))
    r=extract(raw,c,60,[0,1]);p=np.array(r['points'])
    assert len(p)>=24
    assert np.max(abs(np.linalg.norm(p-c,axis=1)-60))<1e-10
    assert all(np.any(np.all(outer==v,axis=1)) for v in p)


def test_repeated_hotspot_and_sparse_noise_rejected():
    raw=events([[160.,100.]]*1000+[[100.,160.],[40.,100.]])
    assert extract(raw,[100.,100.],60,[0,1])['points']==[]


def test_translated_outline_keeps_observed_coordinates():
    actual=np.array([103.,98.]);p=circle(actual,60)
    r=extract(events(p),[100.,100.],60,[0,1]);kept=np.array(r['points'])
    assert len(kept)>=24
    assert max(abs(np.linalg.norm(kept-actual,axis=1)-60))<1e-10
