import pytest
from policy import choose,covers


def test_outer_preferred_when_covered():
    assert choose([0,10000,20000],[0,10000,20000],20000)==('old20','outer20')


def test_fallback_to_old20():
    assert choose([0,10000,20000],[],20000)==('old20','old20')


def test_fallback_to_old10_after_gap():
    assert choose([0,40000,50000],[],50000)==('old10','old10')


def test_no_geometry_no_fabrication():
    assert choose([0,40000],[],30000)==('unavailable','unavailable')


def test_internal_gap_and_invalid_timestamps():
    assert not covers([0,10000,50000,60000],0,60000)
    with pytest.raises(ValueError):covers([0,0,20000],0,20000)


def test_static_prior_is_explicit():
    assert choose([],[],20000,static=True)==('old20','old20')
