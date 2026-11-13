import pytest
from route import covered,select


def test_recover_last10_after_gap():
    assert select([3515000,3555000,3565000,3575000],3565000)=='10ms_fixed_fallback'


def test_internal_gap_not_hidden_by_exact_endpoints():
    assert not covered([0,10000,50000,60000],0,60000)


def test_gap_before_window_does_not_reject():
    assert covered([0,40000,50000,60000],40000,60000)


def test_outside_coverage_and_full_rejection():
    assert select([0,40000,80000],60000)=='geometry_unavailable'
    assert not covered([10000,20000],0,20000)


def test_regular_and_bad_time_inputs():
    assert select([0,10000,20000],20000)=='20ms_median3'
    with pytest.raises(ValueError):covered([0,0,20000],0,20000)
    with pytest.raises(ValueError):covered([0,10000],10000,0)
