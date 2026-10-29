import numpy as np
import pytest
from g1_spin.events import EventBatch, save_npz, load_npz
from g1_spin.synthetic import generate
from g1_spin.baseline import estimate_periodicity


def test_integer_time_and_roundtrip(tmp_path):
    batch = EventBatch([1, 2], [2, 3], np.array([10**16, 10**16 + 1]), [0, 1])
    np.testing.assert_allclose(batch.seconds(), [0, 1e-6], rtol=0, atol=0)
    path = tmp_path / 'events.npz'
    save_npz(path, batch)
    assert np.array_equal(load_npz(path).t_us, batch.t_us)


@pytest.mark.parametrize('times', [[2, 1], [1., 2.], [-1, 0]])
def test_invalid_timestamps(times):
    with pytest.raises(ValueError):
        EventBatch([1, 1], [1, 1], times, [0, 1])


def test_coordinates_and_empty():
    with pytest.raises(ValueError):
        EventBatch([-1], [0], [0], [1])
    batch = generate(rpm=0)
    assert len(batch) == 0
    assert not estimate_periodicity(batch)['valid']
    with pytest.raises(ValueError):
        generate(axis=[0, 0, 0])


def test_synthetic_projection():
    batch = generate(rpm=-2000, noise_fraction=.1)
    assert len(batch) > 1000
    batch.validate(width=256, height=256)
    assert set(vars(batch)) == {'x', 'y', 't_us', 'polarity'}


def test_periodicity_and_noise():
    rng = np.random.default_rng(9)
    bins = np.arange(2000)
    for periodic in [True, False]:
        rates = 25 + 20 * np.sin(2 * np.pi * 20 * bins * .0005) if periodic else np.full(len(bins), 25)
        counts = rng.poisson(rates)
        t = np.repeat(bins * 500, counts).astype(np.int64)
        events = EventBatch(np.ones(len(t)), np.ones(len(t)), t, np.zeros(len(t)))
        result = estimate_periodicity(events)
        if periodic:
            assert result['valid']
            assert abs(result['candidate_frequency_hz'] - 20) < 1
        else:
            assert not result['valid']
