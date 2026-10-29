import numpy as np
import pytest
from g1_spin.events import EventBatch, save_npz
from g1_spin.io import source, windows


def batch(times):
    t = np.asarray(times, dtype=np.int64)
    return EventBatch(np.zeros(len(t)), np.zeros(len(t)), t, np.zeros(len(t)))


def test_chunk_boundaries_and_empty_gaps():
    result = list(windows([batch([100, 100, 109]), batch([110]), batch([140])], 10))
    assert [len(b) for _, _, b in result] == [3, 1, 0, 0, 1]
    assert [s for s, _, _ in result] == [100, 110, 120, 130, 140]
    assert np.array_equal(np.concatenate([b.t_us for _, _, b in result]), [100, 100, 109, 110, 140])


def test_sliding_and_explicit_time_range():
    result = list(windows([batch([100, 109]), batch([110, 120, 125])], 20,
                          start_us=100, duration_us=30, step_us=10))
    assert [b.t_us.tolist() for _, _, b in result] == [[100, 109, 110], [110, 120, 125], [120, 125]]
    assert [(s, e) for s, e, _ in result] == [(100, 120), (110, 130), (120, 130)]
    assert len(list(windows([], 10, start_us=100, duration_us=30))) == 3
    assert list(windows([], 10)) == []


def test_cross_block_regression_and_invalid_arguments():
    with pytest.raises(ValueError, match='regression'):
        list(windows([batch([100, 120]), batch([119])], 10))
    for kwargs in [{'window_us': 0}, {'window_us': 10, 'step_us': 0},
                   {'window_us': 10, 'duration_us': -1}, {'window_us': 1.5}]:
        with pytest.raises(ValueError):
            list(windows([], **kwargs))


def test_npz_source(tmp_path):
    path = tmp_path / 'sample.npz'
    save_npz(path, batch([100, 110, 120]))
    info, chunks = source(path, batch_size=2)
    assert info['timestamp_unit'] == 'us'
    assert [len(b) for b in chunks] == [2, 1]


def test_portable_sdk_manifest_without_sdk(tmp_path):
    import json
    from g1_spin.io import source
    dtype = [('t', '<u8'), ('x', '<u2'), ('y', '<u2'), ('p', 'u1')]
    expected = np.array([(2**54+7, 10, 20, 0), (2**54+8, 11, 21, 1),
                         (2**54+9, 12, 22, 1)], dtype=dtype)
    np.save(tmp_path/'part_000000.npy', expected)
    (tmp_path/'manifest.json').write_text(json.dumps(dict(width=1280, height=720, parts=['part_000000.npy'])))
    info, batches = source(tmp_path, batch_size=2)
    batches = list(batches)
    assert info['width'] == 1280
    assert [len(b) for b in batches] == [2, 1]
    np.testing.assert_array_equal(np.concatenate([b.t_us for b in batches]), expected['t'])
    np.testing.assert_array_equal(batches[0].seconds(), [0, 1e-6])
