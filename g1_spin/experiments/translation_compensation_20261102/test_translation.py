"""Contract checks for the experiment's translation data path."""
import numpy as np
import pytest

from translation import centers_at, compensate, load_centers


META = dict(image_plane="G1 image plane", center_semantics="ball_center",
            timebase="synthetic integer SDK-like us", source="controlled geometry",
            synchronization_verified=True, synchronization_basis="same synthetic clock",
            center_error="exact synthetic construction; real error unknown")


def apply(xy, times, sample_times, samples, **kwargs):
    return compensate(xy, times, sample_times, samples, reference_us=1500,
                      metadata=META, event_timebase=META["timebase"],
                      max_gap_us=1000, **kwargs)


def test_static_and_linear_closed_form():
    times = np.array([1000, 1250, 2000], dtype=np.uint64)
    samples = np.array([1000, 2000], dtype=np.uint64)
    xy = np.array([[12.2, 3.5], [10.1, 7.1], [9.2, 8.3]])
    constant = np.array([[11., 21.], [11., 21.]])
    assert np.array_equal(apply(xy, times, samples, constant)[0], xy)
    centers = np.array([[10., 20.], [12., 19.]])
    corrected, reference = apply(xy, times, samples, centers)
    # Independent expression in seconds, not interpolation.
    expected = xy-(times.astype(float)-1500)[:, None]*np.array([.002, -.001])
    np.testing.assert_allclose(corrected, expected, atol=4e-15)
    np.testing.assert_array_equal(reference, [11., 19.5])


def test_uint64_large_origin_preserves_integer_differences():
    base = 2**63+1234567
    times = np.array([1000, 1250, 2000], dtype=np.uint64)
    samples = np.array([1000, 2000], dtype=np.uint64)
    xy = np.ones((3, 2))
    centers = np.array([[10., 20.], [12., 19.]])
    expected, _ = apply(xy, times, samples, centers)
    actual, _ = compensate(xy, times+np.uint64(base), samples+np.uint64(base), centers,
                           reference_us=base+1500, metadata=META,
                           event_timebase=META["timebase"], max_gap_us=1000)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("times,xy,query,gap", [
    ([1000, 2000], [[1, 2], [3, 4]], [999], 1000),
    ([1000, 2000], [[1, 2], [3, 4]], [2001], 1000),
    ([1000, 2000], [[1, 2], [3, 4]], [1500], 999),
    ([1000, 1000], [[1, 2], [3, 4]], [1000], 1000),
    ([2000, 1000], [[1, 2], [3, 4]], [1500], 1000),
    ([1000, 2000], [[1, 2], [np.nan, 4]], [1500], 1000),
])
def test_bad_trajectory_rejected(times, xy, query, gap):
    with pytest.raises(ValueError):
        centers_at(np.array(times), xy, np.array(query), max_gap_us=gap)


@pytest.mark.parametrize("field,value", [("synchronization_verified", False),
                                        ("synchronization_basis", ""),
                                        ("center_semantics", "marker"),
                                        ("timebase", "another clock")])
def test_unusable_declarations_rejected(field, value):
    metadata = dict(META, **{field: value})
    with pytest.raises(ValueError):
        compensate([[1, 2]], np.array([1500]), np.array([1000, 2000]),
                   [[1, 2], [1, 2]], reference_us=1500, metadata=metadata,
                   event_timebase=META["timebase"], max_gap_us=1000)


def test_csv_roundtrip_and_integer_timestamp_requirement(tmp_path):
    import json
    csv_path, meta_path = tmp_path/"centers.csv", tmp_path/"metadata.json"
    csv_path.write_text("time_us,cx_px,cy_px\n1000,1.25,2.5\n2000,2.25,3.5\n")
    meta_path.write_text(json.dumps(META))
    times, xy, metadata = load_centers(csv_path, meta_path)
    corrected, _ = compensate([[5., 6.]], np.array([1500]), times, xy,
                             reference_us=1500, metadata=metadata,
                             event_timebase=META["timebase"], max_gap_us=1000)
    np.testing.assert_array_equal(corrected, [[5., 6.]])
    with pytest.raises(ValueError):
        centers_at(np.array([1000., 2000.]), xy, np.array([1500]), max_gap_us=1000)
