import numpy as np
import pytest
from g1_spin.spherical_cmax import lift, warp, iwe, estimate


def test_lift_and_bad_roi():
    p = lift([10, 12], [20, 20], 10, 20, 4)
    np.testing.assert_allclose(p[0], [0, 0, -1])
    np.testing.assert_allclose(np.linalg.norm(p, axis=1), 1)
    assert p[1, 2] < 0
    with pytest.raises(ValueError):
        lift([0], [0], 0, 0, 0)
    with pytest.raises(ValueError):
        lift([2], [0], 0, 0, 1)


def test_rotation_inverse_and_sign():
    p = np.array([[1., 0, 0], [0, 1., 0]])
    w = [0., 0., np.pi/2]
    np.testing.assert_allclose(warp(p, 1., w)[0], [0, -1, 0], atol=1e-12)
    np.testing.assert_allclose(warp(warp(p, .1, w), -.1, w), p, atol=1e-12)


def test_iwe_mass_and_longitude_boundary():
    p = np.array([[-1., 1e-12, 0], [-1., -1e-12, 0], [0, 0, -1.]])
    a = iwe(p)
    assert a.sum() == pytest.approx(3)
    np.testing.assert_allclose(iwe(p[:1]), iwe(p[1:2]), atol=1e-10)
    assert iwe(np.empty((0, 3))).sum() == 0


def test_invalid_events_and_times():
    assert not estimate(np.empty((0, 3)), np.empty(0))["valid"]
    p = np.repeat([[0, 0, -1.]], 100, axis=0)
    assert estimate(p, np.zeros(100))["failure_reason"] == "insufficient_time_span"
    assert estimate(p, np.arange(100)[::-1]*1e-6)["failure_reason"] == "nonmonotonic_timestamps"


@pytest.mark.parametrize("rpm,axis", [(500, [0.3, -0.5, 0.8]), (1000, [0.3, -0.5, 0.8]), (2000, [-.5, .2, -.8]), (3000, [0, 0, 1])])
def test_known_rotation_recovery(rpm, axis):
    # Independent forward propagation of asymmetrically placed surface marks.
    from scipy.spatial.transform import Rotation
    anchors = np.array([[.30, .15], [-.40, .20], [.1, -.50], [-.25, -.15]])
    sphere = np.column_stack((anchors, -np.sqrt(1-np.sum(anchors**2, axis=1))))
    omega = np.array(axis, dtype=float)
    omega = omega / np.linalg.norm(omega) * (rpm*2*np.pi/60)
    times = np.repeat(np.linspace(0, .012, 150), len(sphere))
    points = Rotation.from_rotvec(times[:, None]*omega).apply(np.tile(sphere, (150, 1)))
    visible = points[:, 2] < -.1
    result = estimate(points[visible], times[visible], max_rpm=3500)
    assert result["valid"], result["failure_reason"]
    assert abs(result["spin_rpm"]/rpm-1) < .05
    axis = np.array(result["axis_camera"])
    assert np.degrees(np.arccos(np.clip(axis @ (omega/np.linalg.norm(omega)), -1, 1))) < 5


@pytest.mark.parametrize("kind", ["static", "noise", "single"])
def test_negative_controls(kind):
    rng = np.random.default_rng(42)
    times = np.linspace(0, .01, 500)
    n = 500 if kind == "noise" else 4
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    points[:, 2] = -np.abs(points[:, 2])
    if kind == "static":
        points = np.tile(points, (125, 1))
    elif kind == "single":
        points = np.tile([0., 0., -1.], (500, 1))
    result = estimate(points, times)
    assert not result["valid"]
    assert result["spin_rpm"] is None
    assert result["axis_camera"] is None
