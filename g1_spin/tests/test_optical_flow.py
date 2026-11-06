import numpy as np
from g1_spin.events import EventBatch
from g1_spin.optical_flow import solve_constraints, estimate_flow


def test_exact_aperture_constraints_recover_3d_rotation():
    rng = np.random.default_rng(73)
    points = rng.normal(size=(150, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    points[:, 2] = -np.abs(points[:, 2])
    truth = np.array([35., -60., 80.])
    radius = 100.
    projected_velocity = radius * np.cross(truth, points)[:, :2]
    directions = rng.normal(size=(len(points), 2))
    dots = np.einsum('ij,ij->i', directions, projected_velocity)
    keep = np.abs(dots) > 100
    gradients = directions[keep] / dots[keep, None]
    result = solve_constraints(points[keep], gradients, radius)
    assert result['valid']
    np.testing.assert_allclose(result['omega_rad_s'], truth, atol=1e-9)


def test_rank_deficiency_and_empty_input():
    result = solve_constraints(np.tile([0., 0., -1.], (20, 1)), np.tile([.001, 0.], (20, 1)), 100)
    assert not result['valid']
    assert result['failure_reason'] == 'degenerate_flow_system'
    empty = EventBatch([], [], np.array([], dtype=np.int64), [])
    assert not estimate_flow(empty, 128, 128, 100)['valid']
