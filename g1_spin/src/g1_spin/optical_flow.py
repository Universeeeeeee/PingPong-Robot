"""Engineering normal-flow initializer using local event time planes.

The scalar aperture constraint is r*(s cross [a,b,0]) dot omega = 1
for a fitted t_seconds = a*u_pixels + b*v_pixels + c. No full 2D
velocity is invented from a single edge. Thresholds are engineering choices.
"""
import numpy as np
from scipy.spatial import cKDTree


def _invalid(reason, quality=None):
    return dict(valid=False, method='local_plane_normal_flow', omega_rad_s=None,
                axis_camera=None, spin_rad_s=None, spin_rps=None, spin_rpm=None,
                quality=quality or {}, failure_reason=reason)


def solve_constraints(points, gradients, radius):
    """Solve aperture equations from unit sphere points and seconds/pixel slopes."""
    points, gradients = np.asarray(points, float), np.asarray(gradients, float)
    if radius <= 0 or points.ndim != 2 or points.shape[1] != 3 or gradients.shape != (len(points), 2):
        raise ValueError('Expected Nx3 points, Nx2 gradients and positive radius')
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(gradients)):
        raise ValueError('Nonfinite flow constraints')
    if len(points) < 12:
        return _invalid('insufficient_flow_constraints')
    normals = np.column_stack([gradients, np.zeros(len(points))])
    matrix = radius * np.cross(points, normals)
    omega, _, rank, singular = np.linalg.lstsq(matrix, np.ones(len(points)), rcond=None)
    condition = float(singular[0] / singular[-1]) if len(singular) == 3 and singular[-1] > 1e-15 else None
    residual = float(np.sqrt(np.mean((matrix @ omega - 1) ** 2)))
    quality = dict(constraints=len(points), rank=int(rank), condition=condition,
                   constraint_rms=residual, thresholds='engineering choices: rank=3, condition<=100, RMS<=0.35')
    if rank < 3 or condition is None or condition > 100:
        return _invalid('degenerate_flow_system', quality)
    if residual > .35:
        return _invalid('inconsistent_flow_constraints', quality)
    magnitude = float(np.linalg.norm(omega))
    if magnitude < 1e-6:
        return _invalid('no_rotation_evidence', quality)
    return dict(valid=True, method='local_plane_normal_flow', omega_rad_s=omega.tolist(),
                axis_camera=(omega / magnitude).tolist(), spin_rad_s=magnitude,
                spin_rps=magnitude / (2*np.pi), spin_rpm=magnitude * 60 / (2*np.pi),
                quality=quality, failure_reason=None)


def estimate_flow(events, cx, cy, r):
    """Fit up to 300 local planes in 4px/±2ms neighborhoods, CPU only.

    Same-polarity neighborhoods need >=10 samples and 2D spatial support.
    Point-trajectory synthetic events often lack that support and should fail.
    """
    if not np.all(np.isfinite([cx, cy, r])) or r <= 0:
        raise ValueError('Invalid ball ROI')
    if len(events) < 100:
        return _invalid('insufficient_events')
    events.validate()
    xy = np.column_stack([events.x, events.y]).astype(float)
    times = events.seconds()
    tree = cKDTree(xy)
    points, gradients = [], []
    for i in np.linspace(0, len(events)-1, min(300, len(events)), dtype=int):
        indices = np.asarray(tree.query_ball_point(xy[i], 4), dtype=int)
        indices = indices[(np.abs(times[indices] - times[i]) <= .002) &
                          (events.polarity[indices] == events.polarity[i])]
        if len(indices) < 10:
            continue
        coords = xy[indices]
        center = coords.mean(axis=0)
        centered = coords - center
        singular = np.linalg.svd(centered, compute_uv=False)
        if singular[0] < 1e-6 or singular[1] / singular[0] < .08:
            continue
        dt = times[indices] - times[i]
        if np.std(dt) < 1e-6:
            continue
        design = np.column_stack([centered, np.ones(len(indices))])
        coeff, _, rank, _ = np.linalg.lstsq(design, dt, rcond=None)
        rms = float(np.sqrt(np.mean((design @ coeff - dt)**2)))
        if rank < 3 or rms > max(20e-6, .15 * np.std(dt)):
            continue
        sx, sy = (center - [cx, cy]) / r
        if sx*sx + sy*sy >= .9**2:
            continue
        if np.linalg.norm(coeff[:2]) < 1e-8:
            continue
        points.append([sx, sy, -np.sqrt(1-sx*sx-sy*sy)])
        gradients.append(coeff[:2])
    if len(points) < 12:
        return _invalid('insufficient_stable_local_planes', {'accepted_planes': len(points)})
    result = solve_constraints(points, gradients, r)
    result['quality']['accepted_planes'] = len(points)
    result['quality']['model'] = 'local normal flow; aperture constraints; engineering variant'
    return result
