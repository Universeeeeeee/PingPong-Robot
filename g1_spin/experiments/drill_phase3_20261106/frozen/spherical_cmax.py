"""Orthographic sphere lifting and spherical contrast maximization.

Camera axes are x right, y down, z away; the visible hemisphere has z < 0.
All tuning constants below are engineering choices, not published guarantees.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize


def lift(x, y, cx, cy, r):
    """Lift pixels on the visible hemisphere; reject invalid sphere geometry."""
    if not np.isfinite([cx, cy, r]).all() or r <= 0:
        raise ValueError("Ball center and positive radius must be finite")
    x, y = (np.asarray(x, float) - cx) / r, (np.asarray(y, float) - cy) / r
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be equally sized vectors")
    q = x*x + y*y
    if not np.isfinite(q).all() or np.any(q > 1 + 1e-12):
        raise ValueError("Events lie outside the sphere")
    return np.column_stack((x, y, -np.sqrt(np.maximum(0, 1-q))))


def warp(points, dt_s, omega):
    """Apply R(-omega * dt_s) to each unit point using Rodrigues' formula."""
    points = np.asarray(points, float)
    dt_s = np.asarray(dt_s, float)
    omega = np.asarray(omega, float)
    if points.ndim != 2 or points.shape[1] != 3 or omega.shape != (3,):
        raise ValueError("Expected Nx3 points and a three-component omega")
    if dt_s.shape not in [(), (len(points),)]:
        raise ValueError("Times must be scalar or one per point")
    if not all(np.isfinite(a).all() for a in (points, dt_s, omega)):
        raise ValueError("Geometry must be finite")
    speed = np.linalg.norm(omega)
    if speed < 1e-14:
        return points.copy()
    axis = omega / speed
    angle = -speed * dt_s
    c, s = np.cos(angle)[..., None], np.sin(angle)[..., None]
    return points*c + np.cross(axis, points)*s + (points @ axis)[:, None]*axis*(1-c)


def iwe(points, shape=(64, 128)):
    """Bilinear (latitude, longitude) histogram with periodic longitude.

    Each event contributes mass one even at latitude endpoints. Gaussian
    smoothing uses a 5x5 sigma=1 kernel and preserves the histogram mass.
    """
    points = np.asarray(points, float)
    h, w = shape
    if h < 2 or w < 2:
        raise ValueError("IWE dimensions must be at least two")
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Expected finite Nx3 sphere points")
    theta = np.arcsin(np.clip(points[:, 2], -1, 1))
    phi = np.arctan2(points[:, 1], points[:, 0])
    yy = (theta / np.pi + .5) * (h-1)
    xx = ((phi / (2*np.pi) + .5) * w) % w
    y0, x0 = np.floor(yy).astype(int), np.floor(xx).astype(int)
    fy, fx = yy-y0, xx-x0
    flat = np.zeros(h*w)
    for dy, wy in ((0, 1-fy), (1, fy)):
        for dx, wx in ((0, 1-fx), (1, fx)):
            indices = np.minimum(y0+dy, h-1)*w + (x0+dx) % w
            flat += np.bincount(indices, weights=wy*wx, minlength=h*w)
    return gaussian_filter(flat.reshape(h, w), sigma=1, truncate=2,
                           mode=("reflect", "wrap"))


def _invalid(reason, quality=None):
    return dict(valid=False, omega_rad_s=None, axis_camera=None,
                spin_rad_s=None, spin_rps=None, spin_rpm=None,
                method="s_cmax", quality=quality or {}, failure_reason=reason)


def estimate(points, dt_s, max_rpm=3500, initial=None):
    """Fit omega without correspondences or ground truth.

    A deterministic sphere-axis/speed grid supplies multistart Nelder-Mead.
    This is an engineering initializer, not Sony's optical-flow initializer.
    Quality gates are conservative heuristics, not statistical confidence.
    """
    points, dt_s = np.asarray(points, float), np.asarray(dt_s, float)
    if points.ndim != 2 or points.shape[1] != 3 or dt_s.shape != (len(points),):
        raise ValueError("Expected Nx3 points and N times in seconds")
    if not np.isfinite(max_rpm) or max_rpm <= 0:
        raise ValueError("max_rpm must be positive")
    if not np.isfinite(points).all() or not np.isfinite(dt_s).all():
        raise ValueError("Points and times must be finite")
    if len(points) < 100:
        return _invalid("insufficient_events", {"event_count": len(points)})
    if np.any(np.diff(dt_s) < 0):
        return _invalid("nonmonotonic_timestamps")
    span = np.ptp(dt_s)
    if span <= 1e-6:
        return _invalid("insufficient_time_span")
    if np.max(np.abs(np.linalg.norm(points, axis=1)-1)) > 1e-5:
        raise ValueError("Points must lie on the unit sphere")
    # Limit optimization cost deterministically, with identical events per trial.
    if len(points) > 5000:
        ix = np.linspace(0, len(points)-1, 5000).astype(int)
        points, dt_s = points[ix], dt_s[ix]
    dt_s = dt_s - (dt_s[0]+dt_s[-1])/2
    bound = max_rpm * 2*np.pi/60
    history = []

    def contrast(omega):
        return float(np.var(iwe(warp(points, dt_s, omega)) / len(points)))

    before = contrast(np.zeros(3))
    # 64 nearly uniform directed axes, supplemented by coordinate diagonals.
    k = np.arange(64)
    z = 1-2*(k+.5)/64
    az = k*np.pi*(3-np.sqrt(5))
    axes = np.column_stack((np.sqrt(1-z*z)*np.cos(az), np.sqrt(1-z*z)*np.sin(az), z))
    axes = np.vstack((axes, np.eye(3), -np.eye(3)))
    # Resolve angular displacement at about 0.4 rad along the magnitude grid.
    n_speed = int(np.clip(np.ceil(bound*span/.4), 8, 30))
    candidates = [(before, np.zeros(3))]
    scan = []
    for magnitude in np.linspace(bound/n_speed, bound, n_speed):
        scores = [(contrast(axis*magnitude), axis*magnitude) for axis in axes]
        candidates.extend(scores)
        scan.append({"rpm": float(magnitude*60/(2*np.pi)), "contrast": max(s[0] for s in scores)})
    if initial is not None:
        seed = np.asarray(initial, float)
        if seed.shape != (3,) or not np.isfinite(seed).all():
            raise ValueError("initial must be a finite three-vector")
        if np.linalg.norm(seed) <= bound:
            candidates.append((contrast(seed), seed))
    candidates.sort(key=lambda item: item[0], reverse=True)
    seeds = []
    for score, candidate in candidates:
        if all(np.linalg.norm(candidate-other) * span > .45 for other in seeds):
            seeds.append(candidate)
        if len(seeds) == 8:
            break
    solutions = []
    for seed in seeds:
        local_history = []

        def objective(omega):
            speed = np.linalg.norm(omega)
            if speed > bound:
                return 1 + (speed/bound-1)**2
            value = contrast(omega)
            local_history.append(value)
            return -value / max(before, 1e-12)

        # Explicit simplex prevents zero components receiving microscopic steps.
        simplex = np.vstack((seed, seed+np.eye(3)*max(2., bound/40)))
        result = minimize(objective, seed, method="Nelder-Mead", options={
            "initial_simplex": simplex, "maxiter": 350, "xatol": .05, "fatol": 1e-6})
        score = contrast(result.x)
        solutions.append((score, result, seed))
        history.append(local_history)
    solutions.sort(key=lambda item: item[0], reverse=True)
    score, best, seed = solutions[0]
    omega = best.x
    speed = np.linalg.norm(omega)
    gain = score/max(before, 1e-12)-1
    rivals = [float(np.linalg.norm(item[1].x-omega)*span)
              for item in solutions[1:] if item[0] >= score*.98]
    quality = dict(event_count=len(points), contrast_before=before,
                   contrast_after=score, relative_improvement=gain,
                   converged=bool(best.success), iterations=int(best.nit),
                   initializer="engineering_grid_multistart", initial_omega_rad_s=seed.tolist(),
                   candidate_scan=scan, objective_histories=history,
                   competing_solution_displacement_rad=max(rivals, default=0),
                   caveat="Convergence and contrast do not prove physical correctness")
    reason = None
    if speed*span < .03 or gain < .08:
        reason = "no_rotation_evidence"
    elif not best.success:
        reason = "optimization_not_converged"
    elif speed >= bound*.985:
        reason = "search_boundary"
    elif max(rivals, default=0) > .4:
        reason = "inconsistent_multistart_solutions"
    elif np.linalg.eigvalsh(np.cov(points.T))[-2] < 1e-5:
        reason = "spatially_degenerate_features"
    if reason:
        quality["rejected_candidate_omega_rad_s"] = omega.tolist()
        return _invalid(reason, quality)
    return dict(valid=True, omega_rad_s=omega.tolist(), axis_camera=(omega/speed).tolist(),
                spin_rad_s=float(speed), spin_rps=float(speed/(2*np.pi)),
                spin_rpm=float(speed*60/(2*np.pi)), method="s_cmax",
                quality=quality, failure_reason=None)
