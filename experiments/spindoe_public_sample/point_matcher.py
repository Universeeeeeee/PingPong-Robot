"""Order-independent matching of observed sphere features to a known pattern.

Feature coordinates may represent dot centers or other repeatable landmarks.
This module does not detect logos or infer an unknown reference pattern.
"""
from itertools import combinations

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation


def _assignment_score(rotation, observed, reference):
    projected = rotation.apply(reference)
    visible = np.flatnonzero(projected[:, 2] >= 0)
    if len(visible) < len(observed):
        return None
    cost = np.linalg.norm(observed[:, None] - projected[visible][None, :], axis=2)
    rows, columns = linear_sum_assignment(cost)
    return reference[visible[columns]], float(np.sqrt(np.mean(cost[rows, columns] ** 2)))


def match_points(observed, geohasher):
    """Try every nonsingular point-pair basis, with one-to-one verification.

    Returns the best candidate, not a calibrated confidence or a true pose.
    Errors are chord distances on a unit sphere. No threshold is relaxed.
    """
    points = np.asarray(observed, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("observed must have shape (N, 3)")
    if not np.isfinite(points).all():
        raise ValueError("observed must contain finite coordinates")
    if len(points) < 4:
        return None
    # Canonical ordering also makes equal-score tie handling permutation invariant.
    points = points[np.lexsort((points[:, 2], points[:, 1], points[:, 0]))]
    candidates = []
    bases_tested = 0
    for i, j in combinations(range(len(points)), 2):
        if np.linalg.norm(np.cross(points[i], points[j])) < 1e-6:
            continue
        for first, second in ((i, j), (j, i)):
            order = [first, second] + [k for k in range(len(points)) if k not in (first, second)]
            _, rotation, _ = geohasher.identify(points[order])
            bases_tested += 1
            assignment = _assignment_score(rotation, points, geohasher.ref_points)
            if assignment is None:
                continue
            # Refine the rigid rotation with distinct matched reference landmarks.
            for _ in range(2):
                matched, _ = assignment
                rotation, _ = Rotation.align_vectors(points, matched)
                assignment = _assignment_score(rotation, points, geohasher.ref_points)
                if assignment is None:
                    break
            if assignment is not None:
                candidates.append((assignment[1], rotation))
    if not candidates:
        return None
    candidates.sort(key=lambda pair: pair[0])
    best_error, best_rotation = candidates[0]
    # Report competing, geometrically different solutions without inventing a
    # probabilistic confidence. Nearby duplicates are the same pose hypothesis.
    competing = [error for error, rotation in candidates[1:]
                 if (rotation * best_rotation.inv()).magnitude() > np.deg2rad(5)]
    return {"rotation": best_rotation, "rmse": best_error,
            "competing_rmse": min(competing) if competing else None,
            "bases_tested": bases_tested}
