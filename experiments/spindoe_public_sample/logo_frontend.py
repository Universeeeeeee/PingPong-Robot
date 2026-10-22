"""Short-window spin from tracked print features, without a known dot pattern.

Input circles must describe the actual ball in every frame. Sphere lifting uses
weak perspective. This detects dark internal print, not the semantic word/logo;
other dark marks inside the ball can also be tracked. Invalid windows stay invalid.
"""
import cv2
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation


SIZE, CENTER, RADIUS = 128, 64.0, 56.0


def _normalized_frame(frame, circle):
    cx, cy, radius = circle
    scale = RADIUS / radius
    transform = np.array([[scale, 0, CENTER-scale*cx],
                          [0, scale, CENTER-scale*cy]], dtype=float)
    gray = cv2.cvtColor(frame, cv2.COLOR_RGB2GRAY)
    return cv2.warpAffine(gray, transform, (SIZE, SIZE), flags=cv2.INTER_LINEAR,
                         borderMode=cv2.BORDER_CONSTANT, borderValue=0)


def _feature_mask(gray, occupied):
    yy, xx = np.mgrid[:SIZE, :SIZE]
    inside = (xx-CENTER)**2 + (yy-CENTER)**2 < (0.82*RADIUS)**2
    background = float(np.median(gray[inside]))
    dark = np.uint8((gray < background-25) & inside)*255
    mask = cv2.dilate(dark, np.ones((3, 3), np.uint8))
    mask[~inside] = 0
    for point in occupied:
        cv2.circle(mask, tuple(np.rint(point).astype(int)), 4, 0, -1)
    return mask


def _new_points(gray, occupied):
    corners = cv2.goodFeaturesToTrack(gray, maxCorners=max(0, 50-len(occupied)),
                                     qualityLevel=0.01, minDistance=4,
                                     mask=_feature_mask(gray, occupied), blockSize=3)
    if corners is None:
        return np.empty((0, 2), np.float32)
    cv2.cornerSubPix(gray, corners, (3,3), (-1,-1),
                    (cv2.TERM_CRITERIA_EPS|cv2.TERM_CRITERIA_COUNT, 20, 0.01))
    chosen = list(occupied)
    new = []
    for point in corners[:,0]:
        if np.sum((point-CENTER)**2) >= (0.82*RADIUS)**2:
            continue
        if chosen and np.min(np.linalg.norm(np.asarray(chosen)-point,axis=1)) < 3:
            continue
        chosen.append(point)
        new.append(point)
    return np.asarray(new,dtype=np.float32).reshape(-1,2)


def _lift(points):
    xy = (points.astype(float)-CENTER)/RADIUS
    xy[:, 1] *= -1
    return np.c_[xy, np.sqrt(1-np.sum(xy*xy, axis=1))]


def _original_pixels(points, circle):
    return (points-CENTER)*(circle[2]/RADIUS) + circle[:2]


def _fit_motion(a, b, dt):
    result = {"valid": False, "reason": "too_few_tracks", "omega_rad_s": None,
              "track_count": len(a), "inlier_count": 0, "fit_residual": None,
              "fit_residual_unit": "unit_sphere_chord", "reprojection_rms_px": None,
              "observability": None, "inlier_mask": [False]*len(a)}
    if len(a) < 4:
        return result
    p, q = _lift(a), _lift(b)
    # Correspondences already have track identity; no reference point pattern.
    rng = np.random.default_rng(0)
    proposals = [np.arange(len(p))]
    proposals += [rng.choice(len(p), 3, replace=False) for _ in range(32)]
    best = None
    for sample in proposals:
        if min(np.linalg.svd(p[sample],compute_uv=False)[1], np.linalg.svd(q[sample],compute_uv=False)[1]) < 1e-4:
            continue
        rotation, _ = Rotation.align_vectors(q[sample], p[sample])
        errors = np.linalg.norm(rotation.apply(p)-q, axis=1)
        inliers = errors < 0.025
        score = (int(inliers.sum()), -float(np.median(errors)))
        if best is None or score > best[0]:
            best = (score, rotation, inliers)
    if best is None:
        result["reason"] = "poor_feature_spread"
        return result
    _, rotation, inliers = best
    if inliers.sum() < 4:
        result["reason"] = "inconsistent_tracks"
        return result
    fit = least_squares(lambda vector: (Rotation.from_rotvec(vector).apply(p[inliers])-q[inliers]).ravel(),
                        rotation.as_rotvec(), loss="soft_l1", f_scale=0.01)
    rotation = Rotation.from_rotvec(fit.x)
    errors = np.linalg.norm(rotation.apply(p)-q, axis=1)
    inliers = errors < 0.025
    if inliers.sum() < 4:
        result["reason"] = "inconsistent_tracks"
        return result
    information = np.mean(np.eye(3)[None]-p[inliers, :, None]*p[inliers, None, :], axis=0)
    eigenvalues = np.linalg.eigvalsh(information)
    condition = float(eigenvalues[0]/eigenvalues[-1])
    residual = float(np.sqrt(np.mean(errors[inliers]**2)))
    predicted = rotation.apply(p[inliers])
    reprojection = float(RADIUS*np.sqrt(np.mean(np.sum((predicted[:, :2]-q[inliers, :2])**2, axis=1))))
    result.update(inlier_count=int(inliers.sum()), inlier_mask=inliers.tolist(),
                  fit_residual=residual, reprojection_rms_px=reprojection,
                  observability={"eigenvalue_ratio": condition, "eigenvalues": eigenvalues.tolist()})
    if condition < 0.01:
        result["reason"] = "poor_feature_spread"
    elif residual > 0.015:
        result["reason"] = "large_fit_residual"
    elif rotation.magnitude() > np.deg2rad(20):
        result["reason"] = "interframe_rotation_too_large"
    else:
        result.update(valid=True, reason="ok", omega_rad_s=(fit.x/dt).tolist())
    return result


def estimate_sequence(frames, timestamps_s, circles):
    """Return JSON-compatible pair windows and original-image track records.

    frames: RGB uint8 arrays; each may have a different shape.
    circles: (cx,cy,r) in each frame's original pixels, before normalization.
    timestamps_s: capture times, not the playback FPS of a slow-motion video.
    The fit only accepts small estimated steps. True motion must also be small;
    full turns between frames can alias to an apparently stationary texture.
    """
    times = np.asarray(timestamps_s, dtype=float)
    circles = np.asarray(circles, dtype=float)
    if len(frames) < 2 or times.shape != (len(frames),) or circles.shape != (len(frames), 3):
        raise ValueError("Need >=2 frames, matching timestamps, and Nx3 circles")
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Capture timestamps must be finite and strictly increasing")
    if not np.isfinite(circles).all() or np.any(circles[:, 2] <= 3):
        raise ValueError("Ball circles must be finite with radius >3 pixels")
    for frame in frames:
        if frame.dtype != np.uint8 or frame.ndim != 3 or frame.shape[2] != 3:
            raise ValueError("frames must be RGB uint8 images")
    normalized = [_normalized_frame(frame, circle) for frame, circle in zip(frames, circles)]
    active = np.empty((0, 2), np.float32)
    active_ids = np.empty(0, dtype=int)
    tracks, windows = [], []
    next_id = 0

    def record(frame_index, points, identities):
        native = _original_pixels(points, circles[frame_index])
        for identity, point, original in zip(identities, points, native):
            tracks.append({"track_id": int(identity), "frame": frame_index,
                           "xy_px": original.astype(float).tolist(),
                           "xy_normalized_roi": point.astype(float).tolist()})

    for index in range(len(frames)-1):
        new = _new_points(normalized[index], active) if len(active) < 40 else np.empty((0,2), np.float32)
        new_ids = np.arange(next_id, next_id+len(new))
        next_id += len(new)
        record(index, new, new_ids)
        active = np.concatenate([active, new])
        active_ids = np.r_[active_ids, new_ids]
        old, moved, identities = active[:0], active[:0], active_ids[:0]
        if len(active):
            forward, status, error = cv2.calcOpticalFlowPyrLK(normalized[index], normalized[index+1],
                active[:, None], None, winSize=(7,7), maxLevel=2,
                criteria=(cv2.TERM_CRITERIA_EPS|cv2.TERM_CRITERIA_COUNT, 30, 0.01))
            refined = forward.copy()
            safe = status[:,0].astype(bool) & np.isfinite(forward[:,0]).all(axis=1)
            safe &= (forward[:,0] > 3).all(axis=1) & (forward[:,0] < SIZE-4).all(axis=1)
            if safe.any():
                subset=refined[safe].copy()
                cv2.cornerSubPix(normalized[index+1], subset, (3,3), (-1,-1),
                                (cv2.TERM_CRITERIA_EPS|cv2.TERM_CRITERIA_COUNT, 20, 0.01))
                refined[safe]=subset
            refined[~safe]=active[~safe,None]
            refinement = np.linalg.norm(refined[:,0]-forward[:,0],axis=1)
            forward = refined
            backward, back_status, _ = cv2.calcOpticalFlowPyrLK(normalized[index+1], normalized[index],
                forward, None, winSize=(7,7), maxLevel=2,
                criteria=(cv2.TERM_CRITERIA_EPS|cv2.TERM_CRITERIA_COUNT, 30, 0.01))
            moved_all = forward[:, 0]
            fb = np.linalg.norm(backward[:, 0]-active, axis=1)
            inner = np.sum((moved_all-CENTER)**2, axis=1) < (0.82*RADIUS)**2
            keep = safe & back_status[:,0].astype(bool) & inner & (fb < 0.75) & (error[:,0] < 25) & (refinement < 1.0)
            # Refinement must not assign multiple track identities to one corner.
            chosen = []
            for candidate in np.argsort(fb):
                if not keep[candidate]:
                    continue
                if chosen and np.min(np.linalg.norm(moved_all[chosen]-moved_all[candidate],axis=1)) < 2:
                    keep[candidate] = False
                else:
                    chosen.append(candidate)
            old, moved, identities = active[keep], moved_all[keep], active_ids[keep]
        result = _fit_motion(old, moved, times[index+1]-times[index])
        result.update(start_frame=index, end_frame=index+1,
                      start_time_s=float(times[index]), end_time_s=float(times[index+1]),
                      correspondence_ids=identities.tolist(),
                      start_points_px=_original_pixels(old, circles[index]).astype(float).tolist(),
                      end_points_px=_original_pixels(moved, circles[index+1]).astype(float).tolist())
        windows.append(result)
        record(index+1, moved, identities)
        # Keep only forward/backward-confirmed IDs. New detections get new IDs.
        active, active_ids = moved, identities
    valid = sum(window["valid"] for window in windows)
    return {"status": "ok" if valid == len(windows) else ("partial" if valid else "no_valid_windows"),
            "valid_window_count": valid, "window_count": len(windows), "windows": windows,
            "tracks": tracks, "frame_circles": circles.tolist(),
            "coordinate_convention": "x image-right, y image-up, visible sphere z>0; weak perspective",
            "feature_scope": "dark print inside ball; not semantic logo segmentation",
            "validity_scope": "local track and geometry consistency; not independently verified spin",
            "motion_assumption": "small true interframe rotation with maintained feature identity",
            "temporal_aliasing_excluded": False,
            "normalized_roi": {"size": SIZE, "center": CENTER, "radius": RADIUS},
            "independent_ground_truth_available": False}
