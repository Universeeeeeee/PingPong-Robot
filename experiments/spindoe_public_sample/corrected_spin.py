"""Sign-consistent adapter for the pinned upstream SpinRegressor."""

import numpy as np
from scipy.optimize import minimize_scalar

from upstream.python.spin_regressor import SpinRegressor


class CorrectedSpinRegressor(SpinRegressor):
    """RANSAC results use rad/s with R(t)=Exp(omega*dt) R(0).

    The upstream estimator's matching, RANSAC and phase unwrapping are unchanged.
    Its extra minus sign is removed from RANSAC_regress. The unused non-RANSAC
    regress method is not adapted. This does not resolve temporal aliasing or
    incorrect input poses. The quaternion-plane regression degenerates for a
    constant pose sequence, so this adapter rejects it. This is a limitation of
    this backend, not a claim that stationary textured balls cannot be measured.
    """

    def RANSAC_regress(self, t, rots):
        t = np.asarray(t, dtype=float)
        if len(t) != len(rots):
            raise ValueError("Timestamps and poses must have equal lengths")
        # Upstream requires n seed poses and strictly more than d extra inliers.
        if len(t) < self.n + self.d + 1:
            return None, None
        if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
            raise ValueError("Timestamps must be finite and strictly increasing")
        motion = np.array([(rot * rots[0].inv()).magnitude() for rot in rots])
        if np.max(motion) <= 1e-10:
            return None, None
        spin, indices = super().RANSAC_regress(t, rots)
        if spin is None:
            return None, None
        return -spin, indices


class BoundedSpinRegressor(SpinRegressor):
    """Fit circular phase within an explicit speed bound, without unwrapping.

    RANSAC still supplies the quaternion plane and its inliers. Squaring the
    complex quaternion-plane projection removes the q/-q representation choice.
    The supplied bound is a search assumption, not proof against aliasing.
    Diagnostic coherence is an objective value, not a calibrated probability;
    signed frequency is relative to the fitted axis, while spin_rad_s is a vector.
    """

    def __init__(self, max_spin_rps):
        super().__init__()
        if not np.isfinite(max_spin_rps) or max_spin_rps <= 0:
            raise ValueError("max_spin_rps must be finite and positive")
        self.max_spin_rps = float(max_spin_rps)
        self.diagnostics = {}

    def RANSAC_regress(self, t, rots):
        t = np.asarray(t, dtype=float)
        self.diagnostics = {"search_bound_rps": self.max_spin_rps}
        if len(t) != len(rots):
            raise ValueError("Timestamps and poses must have equal lengths")
        if len(t) < self.n + self.d + 1:
            return None, None
        if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
            raise ValueError("Timestamps must be finite and strictly increasing")
        motion = np.array([(rot * rots[0].inv()).magnitude() for rot in rots])
        if np.max(motion) <= 1e-10:
            return None, None
        quats = np.array([self.rot2quat(rot) for rot in rots])
        axis, plane, indices, _ = self.get_RANSAC_spin_axis(quats)
        if axis is None:
            return None, None
        times = t[indices] - t[indices[0]]
        projection = quats[indices] @ plane[0] + 1j * (quats[indices] @ plane[1])
        z = projection**2 / np.abs(projection)**2
        bound = 2 * np.pi * self.max_spin_rps
        # At least 16 grid samples per full phase change over the observation.
        count = max(65, int(np.ceil(2 * bound * times[-1] * 8 / np.pi)) + 1)
        grid = np.linspace(-bound, bound, count)

        def coherence(omega):
            return float(abs(np.mean(z * np.exp(-1j * omega * times))))

        scores = np.abs(np.mean(z[:, None] * np.exp(-1j * times[:, None] * grid), axis=0))
        peaks = (np.flatnonzero((scores[1:-1] >= scores[:-2]) & (scores[1:-1] >= scores[2:])) + 1).tolist()
        if scores[0] >= scores[1]:
            peaks.append(0)
        if scores[-1] >= scores[-2]:
            peaks.append(len(grid) - 1)
        candidates = []
        for index in peaks:
            if index in (0, len(grid) - 1):
                omega = grid[index]
            else:
                fit = minimize_scalar(lambda w: -coherence(w), bounds=(grid[index-1], grid[index+1]), method="bounded", options={"xatol": 1e-10})
                omega = float(fit.x)
            candidates.append((coherence(omega), omega))
        candidates.sort(reverse=True)
        best_score, omega = candidates[0]
        residual = z * np.exp(-1j * omega * times)
        residual *= np.exp(-1j * np.angle(np.mean(residual)))
        phase_errors = np.degrees(np.angle(residual))
        secondary = None
        if len(candidates) > 1:
            second_score, second_omega = candidates[1]
            secondary = {"signed_frequency_rps": float(second_omega / (2*np.pi)), "spin_rad_s": (axis * second_omega).tolist(), "coherence": second_score, "coherence_gap": best_score - second_score}
        self.diagnostics.update({
            "signed_frequency_rps": float(omega / (2*np.pi)),
            "coherence": best_score,
            "secondary_peak": secondary,
            "at_search_boundary": bool(np.isclose(abs(omega), bound, rtol=1e-8)),
            "phase_residual_rms_deg": float(np.sqrt(np.mean(phase_errors**2))),
            "phase_residual_median_abs_deg": float(np.median(abs(phase_errors))),
            "grid_candidates": count,
        })
        return axis * omega, indices
