"""Known-motion tests: these have ground truth unlike the public image sample."""

import os
from pathlib import Path
import unittest

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".cache/matplotlib"))

import numpy as np
from scipy.spatial.transform import Rotation as R

from corrected_spin import BoundedSpinRegressor, CorrectedSpinRegressor
from upstream.python.spin_regressor import SpinRegressor


class CorrectedSpinTests(unittest.TestCase):
    def setUp(self):
        self.times = np.r_[0.0, np.cumsum(np.tile([0.0008, 0.0013, 0.0010], 12))]
        self.initial = R.from_euler("xyz", [0.3, -0.6, 0.8])

    def poses(self, omega):
        return R.from_rotvec(self.times[:, None] * omega) * self.initial

    def test_reproduces_upstream_sign_mismatch(self):
        omega = np.array([20.0, 30.0, -40.0])
        np.random.seed(0)
        raw, _ = SpinRegressor().RANSAC_regress(self.times, self.poses(omega))
        np.testing.assert_allclose(raw, -omega, atol=1e-7)

    def test_axes_directions_and_forward_motion(self):
        cases = [
            [70.0, 0.0, 0.0], [-70.0, 0.0, 0.0],
            [0.0, 110.0, 0.0], [0.0, -110.0, 0.0],
            [0.0, 0.0, 170.0], [0.0, 0.0, -170.0],
            [20.0, 30.0, -40.0], [-20.0, -30.0, 40.0],
        ]
        for values in cases:
            with self.subTest(omega=values):
                omega = np.array(values)
                poses = self.poses(omega)
                np.random.seed(0)
                estimate, indices = CorrectedSpinRegressor().RANSAC_regress(self.times, poses)
                self.assertIsNotNone(estimate)
                np.testing.assert_allclose(estimate, omega, atol=1e-7)
                predicted = R.from_rotvec(self.times[:, None] * estimate) * self.initial
                self.assertLess(np.max((predicted * poses.inv()).magnitude()), 1e-8)
                self.assertEqual(len(indices), len(self.times))

    def test_constant_poses_are_rejected_not_called_zero_spin(self):
        poses = self.poses(np.zeros(3))
        self.assertEqual(CorrectedSpinRegressor().RANSAC_regress(self.times, poses), (None, None))

    def test_insufficient_samples_rejected(self):
        poses = self.poses(np.array([20.0, 30.0, -40.0]))
        for count in range(7):
            with self.subTest(count=count):
                self.assertEqual(CorrectedSpinRegressor().RANSAC_regress(self.times[:count], poses[:count]), (None, None))

    def test_invalid_timestamps_rejected(self):
        poses = self.poses(np.array([20.0, 30.0, -40.0]))
        repeated = self.times.copy()
        repeated[3] = repeated[2]
        nonfinite = self.times.copy()
        nonfinite[3] = np.nan
        for times in [self.times[::-1], repeated, nonfinite]:
            with self.assertRaises(ValueError):
                CorrectedSpinRegressor().RANSAC_regress(times, poses)
        with self.assertRaises(ValueError):
            CorrectedSpinRegressor().RANSAC_regress(self.times[:-1], poses)


class BoundedSpinTests(unittest.TestCase):
    def setUp(self):
        # Missing frames make individual jumps exceed pi at 133 rps.
        self.times = np.arange(27) / 381.0
        self.times = np.delete(self.times, [8, 15, 18, 21, 24])
        self.axis = np.array([0.2, -0.95, 0.24])
        self.axis /= np.linalg.norm(self.axis)
        self.initial = R.from_euler("xyz", [0.3, -0.6, 0.8])

    def check_sequence(self, omega, times, poses):
        regressor = BoundedSpinRegressor(max_spin_rps=190.5)
        np.random.seed(0)
        estimate, indices = regressor.RANSAC_regress(times, poses)
        np.testing.assert_allclose(estimate, omega, atol=2e-4)
        self.assertEqual(len(indices), len(times))
        self.assertLess(regressor.diagnostics["phase_residual_rms_deg"], 1e-4)
        self.assertIn("secondary_peak", regressor.diagnostics)
        self.assertFalse(regressor.diagnostics["at_search_boundary"])
        return estimate

    def test_dropped_frames_both_directions(self):
        for sign in (1, -1):
            with self.subTest(sign=sign):
                omega = self.axis * sign * 133 * 2*np.pi
                poses = R.from_rotvec(self.times[:, None] * omega) * self.initial
                self.check_sequence(omega, self.times, poses)

    def test_quaternion_sign_noise_and_time_origin(self):
        omega = self.axis * 133 * 2*np.pi
        poses = R.from_rotvec(self.times[:, None] * omega) * self.initial
        quats = poses.as_quat()
        quats[::2] *= -1
        noise = np.random.default_rng(7).normal(scale=1e-11, size=(len(poses), 3))
        perturbed = R.from_rotvec(noise) * R.from_quat(quats)
        self.check_sequence(omega, self.times + 10_000.0, perturbed)

    def test_low_speed(self):
        for magnitude in (2.0, -2.0, 20.0, -20.0):
            omega = self.axis * magnitude
            poses = R.from_rotvec(self.times[:, None] * omega) * self.initial
            self.check_sequence(omega, self.times, poses)

    def test_invalid_bound_and_degenerate_sequences(self):
        for bound in (0, -1, np.inf, np.nan):
            with self.assertRaises(ValueError):
                BoundedSpinRegressor(bound)
        regressor = BoundedSpinRegressor(190.5)
        static = R.from_rotvec(np.zeros((len(self.times), 3))) * self.initial
        self.assertEqual(regressor.RANSAC_regress(self.times, static), (None, None))
        self.assertEqual(regressor.RANSAC_regress(self.times[:4], static[:4]), (None, None))


if __name__ == "__main__":
    unittest.main(verbosity=2)
