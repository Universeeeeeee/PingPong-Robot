import os
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))
sys.path.insert(0, str(ROOT / "upstream/python"))
from baygeohasher import BayGeoHasher
from point_matcher import match_points


class PointMatcherTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = np.loadtxt(ROOT / "upstream/cad/ref_pattern.csv", delimiter=",")
        cls.hasher = BayGeoHasher()
        cls.hasher.gen_hash_table(cls.reference)

    def test_known_pose_and_permutation(self):
        expected = Rotation.from_rotvec([0.3, -0.5, 0.2])
        visible = expected.apply(self.reference)
        visible = visible[visible[:, 2] > 0.15]
        baseline = match_points(visible, self.hasher)
        self.assertLess((baseline["rotation"] * expected.inv()).magnitude(), 1e-6)
        rng = np.random.default_rng(123)
        for _ in range(3):
            result = match_points(rng.permutation(visible), self.hasher)
            self.assertLess((result["rotation"] * baseline["rotation"].inv()).magnitude(), 1e-9)
            self.assertAlmostEqual(result["rmse"], baseline["rmse"], places=12)

    def test_too_few_or_invalid_features(self):
        self.assertIsNone(match_points(self.reference[:3], self.hasher))
        with self.assertRaises(ValueError):
            match_points(np.zeros((4, 2)), self.hasher)
        with self.assertRaises(ValueError):
            match_points(np.full((4, 3), np.nan), self.hasher)


if __name__ == "__main__":
    unittest.main()
