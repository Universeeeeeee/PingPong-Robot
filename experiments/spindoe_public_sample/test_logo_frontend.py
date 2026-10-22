"""End-to-end tests on rendered text on a 3D sphere, not real captured data."""
import json
import unittest

import cv2
import numpy as np
from scipy.spatial.transform import Rotation

from logo_frontend import estimate_sequence


def render_text_ball(rotation, circle=(96,96,72), texture=True):
    image = np.full((192,192,3), 25, np.uint8)
    yy, xx = np.mgrid[:192,:192]
    x, y = (xx-circle[0])/circle[2], -(yy-circle[1])/circle[2]
    inside = x*x+y*y <= 1
    points = np.c_[x[inside], y[inside], np.sqrt(1-x[inside]**2-y[inside]**2)]
    shade = 215+20*points[:,2]
    if texture:
        ink = np.full((256,256),255,np.uint8)
        cv2.putText(ink,"R7",(37,158),cv2.FONT_HERSHEY_SIMPLEX,2.6,0,6,cv2.LINE_AA)
        cv2.line(ink,(53,180),(172,174),0,4,cv2.LINE_AA)
        reference = rotation.inv().apply(points)
        u = ((reference[:,0]+0.6)/1.2*255).astype(np.float32)
        v = ((0.6-reference[:,1])/1.2*255).astype(np.float32)
        values = cv2.remap(ink,u[None],v[None],cv2.INTER_LINEAR,borderMode=cv2.BORDER_CONSTANT,borderValue=255)[0]
        mark = (reference[:,2]>0.2) & (values<255)
        shade[mark] *= 0.15+0.85*values[mark]/255
    image[inside] = shade[:,None].astype(np.uint8)
    return image


class LogoFrontendTests(unittest.TestCase):
    def sequence(self, omega, moving=False, texture=True):
        times = np.arange(12)*0.008
        circles = np.tile([96.,96.,72.],(len(times),1))
        if moving:
            circles[:,0] += np.linspace(-7,7,len(times))
            circles[:,1] += np.linspace(4,-4,len(times))
            circles[:,2] += np.linspace(-4,4,len(times))
        frames = [render_text_ball(Rotation.from_rotvec(omega*t),c,texture) for t,c in zip(times,circles)]
        return frames,times,circles

    def test_rendered_text_known_rotation_both_directions(self):
        for values in [[1.4,-2.2,3.0],[-1.4,2.2,-3.0],[3.,0.,0.],[0.,3.,0.],[0.,0.,3.]]:
            truth = np.array(values)
            result = estimate_sequence(*self.sequence(truth))
            estimates = np.array([w["omega_rad_s"] for w in result["windows"] if w["valid"]])
            self.assertGreaterEqual(len(estimates),8)
            self.assertLess(np.linalg.norm(np.median(estimates,axis=0)-truth),0.5)
            json.dumps(result,allow_nan=False)

    def test_disappearing_print_ends_track_identity(self):
        frames,times,circles=self.sequence(np.array([1.4,-2.2,3.]))
        for i in range(4,8):
            frames[i]=render_text_ball(Rotation.identity(),circles[i],texture=False)
        result=estimate_sequence(frames,times,circles)
        self.assertTrue(all(not w["valid"] for w in result["windows"][4:7]))
        before={p["track_id"] for p in result["tracks"] if p["frame"] <= 3}
        after={p["track_id"] for p in result["tracks"] if p["frame"] >= 9}
        self.assertTrue(before.isdisjoint(after))

    def test_untrackable_fast_rotation_rejected(self):
        result=estimate_sequence(*self.sequence(np.array([0.,0.,100.])))
        self.assertEqual(result["valid_window_count"],0)

    def test_translation_and_scale_compensation(self):
        result = estimate_sequence(*self.sequence(np.zeros(3),moving=True))
        estimates = np.array([w["omega_rad_s"] for w in result["windows"] if w["valid"]])
        self.assertGreaterEqual(len(estimates),8)
        # Rasterization and per-frame resampling prevent exact zero at this size.
        self.assertLess(np.linalg.norm(np.median(estimates,axis=0)),0.5)

    def test_blank_ball_is_unknown_not_zero(self):
        result = estimate_sequence(*self.sequence(np.zeros(3),texture=False))
        self.assertEqual(result["valid_window_count"],0)
        self.assertTrue(all(w["omega_rad_s"] is None for w in result["windows"]))

    def test_invalid_capture_times(self):
        frames,times,circles=self.sequence(np.array([1.,2.,3.]))
        for invalid in [np.zeros_like(times),times[::-1],np.full_like(times,np.nan)]:
            with self.assertRaises(ValueError):
                estimate_sequence(frames,invalid,circles)


if __name__=="__main__":
    unittest.main(verbosity=2)
