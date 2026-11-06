"""Idealized moving surface landmarks, NOT a physical contrast-event simulator.

Only x/y/time/polarity leave this generator. Landmark identities and ground truth
are never included in estimator inputs. Subpixel locations deliberately isolate
geometric correctness from sensor quantization and texture/event physics.
"""
import numpy as np
from scipy.spatial.transform import Rotation
from .events import EventBatch


def generate(rpm=1000, axis=(.3, .7, .6), duration_s=.06, seed=0,
             noise_fraction=0, center=(128, 128), radius=100, sample_hz=8400):
    axis = np.asarray(axis, dtype=float)
    if axis.shape != (3,) or not np.all(np.isfinite(axis)) or np.linalg.norm(axis) == 0:
        raise ValueError('axis must be a finite nonzero 3-vector')
    if not np.isfinite(rpm) or duration_s <= 0 or radius <= 0 or sample_hz <= 0 or not 0 <= noise_fraction < 1:
        raise ValueError('Invalid synthetic parameters')
    if len(center) != 2 or min(center) < radius:
        raise ValueError('Ball must lie inside the nonnegative image plane')
    if rpm == 0:
        return EventBatch(np.array([]), np.array([]), np.array([], dtype=np.int64), np.array([], dtype=np.uint8))
    rng = np.random.default_rng(seed)
    landmarks = rng.normal(size=(24, 3))
    landmarks /= np.linalg.norm(landmarks, axis=1, keepdims=True)
    omega = axis / np.linalg.norm(axis) * rpm * (2 * np.pi / 60)
    times = np.linspace(0, duration_s, max(2, int(duration_s * sample_hz)), endpoint=False)
    rotations = Rotation.from_rotvec(times[:, None] * omega).as_matrix()
    points = np.einsum('tij,kj->tki', rotations, landmarks)
    visible = points[:, :, 2] < -0.12
    x = center[0] + radius * points[:, :, 0][visible]
    y = center[1] + radius * points[:, :, 1][visible]
    t_us = np.broadcast_to(np.rint(times * 1e6).astype(np.int64)[:, None], visible.shape)[visible]
    polarity = np.broadcast_to((np.arange(24) % 2).astype(np.uint8), visible.shape)[visible]
    n_noise = int(len(x) * noise_fraction / max(1 - noise_fraction, 1e-6)) if noise_fraction < 1 else len(x)
    if n_noise:
        x = np.concatenate([x, rng.uniform(0, 2 * center[0], n_noise)])
        y = np.concatenate([y, rng.uniform(0, 2 * center[1], n_noise)])
        t_us = np.concatenate([t_us, rng.integers(0, max(1, int(duration_s * 1e6)), n_noise)])
        polarity = np.concatenate([polarity, rng.integers(0, 2, n_noise, dtype=np.uint8)])
    order = np.argsort(t_us, kind='stable')
    return EventBatch(x[order], y[order], t_us[order], polarity[order])
