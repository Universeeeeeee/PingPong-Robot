"""Small, validated event container; time differences precede float conversion."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass
class EventBatch:
    x: np.ndarray
    y: np.ndarray
    t_us: np.ndarray
    polarity: np.ndarray

    def __post_init__(self):
        self.x = np.asarray(self.x)
        self.y = np.asarray(self.y)
        self.t_us = np.asarray(self.t_us)
        self.polarity = np.asarray(self.polarity)
        self.validate()
        self.t_us = self.t_us.astype(np.int64, copy=False)

    def __len__(self):
        return len(self.t_us)

    def seconds(self, reference_us=None):
        if not len(self):
            return np.empty(0, dtype=float)
        reference_us = int(self.t_us[0]) if reference_us is None else int(reference_us)
        return (self.t_us - reference_us).astype(float) * 1e-6

    def select(self, mask):
        return EventBatch(self.x[mask], self.y[mask], self.t_us[mask], self.polarity[mask])

    def validate(self, width=None, height=None):
        arrays = (self.x, self.y, self.t_us, self.polarity)
        if any(a.ndim != 1 for a in arrays) or len({len(a) for a in arrays}) != 1:
            raise ValueError('Event fields must be equal-length one-dimensional arrays')
        if any(not np.all(np.isfinite(a)) for a in arrays):
            raise ValueError('Event fields must be finite')
        if len(self.t_us) and (not np.issubdtype(self.t_us.dtype, np.integer) or
                              np.any(self.t_us < 0) or np.any(self.t_us >= np.iinfo(np.int64).max)):
            raise ValueError('Timestamps must be nonnegative int64 microseconds')
        if np.any(np.diff(self.t_us.astype(np.int64)) < 0):
            raise ValueError('Timestamps must be sorted')
        if not np.all(np.isin(self.polarity, [0, 1])):
            raise ValueError('Polarity must be 0/OFF or 1/ON')
        if np.any(self.x < 0) or np.any(self.y < 0):
            raise ValueError('Event coordinates must be nonnegative')
        if width is not None and np.any(self.x >= width):
            raise ValueError('Event x exceeds sensor width')
        if height is not None and np.any(self.y >= height):
            raise ValueError('Event y exceeds sensor height')
        gaps = np.diff(self.t_us.astype(np.int64))
        return {'sorted': True, 'duplicate_timestamps': int(np.count_nonzero(gaps == 0)),
                'max_gap_us': int(gaps.max()) if gaps.size else 0,
                'note': 'Equal timestamps are normal for event sensors; gaps require recording context.'}


def save_npz(path, batch):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, x=batch.x, y=batch.y, t_us=batch.t_us, polarity=batch.polarity)


def load_npz(path):
    with np.load(path, allow_pickle=False) as data:
        return EventBatch(*(data[key] for key in ('x', 'y', 't_us', 'polarity')))
