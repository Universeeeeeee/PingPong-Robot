"""Bridge to the installed, validated SDK reader; no EVT3 binary decoder here."""
import json
import importlib.util
import os
from pathlib import Path
import numpy as np
from .events import EventBatch, load_npz


def _empty():
    return EventBatch(np.empty(0), np.empty(0), np.empty(0, dtype=np.int64), np.empty(0, dtype=np.uint8))


def source(path, batch_size=250_000):
    """Return metadata and a single-pass iterator; RAW uses existing SDK adapter.

    NPZ is the small demo interchange format and is loaded in memory. RAW stays
    streaming. Event extent in NPZ is not a measured sensor resolution.
    """
    if not isinstance(batch_size, (int, np.integer)) or batch_size <= 0:
        raise ValueError('batch_size must be a positive integer')
    path = Path(path)
    if path.suffix.lower() == '.npz':
        batch = load_npz(path)
        info = {'source': str(path.resolve()), 'timestamp_unit': 'us',
                'width': int(np.floor(batch.x.max())) + 1 if len(batch) else None,
                'height': int(np.floor(batch.y.max())) + 1 if len(batch) else None,
                'resolution_source': 'event_extent_not_sensor_resolution',
                'format': 'g1_spin EventBatch NPZ', 'synthetic': None}
        return info, (batch.select(slice(i, i + batch_size)) for i in range(0, len(batch), batch_size))
    manifest = path / 'manifest.json' if path.is_dir() else path
    if manifest.suffix.lower() == '.json':
        info = json.loads(manifest.read_text())
        def decoded_batches():
            for part in info['parts']:
                events = np.load(manifest.parent / part, mmap_mode='r', allow_pickle=False)
                for i in range(0, len(events), batch_size):
                    chunk = events[i:i + batch_size]
                    yield EventBatch(chunk['x'], chunk['y'], chunk['t'], chunk['p'])
        return info, decoded_batches()
    root = Path(os.environ.get('G1_LEGACY_ROOT', str(Path.home() / 'g1_spin'))).expanduser()
    adapter = root / 'scripts' / 'event_io.py'
    if not adapter.is_file():
        raise RuntimeError(f'Validated official SDK adapter not found: {adapter}; set G1_LEGACY_ROOT')
    spec = importlib.util.spec_from_file_location('_g1_legacy_event_io', adapter)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    info, batches = module.source(path, batch_size=batch_size)

    def converted():
        for events in batches:
            yield EventBatch(events['x'], events['y'], events['t'], events['p'])
    return info, converted()


def windows(batches, window_us, start_us=None, duration_us=None, step_us=None):
    """Yield half-open (start,end,EventBatch) windows, including empty gaps.

    Times are absolute integer microseconds. Without an explicit duration,
    iteration stops at the last observed event. An explicit start and duration
    also allow a completely empty recording to produce empty windows.
    """
    step_us = window_us if step_us is None else step_us
    for name, value in [('window_us', window_us), ('step_us', step_us), ('duration_us', duration_us)]:
        if value is not None and (not isinstance(value, (int, np.integer)) or value <= 0):
            raise ValueError(f'{name} must be a positive integer')
    if start_us is not None and (not isinstance(start_us, (int, np.integer)) or start_us < 0):
        raise ValueError('start_us must be a nonnegative integer')
    start = int(start_us) if start_us is not None else None
    limit = start + duration_us if start is not None and duration_us is not None else None
    buffer = _empty()
    last = None
    for batch in batches:
        batch.validate()
        if not len(batch):
            continue
        if last is not None and int(batch.t_us[0]) < last:
            raise ValueError('Timestamp regression across event batches')
        last = int(batch.t_us[-1])
        if start is None:
            start = int(batch.t_us[0])
            limit = start + duration_us if duration_us is not None else None
        if limit is not None and start >= limit:
            return
        incoming = batch.select(slice(int(np.searchsorted(batch.t_us, start)), None))
        if len(incoming):
            buffer = EventBatch(*(np.concatenate([getattr(buffer, key), getattr(incoming, key)])
                                  for key in ('x', 'y', 't_us', 'polarity')))
        while start <= last:
            end = min(start + window_us, limit) if limit is not None else start + window_us
            if last < end:
                break
            right = int(np.searchsorted(buffer.t_us, end))
            yield start, end, buffer.select(slice(0, right))
            start += step_us
            if limit is not None and start >= limit:
                return
            buffer = buffer.select(slice(int(np.searchsorted(buffer.t_us, start)), None))
    if start is None:
        return
    while (start < limit if limit is not None else last is not None and start <= last):
        end = min(start + window_us, limit) if limit is not None else start + window_us
        right = int(np.searchsorted(buffer.t_us, end))
        yield start, end, buffer.select(slice(0, right))
        start += step_us
        buffer = buffer.select(slice(int(np.searchsorted(buffer.t_us, start)), None))
