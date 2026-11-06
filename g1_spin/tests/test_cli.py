"""Headless integration checks validate delivered video frames and report data."""
import json
import cv2
import numpy as np
import pytest
from g1_spin.cli import analyze, preview
from g1_spin.events import EventBatch, save_npz


def config():
    return {'ball': {'center_x': 128, 'center_y': 128, 'radius_px': 80,
                     'inner_radius_ratio': .9},
            'analysis': {'start_us': 0, 'window_ms': 2, 'step_ms': 2,
                         'duration_ms': 6, 'max_rpm': 3500}}


def test_preview_produces_decodable_frames(tmp_path):
    source = tmp_path / 'events.npz'
    save_npz(source, EventBatch(np.array([100, 101, 102, 103]),
                               np.array([100, 101, 102, 103]),
                               np.array([0, 1000, 2000, 5000]), np.array([0, 1, 0, 1])))
    target = tmp_path / 'preview.mp4'
    assert preview(source, target, window_ms=2, start_us=0, duration_ms=6) == 3
    capture = cv2.VideoCapture(str(target))
    decoded = []
    try:
        while True:
            ok, frame = capture.read()
            if not ok:
                break
            decoded.append(frame)
    finally:
        capture.release()
    assert len(decoded) == 3
    assert all(frame.shape == (256, 256, 3) for frame in decoded)
    assert cv2.imread(str(target.with_suffix('.png'))) is not None
    assert (tmp_path / 'preview_surfaces.png').stat().st_size > 0
    metadata = json.loads(target.with_suffix('.json').read_text())
    assert metadata['frames'] == len(decoded)
    assert metadata['playback_speed'] == .1


def test_analyze_rejects_invalid_roi(tmp_path):
    settings = config()
    settings['ball']['radius_px'] = -1
    with pytest.raises(ValueError, match='radius_px'):
        analyze(tmp_path/'not_needed.npz', tmp_path/'out', settings)


def test_empty_analysis_reports_null_estimates(tmp_path):
    source = tmp_path / 'empty.npz'
    save_npz(source, EventBatch(np.array([]), np.array([]),
                               np.array([], dtype=np.int64), np.array([], dtype=np.uint8)))
    out = tmp_path/'analysis'
    result = analyze(source, out, config())
    assert result['window_count'] == 3
    assert result['valid_fraction'] == 0
    assert result['median_rpm'] is None
    for row in result['windows']:
        assert not row['valid']
        assert row['spin_rpm'] is None and row['axis_camera'] is None
        assert row['failure_reason'] == 'insufficient_events'
    persisted = json.loads((out/'result.json').read_text())
    assert persisted['windows'] == result['windows']
    for filename in ['windows.csv', 'report.md', 'rpm_curve.png',
                     'spherical_before_after.png', 'optimization.png']:
        assert (out/filename).stat().st_size > 0
    video = cv2.VideoCapture(str(out/'event_preview.mp4'))
    try:
        assert video.read()[0], 'Empty event recording must still render its explicit requested interval'
    finally:
        video.release()


def test_zero_frame_preview_rejects_missing_time_origin(tmp_path):
    source = tmp_path/'empty.npz'
    save_npz(source, EventBatch(np.array([]), np.array([]),
                               np.array([], dtype=np.int64), np.array([], dtype=np.uint8)))
    with pytest.raises(ValueError):
        preview(source, tmp_path/'no_frames.mp4')


def test_inspection_includes_silent_millisecond_bins(tmp_path):
    import csv
    from g1_spin.cli import inspect_file
    source = tmp_path/'gaps.npz'
    save_npz(source, EventBatch(np.array([100, 100]), np.array([100, 100]),
                               np.array([10000, 14000]), np.array([0, 1])))
    out = tmp_path/'inspect'
    inspect_file(source, out)
    with (out/'event_rate.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert [(int(row['relative_ms']), int(row['events'])) for row in rows] == [
        (0, 1), (1, 0), (2, 0), (3, 0), (4, 1)]
