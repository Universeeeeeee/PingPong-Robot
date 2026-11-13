"""Offline image-plane translation only; radius and event identities stay fixed."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


def integer_times(values, name):
    values = np.asarray(values)
    if values.ndim != 1 or not np.issubdtype(values.dtype, np.integer):
        raise ValueError(f"{name} must be a vector of integer microseconds")
    # Python integer differences avoid uint wrap and signed integer overflow.
    result = [int(value) for value in values]
    if any(value < 0 for value in result):
        raise ValueError(f"{name} must be nonnegative")
    return result


def validate_metadata(metadata, event_timebase):
    if metadata.get("image_plane") != "G1 image plane":
        raise ValueError("Centers must be in G1 image plane pixel coordinates")
    if metadata.get("center_semantics") != "ball_center":
        raise ValueError("Trajectory must describe ball center, not marker or pivot")
    for key in ("timebase", "source", "synchronization_basis", "center_error"):
        if not isinstance(metadata.get(key), str) or not metadata[key].strip():
            raise ValueError(f"Missing nonempty metadata field: {key}")
    if metadata["timebase"] != event_timebase:
        raise ValueError("Center and event timebase declarations must match")
    if metadata.get("synchronization_verified") is not True:
        raise ValueError("Explicit verified synchronization declaration required")


def centers_at(sample_t_us, sample_xy, query_t_us, *, max_gap_us):
    """Piecewise linear interpolation; no extrapolation or oversized intervals."""
    samples = integer_times(sample_t_us, "sample_t_us")
    queries = integer_times(query_t_us, "query_t_us")
    xy = np.asarray(sample_xy, dtype=float)
    if len(samples) < 2 or xy.shape != (len(samples), 2) or not np.isfinite(xy).all():
        raise ValueError("At least two finite center samples are required")
    if not isinstance(max_gap_us, (int, np.integer)) or max_gap_us <= 0:
        raise ValueError("max_gap_us must be an explicit positive integer")
    gaps = [right-left for left, right in zip(samples, samples[1:])]
    if any(gap <= 0 for gap in gaps):
        raise ValueError("Center sample times must be strictly increasing")
    if any(gap > max_gap_us for gap in gaps):
        raise ValueError("Center sample gap exceeds configured max_gap_us")
    if any(t < samples[0] or t > samples[-1] for t in queries):
        raise ValueError("Center samples do not cover all requested times")
    if samples[-1]-samples[0] > 2**53:
        raise ValueError("Relative microsecond span exceeds exact float integer range")
    origin = samples[0]
    sample_relative = np.array([t-origin for t in samples], dtype=float)
    query_relative = np.array([t-origin for t in queries], dtype=float)
    return np.column_stack([np.interp(query_relative, sample_relative, xy[:, dim])
                            for dim in range(2)])


def compensate(xy, event_t_us, sample_t_us, sample_xy, *, reference_us,
               metadata, event_timebase, max_gap_us):
    """Return u'=u-c(t)+c(tref), and c(tref), in float pixels.

    Declarations are a guard against silent mixing, not independent evidence
    that an externally supplied trajectory is actually synchronized or correct.
    """
    validate_metadata(metadata, event_timebase)
    times = integer_times(event_t_us, "event_t_us")
    if any(right < left for left, right in zip(times, times[1:])):
        raise ValueError("Event times must be nondecreasing")
    if (isinstance(reference_us, (bool, np.bool_)) or
            not isinstance(reference_us, (int, np.integer)) or reference_us < 0):
        raise ValueError("reference_us must be integer microseconds")
    coords = np.asarray(xy, dtype=float)
    if coords.shape != (len(times), 2) or not np.isfinite(coords).all():
        raise ValueError("Event coordinates must be finite Nx2 pixels")
    queries = np.asarray(times + [int(reference_us)], dtype=np.uint64)
    centers = centers_at(sample_t_us, sample_xy, queries, max_gap_us=max_gap_us)
    return coords-(centers[:-1]-centers[-1]), centers[-1]


def load_centers(csv_path, metadata_path):
    """CSV columns: time_us,cx_px,cy_px; JSON carries synchronization evidence."""
    with Path(csv_path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["time_us", "cx_px", "cy_px"]:
            raise ValueError("Expected CSV header time_us,cx_px,cy_px")
        rows = list(reader)
    parsed_times = [int(row["time_us"]) for row in rows]
    if any(t < 0 for t in parsed_times):
        raise ValueError("Center times must be nonnegative")
    times = np.array(parsed_times, dtype=np.uint64)
    xy = np.array([[float(row["cx_px"]), float(row["cy_px"])] for row in rows])
    return times, xy, json.loads(Path(metadata_path).read_text())
