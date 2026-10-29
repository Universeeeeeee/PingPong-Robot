"""Conditional event-rate periodicity baseline; cannot identify spin axis."""
import numpy as np
from scipy.signal import find_peaks


def estimate_periodicity(events, bin_us=500):
    result = {'valid': False, 'method': 'event_rate_periodicity', 'candidate_rpm': None,
              'candidate_frequency_hz': None, 'period_strength': None, 'candidates_hz': [],
              'ambiguity': 'Event repetition may be a harmonic, illumination flicker, or unrelated to spin.',
              'failure_reason': None}
    if bin_us <= 0:
        raise ValueError('bin_us must be positive')
    if len(events) < 50:
        result['failure_reason'] = 'insufficient_events'
        return result
    relative = events.t_us - events.t_us[0]
    n = int(relative[-1] // bin_us) + 1
    if n < 32 or n > 2_000_000:
        result['failure_reason'] = 'unsupported_duration_or_bin_size'
        return result
    counts = np.bincount((relative // bin_us).astype(int), minlength=n).astype(float)
    signal = counts - counts.mean()
    if np.std(signal) < 1e-9:
        result['failure_reason'] = 'no_rate_modulation'
        return result
    spectrum = np.abs(np.fft.rfft(signal * np.hanning(n))) ** 2
    spectrum[:3] = 0  # Require at least three periods in the observation.
    peaks, _ = find_peaks(spectrum)
    if not len(peaks):
        result['failure_reason'] = 'no_periodicity'
        return result
    order = peaks[np.argsort(spectrum[peaks])[::-1]]
    k = int(order[0])
    strength = float(spectrum[max(0, k-1):k+2].sum() / max(spectrum.sum(), 1e-12))
    frequency = k / (n * bin_us * 1e-6)
    # A peak must also recur at a corresponding time lag.
    lag = int(round(n / k))
    corr = float(np.dot(signal[:-lag], signal[lag:]) / max(np.dot(signal, signal), 1e-12)) if 0 < lag < n else 0
    result.update(period_strength=strength, autocorrelation=corr,
                  candidates_hz=[float(j / (n * bin_us * 1e-6)) for j in order[:5]])
    if strength < .35 or corr < .25:
        result['failure_reason'] = 'weak_or_ambiguous_periodicity'
        return result
    result.update(valid=True, candidate_frequency_hz=frequency, candidate_rpm=frequency * 60)
    return result
