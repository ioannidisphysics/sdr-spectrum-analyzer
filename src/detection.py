import numpy as np
from scipy.signal import find_peaks


def detect_spectrum_peaks(freqs, magnitude_db, threshold_db=-30, distance=20):
    """
    Detect peaks in an FFT spectrum.

    Parameters
    ----------
    freqs : np.ndarray
        Frequency axis in Hz.

    magnitude_db : np.ndarray
        Spectrum magnitude in dB.

    threshold_db : float
        Minimum peak height in dB.

    distance : int
        Minimum distance between detected peaks in FFT bins.

    Returns
    -------
    peak_freqs : np.ndarray
        Detected peak frequencies in Hz.

    peak_levels : np.ndarray
        Detected peak levels in dB.
    """

    peaks, properties = find_peaks(
        magnitude_db,
        height=threshold_db,
        distance=distance
    )

    peak_freqs = freqs[peaks]
    peak_levels = magnitude_db[peaks]

    return peak_freqs, peak_levels


def estimate_noise_floor(freqs, magnitude_db, peak_freqs, exclusion_bw=5_000):
    """
    Estimate the noise floor by excluding frequency regions around detected peaks.

    Parameters
    ----------
    freqs : np.ndarray
        Frequency axis in Hz.

    magnitude_db : np.ndarray
        Spectrum magnitude in dB.

    peak_freqs : np.ndarray
        Detected peak frequencies in Hz.

    exclusion_bw : float
        Frequency range around each peak to exclude from noise estimation.

    Returns
    -------
    noise_floor_db : float
        Estimated noise floor in dB.
    """

    noise_mask = np.ones_like(freqs, dtype=bool)

    for peak_freq in peak_freqs:
        noise_mask &= np.abs(freqs - peak_freq) > exclusion_bw

    noise_floor_db = np.median(magnitude_db[noise_mask])

    return noise_floor_db


def build_peak_report(peak_freqs, peak_levels, noise_floor_db):
    """
    Build a structured peak report with frequency, level, and SNR.

    Parameters
    ----------
    peak_freqs : np.ndarray
        Detected peak frequencies in Hz.

    peak_levels : np.ndarray
        Detected peak levels in dB.

    noise_floor_db : float
        Estimated noise floor in dB.

    Returns
    -------
    report : list[dict]
        List of detected peak measurements.
    """

    report = []

    for freq, level in zip(peak_freqs, peak_levels):
        snr_db = level - noise_floor_db

        report.append({
            "frequency_hz": freq,
            "frequency_khz": freq / 1e3,
            "level_db": level,
            "snr_db": snr_db,
        })

    return report