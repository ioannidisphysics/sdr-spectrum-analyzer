"""
Power spectral density estimation for the SDR Spectrum Analyzer.

The existing compute_fft_spectrum() in spectrum.py normalises by the coherent
gain of the window, which is the right choice for reading the level of a tone.
A noise measurement needs the opposite normalisation: power per unit bandwidth,
using the sum of the squared window. Welch's method is used so the noise floor
is an average over many segments instead of a single noisy realisation.
"""

import numpy as np

WINDOWS = {
    "hann": np.hanning,
    "hamming": np.hamming,
    "blackman": np.blackman,
    "rectangular": np.ones,
}


def window_metrics(window, fs):
    """
    Noise bandwidth of a window.

    The equivalent noise bandwidth (ENBW) is the width of the ideal rectangular
    filter that would pass the same noise power as this window. It is what the
    resolution bandwidth of the analyser actually is, and it is wider than one
    FFT bin for every window except the rectangular one: 1.5 bins for Hann.

    Parameters
    ----------
    window : np.ndarray
        Window coefficients.

    fs : float
        Sampling frequency in Hz.

    Returns
    -------
    enbw_bins : float
        Equivalent noise bandwidth in FFT bins.

    rbw_hz : float
        Resolution bandwidth in Hz.
    """

    s1 = np.sum(window)
    s2 = np.sum(window ** 2)

    enbw_bins = len(window) * s2 / s1 ** 2
    rbw_hz = enbw_bins * fs / len(window)

    return enbw_bins, rbw_hz


def welch_psd(x, fs, nfft=4096, overlap=0.5, window_name="hann", remove_dc=True):
    """
    Welch power spectral density of complex IQ samples.

    Parameters
    ----------
    x : np.ndarray
        Complex IQ samples, full scale roughly +/-1.

    fs : float
        Sampling frequency in Hz.

    nfft : int
        FFT size of each segment. Sets the resolution bandwidth together with
        the window: RBW = ENBW * fs / nfft.

    overlap : float
        Fractional overlap between segments, 0 to 0.9. Hann segments are
        normally overlapped by 0.5 so no samples are wasted under the taper.

    window_name : str
        Key into WINDOWS.

    remove_dc : bool
        Subtract the mean of the capture first. The RTL-SDR has a DC offset at
        the centre of the band that is an artefact of the receiver, not a
        signal. Removing it keeps it from dominating the centre bins; the
        residual is excluded by band edge trimming in sweep.py anyway.

    Returns
    -------
    freqs : np.ndarray
        Baseband frequency axis from -fs/2 to +fs/2, in Hz.

    psd_db : np.ndarray
        Power spectral density in dBFS per Hz.

    info : dict
        Segment count, ENBW in bins and RBW in Hz.
    """

    if window_name not in WINDOWS:
        raise ValueError(f"unknown window {window_name!r}, expected one of {sorted(WINDOWS)}")

    if not 0.0 <= overlap < 0.95:
        raise ValueError("overlap must be in [0, 0.95)")

    x = np.asarray(x)

    if len(x) < nfft:
        raise ValueError(f"capture has {len(x)} samples, needs at least nfft = {nfft}")

    if remove_dc:
        x = x - np.mean(x)

    window = WINDOWS[window_name](nfft)
    s2 = np.sum(window ** 2)
    step = max(1, int(round(nfft * (1.0 - overlap))))

    # -----------------------------
    # Average the periodograms
    # -----------------------------

    accumulator = np.zeros(nfft)
    segments = 0

    for start in range(0, len(x) - nfft + 1, step):
        segment = x[start:start + nfft] * window
        spectrum = np.fft.fftshift(np.fft.fft(segment))
        accumulator += np.abs(spectrum) ** 2
        segments += 1

    # Divide by fs * sum(w^2) so the result is power per Hz and does not change
    # when nfft changes. That is what makes two sweeps with different settings
    # comparable.
    psd = accumulator / (segments * fs * s2)

    freqs = np.fft.fftshift(np.fft.fftfreq(nfft, d=1 / fs))
    psd_db = 10 * np.log10(psd + 1e-20)

    enbw_bins, rbw_hz = window_metrics(window, fs)

    info = {
        "segments": segments,
        "enbw_bins": enbw_bins,
        "rbw_hz": rbw_hz,
        "window": window_name,
        "nfft": nfft,
    }

    return freqs, psd_db, info


def moving_average_db(values, bins):
    """
    Smooth a dB trace with a centred moving average.

    Averaging is done on the dB values, which is a median-like smoothing of the
    trace rather than a true power average. That is what a spectrum analyser's
    video filter does, and it is the right choice here because the aim is to
    see the shape of the noise floor, not to preserve absolute peak levels.

    Parameters
    ----------
    values : np.ndarray
        Trace in dB.

    bins : int
        Window length in bins. Values below 2 return the trace unchanged.

    Returns
    -------
    np.ndarray
        Smoothed trace, same length as the input.
    """

    if bins < 2:
        return values

    kernel = np.ones(bins) / bins
    padded = np.pad(values, (bins // 2, bins - 1 - bins // 2), mode="edge")

    return np.convolve(padded, kernel, mode="valid")
