import numpy as np


def compute_fft_spectrum(x, fs):
    """
    Compute FFT spectrum of complex IQ samples using a Hann window.

    Parameters
    ----------
    x : np.ndarray
        Complex IQ samples.

    fs : float
        Sampling frequency in Hz.

    Returns
    -------
    freqs_shifted : np.ndarray
        Frequency axis from -fs/2 to +fs/2.

    magnitude_db : np.ndarray
        Normalized FFT magnitude in dB.
    """

    N = len(x)

    window = np.hanning(N)
    x_windowed = x * window

    X = np.fft.fft(x_windowed)
    freqs = np.fft.fftfreq(N, d=1 / fs)

    X_shifted = np.fft.fftshift(X)
    freqs_shifted = np.fft.fftshift(freqs)

    magnitude = np.abs(X_shifted) / np.sum(window)
    magnitude_db = 20 * np.log10(magnitude + 1e-12)

    return freqs_shifted, magnitude_db


def compute_spectrogram(x, fs, nfft=1024, hop=512):
    """
    Compute a simple STFT spectrogram.

    Parameters
    ----------
    x : np.ndarray
        Complex IQ samples.

    fs : float
        Sampling frequency in Hz.

    nfft : int
        FFT size for each frame.

    hop : int
        Hop size between frames.

    Returns
    -------
    spectrogram : np.ndarray
        Spectrogram matrix in dB.

    spectrogram_freqs : np.ndarray
        Frequency axis in Hz.

    spectrogram_times : np.ndarray
        Time axis in seconds.
    """

    stft_window = np.hanning(nfft)
    spectrogram_frames = []

    for start in range(0, len(x) - nfft, hop):
        frame = x[start:start + nfft]
        frame_windowed = frame * stft_window

        X_frame = np.fft.fft(frame_windowed)
        X_frame_shifted = np.fft.fftshift(X_frame)

        frame_magnitude = np.abs(X_frame_shifted) / np.sum(stft_window)
        frame_magnitude_db = 20 * np.log10(frame_magnitude + 1e-12)

        spectrogram_frames.append(frame_magnitude_db)

    spectrogram = np.array(spectrogram_frames).T

    spectrogram_freqs = np.fft.fftshift(np.fft.fftfreq(nfft, d=1 / fs))
    spectrogram_times = np.arange(spectrogram.shape[1]) * hop / fs

    return spectrogram, spectrogram_freqs, spectrogram_times    