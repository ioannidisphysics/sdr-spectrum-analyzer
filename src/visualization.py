import matplotlib.pyplot as plt
import numpy as np


def plot_time_domain(t, x, save_path=None):
    """
    Plot I and Q components of the complex IQ signal.

    Parameters
    ----------
    t : np.ndarray
        Time vector in seconds.

    x : np.ndarray
        Complex IQ samples.

    save_path : Path or None
        If provided, save the figure to this path.
    """

    plt.figure(figsize=(10, 4))
    plt.plot(t[:300] * 1e6, np.real(x[:300]), label="I / Real")
    plt.plot(t[:300] * 1e6, np.imag(x[:300]), label="Q / Imaginary")
    plt.xlabel("Time [μs]")
    plt.ylabel("Amplitude")
    plt.title("Synthetic IQ Signal with Tones, AM and Chirp")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=200)

    plt.show()


def plot_fft_spectrum(freqs, magnitude_db, peak_freqs, peak_levels, save_path=None):
    """
    Plot FFT spectrum with detected peaks.

    Parameters
    ----------
    freqs : np.ndarray
        Frequency axis in Hz.

    magnitude_db : np.ndarray
        Spectrum magnitude in dB.

    peak_freqs : np.ndarray
        Detected peak frequencies in Hz.

    peak_levels : np.ndarray
        Detected peak levels in dB.

    save_path : Path or None
        If provided, save the figure to this path.
    """

    plt.figure(figsize=(10, 4))
    plt.plot(freqs / 1e3, magnitude_db)
    plt.scatter(
        peak_freqs / 1e3,
        peak_levels,
        marker="x",
        label="Detected peaks"
    )
    plt.xlabel("Frequency [kHz]")
    plt.ylabel("Magnitude [dB]")
    plt.title("FFT Spectrum with Detected Peaks")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=200)

    plt.show()


def plot_spectrogram(spectrogram, spectrogram_freqs, spectrogram_times, save_path=None):
    """
    Plot STFT spectrogram.

    Parameters
    ----------
    spectrogram : np.ndarray
        Spectrogram matrix in dB.

    spectrogram_freqs : np.ndarray
        Frequency axis in Hz.

    spectrogram_times : np.ndarray
        Time axis in seconds.

    save_path : Path or None
        If provided, save the figure to this path.
    """

    plt.figure(figsize=(10, 5))
    plt.imshow(
        spectrogram,
        aspect="auto",
        origin="lower",
        extent=[
            spectrogram_times[0],
            spectrogram_times[-1],
            spectrogram_freqs[0] / 1e3,
            spectrogram_freqs[-1] / 1e3,
        ],
    )
    plt.xlabel("Time [s]")
    plt.ylabel("Frequency [kHz]")
    plt.title("Spectrogram of Synthetic IQ Signal")
    plt.colorbar(label="Magnitude [dB]")
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=200)

    plt.show()