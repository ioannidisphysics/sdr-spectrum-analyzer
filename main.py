from src.config import (
    SAMPLE_RATE_HZ,
    DURATION_SECONDS,
    RANDOM_SEED,
    PEAK_THRESHOLD_DB,
    PEAK_MIN_DISTANCE_BINS,
    NOISE_EXCLUSION_BW_HZ,
    SPECTROGRAM_NFFT,
    SPECTROGRAM_HOP,
)

from src.synthetic import generate_synthetic_iq
from src.spectrum import compute_fft_spectrum, compute_spectrogram
from src.detection import (
    detect_spectrum_peaks,
    estimate_noise_floor,
    build_peak_report,
)
from src.io_utils import create_output_dirs, save_peak_report_csv
from src.visualization import (
    plot_time_domain,
    plot_fft_spectrum,
    plot_spectrogram,
)


def print_analyzer_summary(fs, duration, num_samples, metadata):
    """
    Print analyzer settings and synthetic signal information.
    """

    frequency_resolution = fs / num_samples
    bandwidth = fs
    f_min = -fs / 2
    f_max = fs / 2

    print("Sampling frequency:", fs, "Hz")
    print("Duration:", duration, "seconds")
    print("Number of samples:", num_samples)

    print("\nAnalyzer settings:")
    print(f"Bandwidth: {bandwidth / 1e6:.2f} MHz")
    print(f"Frequency span: {f_min / 1e3:.1f} kHz to {f_max / 1e3:.1f} kHz")
    print(f"FFT frequency resolution: {frequency_resolution:.2f} Hz/bin")

    print("\nSynthetic signals:")
    print("Signal 1 frequency:", metadata["f_signal_1"], "Hz")
    print("Signal 2 frequency:", metadata["f_signal_2"], "Hz")
    print("AM carrier frequency:", metadata["f_carrier"], "Hz")
    print("AM modulation frequency:", metadata["f_mod"], "Hz")
    print("AM modulation index:", metadata["mod_index"])
    print("Chirp start frequency:", metadata["f_chirp_start"], "Hz")
    print("Chirp end frequency:", metadata["f_chirp_end"], "Hz")
    print("Noise power:", metadata["noise_power"])
    print("Random seed:", metadata["random_seed"])

def print_peak_report(report, noise_floor_db):
    """
    Print detected peaks and SNR values.
    """

    print("\nEstimated noise floor:")
    print(f"{noise_floor_db:.2f} dB")

    print("\nDetected peaks:")

    for row in report:
        print(
            f"{row['frequency_khz']:8.2f} kHz | "
            f"Level: {row['level_db']:6.2f} dB | "
            f"SNR: {row['snr_db']:6.2f} dB"
        )


def main():
    print("SDR Spectrum Analyzer started")

    # -----------------------------
    # Output directories
    # -----------------------------

    output_dir, figures_dir = create_output_dirs()

    # -----------------------------
    # Analyzer configuration
    # -----------------------------

    fs = SAMPLE_RATE_HZ
    duration = DURATION_SECONDS
    # -----------------------------
    # Generate synthetic IQ signal
    # -----------------------------

    t, x, metadata = generate_synthetic_iq(
        fs=fs,
        duration=duration,
        random_seed=RANDOM_SEED,
    )
    print_analyzer_summary(
        fs=fs,
        duration=duration,
        num_samples=len(x),
        metadata=metadata,
    )

    # -----------------------------
    # FFT spectrum
    # -----------------------------

    freqs, magnitude_db = compute_fft_spectrum(x, fs)

    # -----------------------------
    # Peak detection and SNR
    # -----------------------------

    peak_freqs, peak_levels = detect_spectrum_peaks(
        freqs=freqs,
        magnitude_db=magnitude_db,
        threshold_db=PEAK_THRESHOLD_DB,
        distance=PEAK_MIN_DISTANCE_BINS,
    )

    noise_floor_db = estimate_noise_floor(
        freqs=freqs,
        magnitude_db=magnitude_db,
        peak_freqs=peak_freqs,
        exclusion_bw=NOISE_EXCLUSION_BW_HZ,
    )

    report = build_peak_report(
        peak_freqs=peak_freqs,
        peak_levels=peak_levels,
        noise_floor_db=noise_floor_db,
    )

    print_peak_report(report, noise_floor_db)

    # -----------------------------
    # Save CSV report
    # -----------------------------

    csv_path = output_dir / "detected_peaks.csv"
    save_peak_report_csv(report, csv_path)

    print(f"\nSaved peak report to: {csv_path}")

    # -----------------------------
    # Spectrogram
    # -----------------------------

    spectrogram, spectrogram_freqs, spectrogram_times = compute_spectrogram(
        x=x,
        fs=fs,
        nfft=SPECTROGRAM_NFFT,
        hop=SPECTROGRAM_HOP,
    )

    # -----------------------------
    # Plots
    # -----------------------------

    plot_time_domain(
        t=t,
        x=x,
        save_path=figures_dir / "time_domain_iq.png",
    )

    plot_fft_spectrum(
        freqs=freqs,
        magnitude_db=magnitude_db,
        peak_freqs=peak_freqs,
        peak_levels=peak_levels,
        save_path=figures_dir / "fft_spectrum.png",
    )

    plot_spectrogram(
        spectrogram=spectrogram,
        spectrogram_freqs=spectrogram_freqs,
        spectrogram_times=spectrogram_times,
        save_path=figures_dir / "spectrogram.png",
    )

    print(f"Saved figures to: {figures_dir}")


if __name__ == "__main__":
    main()