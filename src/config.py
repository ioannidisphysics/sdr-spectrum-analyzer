"""
Configuration parameters for the SDR Spectrum Analyzer.
"""

# -----------------------------
# Analyzer settings
# -----------------------------

SAMPLE_RATE_HZ = 1_000_000
DURATION_SECONDS = 0.05
RANDOM_SEED = 42
# -----------------------------
# FFT / spectrum settings
# -----------------------------

PEAK_THRESHOLD_DB = -30
PEAK_MIN_DISTANCE_BINS = 20
NOISE_EXCLUSION_BW_HZ = 5_000

# -----------------------------
# Spectrogram / STFT settings
# -----------------------------

SPECTROGRAM_NFFT = 1024
SPECTROGRAM_HOP = 512