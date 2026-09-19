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

# -----------------------------
# RTL-SDR hardware settings
# -----------------------------

# 2.4 MS/s is the highest rate the RTL-SDR streams without dropping samples on
# most machines. Drop to 2.048e6 if the capture reports gaps.
SDR_SAMPLE_RATE_HZ = 2_400_000

# Fixed manual gain. Must be identical between an antenna sweep and its
# reference sweep, otherwise the receiver's own contribution does not cancel in
# the difference. The nearest value the tuner supports is selected and
# reported. Start here and lower it if clipping is reported.
SDR_GAIN_DB = 30.0

# Samples captured per tuning step. At 2.4 MS/s and nfft 4096 this gives about
# 120 Welch segments, enough for a noise floor steady to a few tenths of a dB.
SDR_SAMPLES_PER_STEP = 262_144

# -----------------------------
# Welch PSD settings
# -----------------------------

SDR_NFFT = 4096
SDR_OVERLAP = 0.5
SDR_WINDOW = "hann"

# Resolution bandwidth follows from these: RBW = ENBW * fs / nfft, and ENBW is
# 1.5 bins for a Hann window. At 2.4 MS/s and nfft 4096 that is about 879 Hz.

# -----------------------------
# Sweep stitching
# -----------------------------

# Fraction of each capture kept, centred, to stay clear of the analogue filter
# skirts at the band edges.
SDR_USABLE_FRACTION = 0.8

# Width discarded at the centre of each capture, where the DC offset and 1/f
# noise of the direct conversion receiver sit.
SDR_DC_EXCLUSION_HZ = 40_000

# Samples read and thrown away after each retune, while the PLL locks.
SDR_SETTLE_SAMPLES = 65_536

# -----------------------------
# Antenna noise comparison
# -----------------------------

# Moving average applied to the difference of two sweeps. Must stay well below
# the width of the antenna's match, which is 8 to 12 percent of centre
# frequency for the dipoles in this measurement.
SDR_SMOOTH_BINS = 201
