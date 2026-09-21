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

# -----------------------------
# Wideband survey
# -----------------------------

# Full tuning range of the RTL-SDR Blog V4. The tuner covers 24 to 1766 MHz;
# anything above that, including the 2.45 GHz the patch antenna is designed
# for, needs a downconverter in front of the dongle.
SURVEY_START_HZ = 24e6
SURVEY_STOP_HZ = 1766e6

# Survey grid. 10 kHz resolves the 12.5 and 25 kHz land mobile channels and
# keeps the full range to about 175 000 points, which plots and stores
# comfortably. The Welch estimate inside each step is still about 880 Hz; the
# grid averages it down in power, which also lowers the variance.
SURVEY_BIN_HZ = 10e3

# Passes over the band. One gives a spectrum, several give occupancy. Three is
# enough to separate an always on carrier from a one off burst; ten or more
# makes the waterfall worth looking at.
SURVEY_PASSES = 3

# Samples per tuning step during a survey. Lower than the antenna noise
# measurement uses: a survey trades a few tenths of a dB on the noise floor
# for covering the band more often.
SURVEY_SAMPLES_PER_STEP = 131_072

# A bin counts as occupied in a pass when it stands this far above the local
# noise floor of that pass.
SURVEY_OCCUPANCY_SNR_DB = 6.0

# Threshold for the detection list, applied to the max hold trace. Below about
# 6 dB the noise itself starts producing detections at this grid size.
SURVEY_DETECT_SNR_DB = 8.0

# Sliding window noise floor, given in Hz so it keeps its meaning when the
# grid changes. 40 MHz leaves the widest emission worth separating from the
# floor, a 10 MHz LTE block or a 7.6 MHz television multiplex, taking up a
# quarter of the window at most, so the 10th percentile comes from bins
# outside it. The tuner's own response changes over hundreds of megahertz.
SURVEY_FLOOR_WINDOW_HZ = 40e6
SURVEY_FLOOR_PERCENTILE = 10.0

# Runs of bins separated by less than this are reported as one signal. Five
# bins is 50 kHz at the default grid, just wider than the 40 kHz thrown away
# at the centre of each tuning step, so a single pass survey still reports a
# wide signal once even where a DC hole falls inside it. With more than one
# pass the staggered grid fills those holes in and this only has to bridge the
# dips in a modulated signal.
SURVEY_MERGE_GAP_BINS = 5

# Report layout.
SURVEY_OVERVIEW_ROWS = 4
SURVEY_ZOOM_BANDS = 6
