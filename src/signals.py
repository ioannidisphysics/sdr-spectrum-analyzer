"""
Finding and identifying signals in a wideband survey.

Three problems have to be solved before a trace becomes a list of signals.

The noise floor is not a number. Over 24 to 1766 MHz the floor of an RTL-SDR
moves by more than 20 dB: the tuner's gain is not flat, the antenna is
resonant somewhere and deaf elsewhere, and the man made noise that dominates
the low VHF end is absent at 1.5 GHz. A single threshold marks the whole low
end as occupied and misses everything at the top. The floor used here is a
sliding low percentile of the trace, which follows the shape of the receiver
while ignoring the signals sitting on it.

A signal is not a peak. detection.py finds local maxima, which is right for
tones in a laboratory spectrum. A DVB-T multiplex is 7.6 MHz of flat noise
like power with no peak at all, and a peak finder reports it as a few dozen
unrelated bumps. What is wanted is the contiguous run of bins that stand
above the floor, reported once, with its edges.

Not every detection is a signal. An eight bit receiver overloaded by a nearby
transmitter produces harmonics of it, and the 28.8 MHz reference oscillator
leaks into the tuner at its own multiples. Both look exactly like carriers.
They cannot be removed from the trace, but they can be named, and a detection
that carries a plausible explanation is worth more than a silently deleted
one.
"""

import csv

import numpy as np
from scipy.ndimage import percentile_filter

from . import bands

# Reference oscillator of the RTL2832U. Its multiples and the multiples of the
# 28.8 MHz derived clocks appear as narrow birdies with the antenna
# disconnected, which is the test that confirms them.
REFERENCE_CLOCK_HZ = 28.8e6


def window_bins_for(window_hz, bin_hz, n_points):
    """
    Noise floor window in bins, from a width in Hz.

    Specifying the window in Hz rather than in bins keeps it meaningful when
    the survey grid changes. It is forced odd, kept at 21 bins or more so the
    percentile means something, and capped at half the trace so a narrow
    survey does not end up with a window wider than its own data.
    """

    bins = int(round(window_hz / bin_hz)) | 1
    bins = max(21, bins)

    return min(bins, max(21, (n_points // 2) | 1))


def local_noise_floor(trace_db, window_bins=4001, percentile=10.0):
    """
    Noise floor that follows the shape of the receiver.

    A low percentile in a sliding window is used rather than a median. The
    median is pulled upwards wherever signals occupy more than half of the
    window, and across the FM band or inside the UHF television band they do.

    Two things set the window. It has to be wide enough that the widest
    signal to be separated from the floor takes up well under the percentile
    being asked for, and narrow enough to still follow the receiver. At a
    10 kHz grid, 4001 bins is 40 MHz, which leaves a 10 MHz LTE block or a
    7.6 MHz television multiplex occupying a quarter of the window at most,
    so the 10th percentile is taken from bins outside it. The tuner's own gain
    changes over hundreds of megahertz, far more slowly than that.

    The failure case is a band that is continuously occupied over more than
    about nine tenths of the window. There the floor is not visible anywhere
    inside it, and the estimate rises with the signals: wide emissions are
    then reported by their edges instead of as one block. Widening the window
    is the first thing to try.

    Parameters
    ----------
    trace_db : np.ndarray
        Trace in dB. NaN is allowed and is filled by interpolation before
        filtering, so gaps in the grid do not punch holes in the floor.

    window_bins : int
        Sliding window length in bins, forced odd.

    percentile : float
        Percentile taken inside the window.

    Returns
    -------
    np.ndarray
        Noise floor in dB, same length as the input, finite everywhere.
    """

    trace = np.asarray(trace_db, dtype=float)
    window = max(3, int(window_bins) | 1)

    valid = np.isfinite(trace)
    if not valid.any():
        return np.full_like(trace, np.nan)

    filled = trace.copy()
    if not valid.all():
        positions = np.arange(len(trace))
        filled[~valid] = np.interp(positions[~valid], positions[valid], trace[valid])

    return percentile_filter(filled, percentile=percentile, size=window, mode="nearest")


def _parabolic_peak(x, y, index):
    """Refine the position of a maximum by fitting the three points around it."""

    if 0 < index < len(y) - 1:
        y0, y1, y2 = y[index - 1 : index + 2]
        denominator = y0 - 2 * y1 + y2
        if denominator != 0:
            return x[index] + 0.5 * (y0 - y2) / denominator * (x[index + 1] - x[index])

    return x[index]


def _occupied_bandwidth(freqs_hz, power, fraction=0.99):
    """
    Width holding a given fraction of a signal's power, centred on the median.

    This is the occupied bandwidth as ITU-R SM.443 defines it, with the
    remainder split equally between the two sides. It is more stable than the
    width at a fixed level, because it does not depend on where the threshold
    was put.
    """

    total = np.sum(power)
    if total <= 0 or len(freqs_hz) < 2:
        return np.nan

    cumulative = np.cumsum(power) / total
    tail = (1.0 - fraction) / 2.0

    lower = np.interp(tail, cumulative, freqs_hz)
    upper = np.interp(1.0 - tail, cumulative, freqs_hz)

    return float(upper - lower)


def detect_signals(
    centers_hz,
    trace_db,
    floor_db=None,
    snr_db=8.0,
    hysteresis_db=3.0,
    min_width_bins=1,
    merge_gap_bins=3,
    occupancy=None,
    window_bins=4001,
    percentile=10.0,
):
    """
    Contiguous runs of the trace that stand above the local noise floor.

    Parameters
    ----------
    centers_hz : np.ndarray
        Survey grid.

    trace_db : np.ndarray
        Trace to search, normally the max hold.

    floor_db : np.ndarray or None
        Noise floor. Computed with local_noise_floor() if not given.

    snr_db : float
        How far above the floor a bin has to be. 8 dB is a compromise: below
        about 6 dB the noise itself starts producing detections at this grid
        size, above about 12 dB weak but real carriers are lost.

    hysteresis_db : float
        A run is started where the trace crosses floor + snr_db, and then
        extended outwards for as long as it stays above floor + snr_db minus
        this. Without it the flat top of a wide signal, which wanders by a
        couple of dB, is reported as several signals wherever it happens to
        dip through a single threshold.

    min_width_bins : int
        Runs narrower than this are discarded.

    merge_gap_bins : int
        Runs separated by fewer bins than this are joined. A modulated
        carrier dips below the threshold between its components, and without
        merging one transmitter is reported as several signals.

    occupancy : np.ndarray or None
        Per bin occupancy from the survey, summarised per detection.

    window_bins, percentile :
        Passed to local_noise_floor() when floor_db is not given.

    Returns
    -------
    list of dict
        One entry per detection, strongest first. Keys are documented in
        signals_to_csv().
    """

    centers_hz = np.asarray(centers_hz, dtype=float)
    trace_db = np.asarray(trace_db, dtype=float)

    if floor_db is None:
        floor_db = local_noise_floor(trace_db, window_bins=window_bins, percentile=percentile)

    bin_hz = float(np.median(np.diff(centers_hz)))

    finite = np.isfinite(trace_db)
    strong = finite & (trace_db > floor_db + snr_db)
    weak = finite & (trace_db > floor_db + snr_db - hysteresis_db)

    if not strong.any():
        return []

    # ---------------------------------------------
    # Runs of the low threshold that contain at least one high threshold bin,
    # then close the small gaps
    # ---------------------------------------------

    edges = np.diff(weak.astype(np.int8))
    starts = list(np.flatnonzero(edges == 1) + 1)
    stops = list(np.flatnonzero(edges == -1) + 1)

    if weak[0]:
        starts.insert(0, 0)
    if weak[-1]:
        stops.append(len(weak))

    runs = [(start, stop) for start, stop in zip(starts, stops) if strong[start:stop].any()]

    merged = []
    for start, stop in runs:
        if merged and start - merged[-1][1] <= merge_gap_bins:
            merged[-1][1] = stop
        else:
            merged.append([start, stop])

    # ---------------------------------------------
    # Measure each run
    # ---------------------------------------------

    power = np.power(10.0, np.where(np.isfinite(trace_db), trace_db, -300.0) / 10.0)
    floor_power = np.power(10.0, floor_db / 10.0)

    detections = []

    for start, stop in merged:
        if stop - start < min_width_bins:
            continue

        segment = slice(start, stop)
        f_segment = centers_hz[segment]
        db_segment = trace_db[segment]

        local_index = int(np.nanargmax(db_segment))
        peak_db = float(db_segment[local_index])
        f_peak_hz = float(_parabolic_peak(f_segment, db_segment, local_index))

        # Power above the floor only, so a wide signal is not credited with
        # the noise underneath it.
        excess = np.clip(power[segment] - floor_power[segment], 0.0, None)
        channel_power_db = float(10 * np.log10(np.sum(excess) * bin_hz + 1e-30))

        # Power weighted centre. For a modulated signal with no peak of its
        # own, a DVB-T multiplex or a DAB block, the largest bin lands
        # anywhere inside the emission and the centroid is the only stable
        # estimate of where the channel actually is.
        weight = np.sum(excess)
        f_center_hz = float(np.sum(f_segment * excess) / weight) if weight > 0 else f_peak_hz

        entry = {
            "f_peak_hz": f_peak_hz,
            "f_center_hz": f_center_hz,
            "f_lo_hz": float(f_segment[0] - bin_hz / 2),
            "f_hi_hz": float(f_segment[-1] + bin_hz / 2),
            "width_hz": float(len(f_segment) * bin_hz),
            "obw_hz": _occupied_bandwidth(f_segment, excess),
            "peak_db": peak_db,
            "floor_db": float(floor_db[segment][local_index]),
            "snr_db": float(peak_db - floor_db[segment][local_index]),
            "channel_power_db": channel_power_db,
            "bins": int(stop - start),
            "flags": [],
        }

        if occupancy is not None:
            window = np.asarray(occupancy)[segment]
            entry["occupancy_max"] = float(np.nanmax(window))
            entry["occupancy_mean"] = float(np.nanmean(window))
        else:
            entry["occupancy_max"] = np.nan
            entry["occupancy_mean"] = np.nan

        band = bands.lookup(f_peak_hz)
        entry["band"] = band["short"] if band else "unallocated"
        entry["service"] = band["description"] if band else ""
        entry["category"] = band["category"] if band else "unknown"

        channel = bands.fm_channel(f_center_hz)
        entry["fm_channel_hz"] = channel["channel_hz"] if channel else np.nan
        entry["fm_offset_hz"] = channel["offset_hz"] if channel else np.nan

        detections.append(entry)

    detections.sort(key=lambda row: row["peak_db"], reverse=True)

    return detections


def flag_artifacts(detections, bin_hz, clock_hz=REFERENCE_CLOCK_HZ, max_harmonic=6,
                   harmonic_margin_db=15.0):
    """
    Name the detections that are probably the receiver's own doing.

    Two mechanisms produce signals that are not on the air:

    Harmonics. The ADC has eight bits, so a transmitter strong enough to use
    most of that range produces distortion products at integer multiples of
    its frequency. A candidate is flagged when a much stronger detection
    exists at f divided by an integer.

    Reference spurs. The 28.8 MHz oscillator and the clocks derived from it
    leak into the tuner, and appear as narrow birdies at their multiples.

    Both flags are hypotheses, not conclusions. The test that settles a spur
    is to disconnect the antenna: a spur stays, a signal disappears. The test
    for a harmonic is to reduce the gain by 10 dB, which drops a harmonic by
    about 20 or 30 dB and a real signal by 10.

    Parameters
    ----------
    detections : list of dict
        From detect_signals(), modified in place.

    bin_hz : float
        Survey grid cell width, which sets how close a match has to be.

    clock_hz : float
        Reference oscillator frequency.

    max_harmonic : int
        Highest multiple tested.

    harmonic_margin_db : float
        How much stronger the fundamental has to be for the flag to be added.

    Returns
    -------
    list of dict
        The same list, for chaining.
    """

    tolerance_hz = max(3 * bin_hz, 5e3)

    for entry in detections:
        f = entry["f_peak_hz"]

        # Reference oscillator multiples.
        nearest_multiple = round(f / clock_hz)
        if nearest_multiple >= 1 and abs(f - nearest_multiple * clock_hz) <= tolerance_hz:
            entry["flags"].append(f"possible reference spur, {nearest_multiple} x 28.8 MHz")

        # Harmonics of a stronger detection.
        for n in range(2, max_harmonic + 1):
            fundamental = f / n
            for other in detections:
                if other is entry:
                    continue
                if abs(other["f_peak_hz"] - fundamental) > max(tolerance_hz, 1e-4 * fundamental):
                    continue
                if other["peak_db"] >= entry["peak_db"] + harmonic_margin_db:
                    entry["flags"].append(
                        f"possible harmonic {n} of {other['f_peak_hz']/1e6:.3f} MHz"
                    )
                    break
            else:
                continue
            break

    return detections


def summarise_by_band(detections):
    """
    Group detections by band plan entry, for the report.

    Returns
    -------
    list of dict
        band, service, category, count, strongest level, total width and the
        largest occupancy seen, ordered by the lowest frequency in the group.
    """

    groups = {}

    for entry in detections:
        key = (entry["band"], entry["service"], entry["category"])
        groups.setdefault(key, []).append(entry)

    rows = []

    for (band, service, category), members in groups.items():
        rows.append(
            {
                "band": band,
                "service": service,
                "category": category,
                "count": len(members),
                "f_lo_hz": min(m["f_lo_hz"] for m in members),
                "f_hi_hz": max(m["f_hi_hz"] for m in members),
                "peak_db": max(m["peak_db"] for m in members),
                "best_snr_db": max(m["snr_db"] for m in members),
                "total_width_hz": sum(m["width_hz"] for m in members),
                "occupancy_max": np.nanmax([m["occupancy_max"] for m in members]),
            }
        )

    rows.sort(key=lambda row: row["f_lo_hz"])

    return rows


# -----------------------------
# Output
# -----------------------------

CSV_COLUMNS = [
    "f_peak_mhz",
    "f_center_mhz",
    "f_lo_mhz",
    "f_hi_mhz",
    "width_khz",
    "obw_khz",
    "peak_dbfs_per_hz",
    "floor_dbfs_per_hz",
    "snr_db",
    "channel_power_dbfs",
    "occupancy_max",
    "occupancy_mean",
    "band",
    "service",
    "category",
    "fm_channel_mhz",
    "fm_offset_khz",
    "flags",
]


def signals_to_csv(detections, path):
    """
    Write the detection list.

    Levels are in dBFS per Hz. There is no absolute calibration in this
    receiver, so the numbers are comparable within a survey and between
    surveys taken at the same gain, and mean nothing in dBm.
    """

    def mhz(value):
        return "" if value is None or not np.isfinite(value) else f"{value/1e6:.6f}"

    def khz(value):
        return "" if value is None or not np.isfinite(value) else f"{value/1e3:.3f}"

    def db(value):
        return "" if value is None or not np.isfinite(value) else f"{value:.2f}"

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)

        for entry in sorted(detections, key=lambda row: row["f_peak_hz"]):
            writer.writerow(
                [
                    mhz(entry["f_peak_hz"]),
                    mhz(entry["f_center_hz"]),
                    mhz(entry["f_lo_hz"]),
                    mhz(entry["f_hi_hz"]),
                    khz(entry["width_hz"]),
                    khz(entry["obw_hz"]),
                    db(entry["peak_db"]),
                    db(entry["floor_db"]),
                    db(entry["snr_db"]),
                    db(entry["channel_power_db"]),
                    db(entry["occupancy_max"]),
                    db(entry["occupancy_mean"]),
                    entry["band"],
                    entry["service"],
                    entry["category"],
                    mhz(entry["fm_channel_hz"]),
                    khz(entry["fm_offset_hz"]),
                    "; ".join(entry["flags"]),
                ]
            )
