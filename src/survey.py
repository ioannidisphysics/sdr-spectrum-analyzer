"""
Repeated wideband survey of everything the receiver can hear.

sweep.py covers a band once and returns a stitched trace. That is the right
tool for a noise floor measurement, where the thing being measured is always
there. It is the wrong tool for finding signals, because most of the radio
spectrum is not transmitting at the instant it is swept: a TETRA channel is
idle between calls, a handset transmits in bursts, a pager fires for a second
a day. A single pass of 24 to 1766 MHz takes minutes, so every channel is
sampled for a few tens of milliseconds out of that, and a single trace says
almost nothing about occupancy.

This module sweeps the band repeatedly and keeps three reductions per
frequency bin, which is what a bench spectrum analyser offers as detectors:

    mean        average power over all passes, the right statistic for noise
    max hold    largest power seen in any pass, which is what catches bursts
    occupancy   fraction of passes in which the bin stood above its own local
                noise floor, which separates a carrier that is always on from
                one that keys up occasionally

Frequency resolution is deliberately reduced before any of this. The Welch
estimate inside a step has a resolution bandwidth of about 880 Hz, and over
1.7 GHz that is two million points: too many to plot, too many to hold for
many passes, and finer than any of the services being looked for. Binning
into 10 kHz cells averages power inside each cell, which lowers the variance
of the noise estimate at the same time.

Averaging is done on power, never on decibels. The mean of two traces in dB
is the geometric mean of the powers, which reads several dB low on a noisy
trace and would make the noise floor depend on how many passes were run.

The tuning grid is staggered between passes. Every step throws away the
40 kHz at its own centre, where the direct conversion receiver's DC offset
sits, so a sweep taken on a fixed grid is blind at the same 900 or so
frequencies on every pass, and a signal unlucky enough to sit on one of them
is never seen at all. Shifting the grid by a fraction of the step width on
each pass moves those holes somewhere else each time, and the accumulated
survey covers everything. It costs one extra tuning step per pass.
"""

import json
import time

import numpy as np

from .psd import welch_psd
from .sdr_source import V4_FREQ_MAX_HZ, V4_FREQ_MIN_HZ, capture_at, check_clipping
from .sweep import plan_sweep

# Columns kept in the stored waterfall image. The full grid is far wider than
# any screen, and the waterfall is read as a picture rather than measured, so
# it is decimated by taking the largest power in each group of columns.
WATERFALL_MAX_COLUMNS = 3000


def make_grid(f_start_hz, f_stop_hz, bin_hz):
    """
    Common frequency grid that every pass is reduced onto.

    Fixing the grid before the first pass is what makes the passes
    comparable. Each pass retunes to the same nominal centres, but the
    stitched points do not land on identical frequencies from one pass to the
    next, and comparing them bin by bin requires a shared axis.

    Parameters
    ----------
    f_start_hz, f_stop_hz : float
        Range to cover.

    bin_hz : float
        Width of each cell.

    Returns
    -------
    edges_hz : np.ndarray
        Bin edges, length n + 1.

    centers_hz : np.ndarray
        Bin centres, length n.
    """

    if f_stop_hz <= f_start_hz:
        raise ValueError("f_stop must be above f_start")

    n_bins = int(np.ceil((f_stop_hz - f_start_hz) / bin_hz))
    edges_hz = f_start_hz + bin_hz * np.arange(n_bins + 1)
    centers_hz = 0.5 * (edges_hz[:-1] + edges_hz[1:])

    return edges_hz, centers_hz


def rebin(freqs_hz, psd_db, edges_hz):
    """
    Reduce a stitched trace onto the survey grid.

    Returns both detectors a spectrum analyser would call average and peak.
    The average detector is the mean power in the cell and is what the noise
    floor should be read from. The peak detector is the largest single point
    in the cell and is what keeps a narrow carrier from being averaged away
    when the cell is much wider than the carrier.

    Parameters
    ----------
    freqs_hz : np.ndarray
        Frequency axis of the stitched trace, ascending.

    psd_db : np.ndarray
        Power spectral density in dBFS per Hz.

    edges_hz : np.ndarray
        Grid edges from make_grid().

    Returns
    -------
    avg_power : np.ndarray
        Mean power per cell, linear, NaN where the cell had no points.

    peak_power : np.ndarray
        Largest power per cell, linear, NaN where the cell had no points.

    counts : np.ndarray
        Points that landed in each cell.
    """

    n_bins = len(edges_hz) - 1
    index = np.searchsorted(edges_hz, freqs_hz, side="right") - 1

    inside = (index >= 0) & (index < n_bins)
    index = index[inside]
    power = np.power(10.0, psd_db[inside] / 10.0)

    counts = np.bincount(index, minlength=n_bins)
    sums = np.bincount(index, weights=power, minlength=n_bins)

    avg_power = np.full(n_bins, np.nan)
    peak_power = np.full(n_bins, np.nan)

    occupied = counts > 0
    avg_power[occupied] = sums[occupied] / counts[occupied]

    # The stitched trace is sorted, so the cell indices are non decreasing and
    # the maximum per cell is a reduceat over the first occurrence of each
    # non empty cell. Empty cells contribute no elements, so consecutive
    # starts delimit exactly one cell's points.
    if index.size:
        starts = np.searchsorted(index, np.arange(n_bins), side="left")
        peak_power[occupied] = np.maximum.reduceat(power, starts[occupied])

    return avg_power, peak_power, counts


class SurveyAccumulator:
    """
    Running mean, max hold and occupancy on a fixed grid.

    Kept as a class so a survey can be interrupted and still written out: the
    accumulator holds everything needed for the report after any number of
    completed passes.
    """

    def __init__(self, centers_hz, occupancy_snr_db=6.0, waterfall_columns=WATERFALL_MAX_COLUMNS):
        self.centers_hz = centers_hz
        self.occupancy_snr_db = occupancy_snr_db

        n = len(centers_hz)
        self.sum_power = np.zeros(n)
        self.sum_counts = np.zeros(n)
        self.max_power = np.full(n, np.nan)
        self.hits = np.zeros(n)
        self.passes = 0
        self.pass_times = []

        self.decimation = max(1, int(np.ceil(n / waterfall_columns)))
        self.waterfall_rows = []

    def add_pass(self, avg_power, peak_power, floor_power, timestamp=None):
        """
        Fold one pass into the accumulators.

        Parameters
        ----------
        avg_power, peak_power : np.ndarray
            Output of rebin() for this pass, linear power, NaN where empty.

        floor_power : np.ndarray
            Local noise floor of this pass, linear power, on the same grid.
            Occupancy is counted against this rather than against a fixed
            level, because the noise floor of the receiver changes by more
            than 20 dB across the tuning range and a fixed threshold would
            simply mark the low frequency end as busy.

        timestamp : float or None
            Unix time at the start of the pass.
        """

        valid = np.isfinite(avg_power)

        self.sum_power[valid] += avg_power[valid]
        self.sum_counts[valid] += 1

        seen = np.isfinite(peak_power)
        self.max_power[seen] = np.fmax(self.max_power[seen], peak_power[seen])

        threshold = floor_power * 10.0 ** (self.occupancy_snr_db / 10.0)
        above = seen & np.isfinite(floor_power) & (peak_power > threshold)
        self.hits[above] += 1

        self.passes += 1
        self.pass_times.append(time.time() if timestamp is None else timestamp)

        self.waterfall_rows.append(self._decimate(peak_power))

    def _decimate(self, power):
        """Largest power in each group of columns, for the waterfall image."""

        d = self.decimation
        if d == 1:
            return power.astype(np.float32)

        padded = np.full(int(np.ceil(len(power) / d)) * d, -np.inf)
        finite = np.isfinite(power)
        padded[: len(power)][finite] = power[finite]

        # -inf rather than NaN so a group with no measured bins does not make
        # numpy warn about an all-NaN slice; it is put back afterwards.
        largest = padded.reshape(-1, d).max(axis=1)
        largest[~np.isfinite(largest)] = np.nan

        return largest.astype(np.float32)

    def result(self):
        """
        Finished traces in dB.

        Returns
        -------
        dict
            mean_db, maxhold_db, occupancy, waterfall_db, waterfall_freqs_hz,
            pass_times.
        """

        with np.errstate(divide="ignore", invalid="ignore"):
            mean_power = np.where(self.sum_counts > 0, self.sum_power / np.maximum(self.sum_counts, 1), np.nan)
            mean_db = 10 * np.log10(mean_power)
            maxhold_db = 10 * np.log10(self.max_power)

            waterfall = np.array(self.waterfall_rows, dtype=np.float32)
            waterfall_db = 10 * np.log10(waterfall) if waterfall.size else waterfall

        occupancy = self.hits / max(self.passes, 1)

        d = self.decimation
        n_cols = waterfall_db.shape[1] if waterfall_db.ndim == 2 else 0
        waterfall_freqs = np.array(
            [np.mean(self.centers_hz[i * d : min((i + 1) * d, len(self.centers_hz))]) for i in range(n_cols)]
        )

        return {
            "mean_db": mean_db,
            "maxhold_db": maxhold_db,
            "occupancy": occupancy,
            "waterfall_db": waterfall_db,
            "waterfall_freqs_hz": waterfall_freqs,
            "pass_times": np.array(self.pass_times),
        }


def estimate_pass_seconds(f_start_hz, f_stop_hz, sample_rate_hz, samples_per_step,
                          settle_samples, usable_fraction, retune_overhead_s=0.04):
    """
    Rough wall clock time for one pass, printed before a survey starts.

    Only the sample counts are known exactly. The retune overhead is the USB
    round trip for setting the centre frequency and is measured at a few tens
    of milliseconds on a typical machine, so the estimate is good to about
    twenty percent, which is enough to decide whether to leave it running.
    """

    centers, _ = plan_sweep(f_start_hz, f_stop_hz, sample_rate_hz, usable_fraction)
    per_step = (samples_per_step + settle_samples) / sample_rate_hz + retune_overhead_s

    return len(centers) * per_step, len(centers)


def run_survey(
    sdr,
    f_start_hz,
    f_stop_hz,
    sample_rate_hz,
    passes,
    bin_hz,
    samples_per_step,
    nfft=4096,
    overlap=0.5,
    window_name="hann",
    usable_fraction=0.8,
    dc_exclusion_hz=40e3,
    settle_samples=65536,
    occupancy_snr_db=6.0,
    floor_window_hz=40e6,
    floor_percentile=10.0,
    progress=True,
):
    """
    Sweep the band `passes` times and accumulate the three detectors.

    Parameters
    ----------
    sdr : object
        Anything with the RtlSdr interface used by capture_at(): a real
        device from sdr_source.open_sdr(), or the simulated source.

    f_start_hz, f_stop_hz : float
        Range to cover.

    sample_rate_hz : float
        Must match what the device was opened with.

    passes : int
        Number of complete sweeps. One pass gives a spectrum; several give
        occupancy, which is the part that distinguishes a live channel from
        a momentary one.

    bin_hz : float
        Survey grid cell width. 10 kHz over the full range is a good default:
        narrow enough to resolve 25 kHz land mobile channels, wide enough to
        keep the trace to a couple of hundred thousand points.

    samples_per_step : int
        Samples per tuning step. Lower than the antenna noise measurement
        uses, because a survey trades noise floor precision for covering the
        band more often.

    nfft, overlap, window_name, usable_fraction, dc_exclusion_hz, settle_samples :
        Passed through to the per step capture and PSD, same meaning as in
        sweep.py.

    occupancy_snr_db : float
        How far above its local noise floor a bin must be, in a given pass,
        to count as occupied in that pass.

    floor_window_hz, floor_percentile :
        Sliding window noise floor used for the occupancy test, see
        signals.local_noise_floor().

    progress : bool
        Print a line per pass.

    Returns
    -------
    centers_hz : np.ndarray
        Survey grid.

    result : dict
        Output of SurveyAccumulator.result().

    meta : dict
        Settings, timing and any clipping warnings.
    """

    from .signals import local_noise_floor, window_bins_for

    edges_hz, centers_hz = make_grid(f_start_hz, f_stop_hz, bin_hz)
    floor_window_bins = window_bins_for(floor_window_hz, bin_hz, len(centers_hz))
    accumulator = SurveyAccumulator(centers_hz, occupancy_snr_db=occupancy_snr_db)

    _, usable_bw_hz = plan_sweep(f_start_hz, f_stop_hz, sample_rate_hz, usable_fraction)

    warnings = []
    rbw_hz = None
    segments = None
    started = time.time()

    for pass_index in range(1, passes + 1):
        pass_started = time.time()
        pass_freqs = []
        pass_psd = []
        clipped_steps = 0

        # Move the tuning grid so this pass is blind in different places than
        # the last one. Starting lower keeps the requested range covered.
        stagger_hz = (pass_index - 1) / passes * usable_bw_hz
        step_centers, _ = plan_sweep(
            f_start_hz - stagger_hz, f_stop_hz, sample_rate_hz, usable_fraction
        )

        # The plan puts half a step of margin beyond each edge of the
        # requested range, which at the ends of the tuner's range is a
        # frequency it cannot reach. Clamping rather than dropping the step
        # keeps the edge covered: the capture simply overlaps its neighbour,
        # and the points that fall outside the requested range are discarded
        # by the rebinning anyway.
        step_centers = np.clip(step_centers, V4_FREQ_MIN_HZ, V4_FREQ_MAX_HZ)
        step_centers = step_centers[np.concatenate(([True], np.diff(step_centers) > 0))]

        for fc in step_centers:
            x = capture_at(sdr, fc, samples_per_step, settle_samples=settle_samples)

            clipped, fraction = check_clipping(x)
            if clipped:
                clipped_steps += 1
                warnings.append(
                    f"pass {pass_index}, {fc/1e6:.1f} MHz: {100*fraction:.2f}% of samples "
                    f"at full scale, reduce gain"
                )

            f_rel, psd_db_step, info = welch_psd(
                x, sample_rate_hz, nfft=nfft, overlap=overlap, window_name=window_name
            )
            rbw_hz = info["rbw_hz"]
            segments = info["segments"]

            keep = (np.abs(f_rel) <= usable_bw_hz / 2) & (np.abs(f_rel) >= dc_exclusion_hz / 2)
            pass_freqs.append(fc + f_rel[keep])
            pass_psd.append(psd_db_step[keep])

        freqs_hz = np.concatenate(pass_freqs)
        psd_db = np.concatenate(pass_psd)

        order = np.argsort(freqs_hz)
        avg_power, peak_power, _ = rebin(freqs_hz[order], psd_db[order], edges_hz)

        with np.errstate(divide="ignore", invalid="ignore"):
            avg_db = 10 * np.log10(avg_power)

        floor_db = local_noise_floor(avg_db, window_bins=floor_window_bins, percentile=floor_percentile)
        floor_power = np.power(10.0, floor_db / 10.0)

        accumulator.add_pass(avg_power, peak_power, floor_power, timestamp=pass_started)

        if progress:
            elapsed = time.time() - pass_started
            occupied_now = 100 * np.mean(accumulator.hits > 0)
            print(
                f"  pass {pass_index:2d}/{passes}  {elapsed:6.1f} s  "
                f"{len(step_centers)} steps  "
                f"{occupied_now:5.1f}% of bins have been active"
                + (f"  CLIPPING on {clipped_steps} steps" if clipped_steps else "")
            )

    result = accumulator.result()

    meta = {
        "f_start_hz": float(f_start_hz),
        "f_stop_hz": float(f_stop_hz),
        "sample_rate_hz": float(sample_rate_hz),
        "bin_hz": float(bin_hz),
        "passes": int(passes),
        "samples_per_step": int(samples_per_step),
        "steps_per_pass": int(len(step_centers)),
        "tuner_range_hz": [float(V4_FREQ_MIN_HZ), float(V4_FREQ_MAX_HZ)],
        "usable_bw_hz": float(usable_bw_hz),
        "staggered": True,
        "dc_exclusion_hz": float(dc_exclusion_hz),
        "nfft": int(nfft),
        "overlap": float(overlap),
        "window": window_name,
        "rbw_hz": float(rbw_hz) if rbw_hz else None,
        "segments_per_step": int(segments) if segments else None,
        "occupancy_snr_db": float(occupancy_snr_db),
        "floor_window_hz": float(floor_window_hz),
        "floor_window_bins": int(floor_window_bins),
        "floor_percentile": float(floor_percentile),
        "waterfall_decimation": int(accumulator.decimation),
        "duration_s": float(time.time() - started),
        "warnings": warnings,
    }

    return centers_hz, result, meta


def save_survey(path, centers_hz, result, meta, settings):
    """Write a survey so the analysis can be repeated without the hardware."""

    np.savez_compressed(
        path,
        centers_hz=centers_hz,
        mean_db=result["mean_db"],
        maxhold_db=result["maxhold_db"],
        occupancy=result["occupancy"],
        waterfall_db=result["waterfall_db"],
        waterfall_freqs_hz=result["waterfall_freqs_hz"],
        pass_times=result["pass_times"],
        meta_json=json.dumps(meta),
        settings_json=json.dumps(settings),
    )


def load_survey(path):
    """
    Read a survey written by save_survey().

    Returns
    -------
    centers_hz, result, meta, settings
    """

    data = np.load(path, allow_pickle=False)

    result = {
        "mean_db": data["mean_db"],
        "maxhold_db": data["maxhold_db"],
        "occupancy": data["occupancy"],
        "waterfall_db": data["waterfall_db"],
        "waterfall_freqs_hz": data["waterfall_freqs_hz"],
        "pass_times": data["pass_times"],
    }

    return (
        data["centers_hz"],
        result,
        json.loads(str(data["meta_json"])),
        json.loads(str(data["settings_json"])),
    )
