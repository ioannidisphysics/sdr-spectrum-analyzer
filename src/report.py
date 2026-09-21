"""
Markdown report for a wideband survey.

Written as a file rather than printed, because the useful form of this output
is something that can be committed next to the figures and read later. The
terminal gets a short summary; the file gets everything, including the
settings, so a result can be traced back to the run that produced it.
"""

import datetime as _dt

import numpy as np

from . import bands
from .signals import summarise_by_band


def _mhz(value):
    return f"{value/1e6:.4f}"


def _percent(value):
    return "" if not np.isfinite(value) else f"{100*value:.0f}%"


def build_markdown(
    centers_hz,
    result,
    meta,
    settings,
    detections,
    floor_db,
    detect_snr_db,
    figures=(),
    label="survey",
    top_n=30,
):
    """
    Assemble the report text.

    Parameters
    ----------
    centers_hz : np.ndarray
        Survey grid.

    result : dict
        From survey.load_survey().

    meta, settings : dict
        Run settings.

    detections : list of dict
        From signals.detect_signals(), already passed through
        signals.flag_artifacts().

    floor_db : np.ndarray
        Noise floor used for the detection.

    detect_snr_db : float
        Threshold above the floor.

    figures : iterable of str
        Paths, relative to the repository root, to link.

    label : str
        Name of the survey.

    top_n : int
        Rows in the strongest signals table.

    Returns
    -------
    str
        Markdown text.
    """

    simulated = settings.get("source") == "simulated"
    maxhold = result["maxhold_db"]
    occupancy = result["occupancy"]

    finite = np.isfinite(maxhold)
    above = finite & (maxhold > floor_db + detect_snr_db)

    lines = []
    add = lines.append

    add(f"# Wideband survey: {label}")
    add("")

    if simulated:
        add("> **Simulated source.** These numbers come from `src/simulated.py`, not from")
        add("> the radio. They are here to exercise the analysis, and nothing in this")
        add("> report describes the actual radio environment.")
        add("")

    add(f"Generated {_dt.datetime.now().strftime('%Y-%m-%d %H:%M')}.")
    add("")

    # -----------------------------
    # Settings
    # -----------------------------

    add("## Run")
    add("")
    add("| Setting | Value |")
    add("| --- | --- |")
    add(f"| Range | {meta['f_start_hz']/1e6:.1f} to {meta['f_stop_hz']/1e6:.1f} MHz |")
    add(f"| Passes | {meta['passes']} |")
    add(f"| Tuning steps per pass | {meta['steps_per_pass']} |")
    add(f"| Sample rate | {meta['sample_rate_hz']/1e6:.3f} MS/s |")
    add(f"| Samples per step | {meta['samples_per_step']} |")
    add(f"| Welch FFT size | {meta['nfft']}, {meta['window']} window, {meta['overlap']:.0%} overlap |")
    rbw = meta.get("rbw_hz")
    add(f"| Resolution bandwidth | {rbw:.0f} Hz |" if rbw else "| Resolution bandwidth | unknown |")
    add(f"| Survey grid | {meta['bin_hz']/1e3:.1f} kHz, {len(centers_hz)} bins |")
    add(f"| Tuner gain | {settings.get('gain_db_applied', float('nan')):.1f} dB, AGC off |")
    add(f"| Occupancy threshold | floor + {meta['occupancy_snr_db']:.0f} dB |")
    add(f"| Detection threshold | floor + {detect_snr_db:.0f} dB on the max hold |")
    add(f"| Noise floor estimator | {meta['floor_percentile']:.0f}th percentile over "
        f"{meta['floor_window_bins'] * meta['bin_hz']/1e6:.0f} MHz |")
    add(f"| Duration | {meta['duration_s']/60:.1f} min |")
    add("")

    if meta.get("warnings"):
        add("### Warnings")
        add("")
        for warning in meta["warnings"][:20]:
            add(f"- {warning}")
        if len(meta["warnings"]) > 20:
            add(f"- ... and {len(meta['warnings']) - 20} more")
        add("")

    # -----------------------------
    # Headline numbers
    # -----------------------------

    add("## Result")
    add("")
    add(f"- Signals detected: **{len(detections)}**")
    add(f"- Spectrum above the detection threshold: **{100*np.mean(above[finite]):.2f}%** of the range")
    add(f"- Bins active in at least one pass: **{100*np.mean(occupancy > 0):.2f}%**")
    add(f"- Bins active in every pass: **{100*np.mean(occupancy >= 1.0):.2f}%**")

    if detections:
        strongest = max(detections, key=lambda d: d["peak_db"])
        add(
            f"- Strongest signal: **{_mhz(strongest['f_peak_hz'])} MHz**, "
            f"{strongest['peak_db']:.1f} dBFS/Hz, {strongest['snr_db']:.1f} dB above the floor "
            f"({strongest['band']})"
        )

    flagged = [d for d in detections if d["flags"]]
    add(f"- Detections flagged as probable receiver artefacts: **{len(flagged)}**")
    add("")

    # -----------------------------
    # By band
    # -----------------------------

    add("## By band")
    add("")
    add("| Band | Service | Signals | Strongest (dBFS/Hz) | Best SNR (dB) | Max occupancy |")
    add("| --- | --- | --- | --- | --- | --- |")

    for row in summarise_by_band(detections):
        add(
            f"| {row['band']} | {row['service']} | {row['count']} | "
            f"{row['peak_db']:.1f} | {row['best_snr_db']:.1f} | {_percent(row['occupancy_max'])} |"
        )

    add("")

    # -----------------------------
    # Strongest signals
    # -----------------------------

    add(f"## Strongest {min(top_n, len(detections))} signals")
    add("")
    add("| Frequency (MHz) | Width (kHz) | OBW 99% (kHz) | Peak (dBFS/Hz) | SNR (dB) | Occupancy | Band | Notes |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- |")

    for entry in sorted(detections, key=lambda d: d["peak_db"], reverse=True)[:top_n]:
        obw = "" if not np.isfinite(entry["obw_hz"]) else f"{entry['obw_hz']/1e3:.1f}"
        notes = "; ".join(entry["flags"])

        if np.isfinite(entry["fm_offset_hz"]):
            offset = f"{entry['fm_offset_hz']/1e3:+.1f} kHz from the {_mhz(entry['fm_channel_hz'])} channel"
            notes = f"{notes}; {offset}" if notes else offset

        add(
            f"| {_mhz(entry['f_peak_hz'])} | {entry['width_hz']/1e3:.1f} | {obw} | "
            f"{entry['peak_db']:.1f} | {entry['snr_db']:.1f} | "
            f"{_percent(entry['occupancy_max'])} | {entry['band']} | {notes} |"
        )

    add("")

    # -----------------------------
    # Artefacts
    # -----------------------------

    if flagged:
        add("## Probable receiver artefacts")
        add("")
        add("These are hypotheses, not conclusions. A reference spur stays when the")
        add("antenna is disconnected and a real signal does not. A harmonic falls by")
        add("about 20 to 30 dB when the gain is reduced by 10 dB, a real signal by 10.")
        add("")
        add("| Frequency (MHz) | Peak (dBFS/Hz) | Explanation |")
        add("| --- | --- | --- |")

        for entry in sorted(flagged, key=lambda d: d["f_peak_hz"]):
            add(f"| {_mhz(entry['f_peak_hz'])} | {entry['peak_db']:.1f} | {'; '.join(entry['flags'])} |")

        add("")

    # -----------------------------
    # Frequency axis check
    # -----------------------------

    fm_offsets = [d["fm_offset_hz"] for d in detections if np.isfinite(d.get("fm_offset_hz", np.nan))]

    if len(fm_offsets) >= 3:
        offsets = np.array(fm_offsets)
        median_hz = float(np.median(offsets))
        reference_hz = float(np.median([d["f_peak_hz"] for d in detections
                                        if np.isfinite(d.get("fm_offset_hz", np.nan))]))
        ppm = 1e6 * median_hz / reference_hz

        add("## Frequency axis check")
        add("")
        add("FM broadcast carriers in Region 1 sit on odd multiples of 100 kHz, which")
        add("makes the band a free frequency reference. The offset between the detected")
        add("carriers and their channels is the crystal error of the receiver plus what")
        add("the sweep stitching contributes.")
        add("")
        add(f"- Carriers matched: **{len(offsets)}**")
        add(f"- Median offset: **{median_hz/1e3:+.2f} kHz**, which at {reference_hz/1e6:.0f} MHz "
            f"is **{ppm:+.1f} ppm**")
        add(f"- Spread: {np.std(offsets, ddof=1)/1e3:.2f} kHz standard deviation")
        add("")
        add(f"Passing `--ppm {ppm:+.0f}` to the capture removes most of this.")
        add("")

    # -----------------------------
    # Figures
    # -----------------------------

    if figures:
        add("## Figures")
        add("")
        for path in figures:
            add(f"![{path}]({path})")
            add("")

    # -----------------------------
    # What the numbers are not
    # -----------------------------

    add("## Limits")
    add("")
    add("- Levels are dBFS per Hz. The receiver has no absolute calibration, so they")
    add("  compare within a survey and against another survey at the same gain, and")
    add("  are not dBm. Converting them would need a source of known power.")
    add("- The trace is the product of the antenna, the cable and the receiver as much")
    add("  as of what is transmitting. A band that looks empty may be a band the")
    add("  antenna cannot hear.")
    add("- Occupancy is sampled, not monitored. Each bin is observed for a few tens of")
    add("  milliseconds per pass, so a channel that transmits rarely can be missed")
    add(f"  entirely: with {meta['passes']} passes the survey saw each bin {meta['passes']} times.")
    add("- An eight bit receiver overloaded by a strong transmitter invents harmonics.")
    add("  The flags above are the first check, reducing the gain is the second.")
    add(f"- The tuner reaches {bands.RX_MAX_HZ/1e6:.0f} MHz. Wi-Fi, Bluetooth and the")
    add("  2.45 GHz patch antenna are above that and need a downconverter.")
    add("")

    return "\n".join(lines)
