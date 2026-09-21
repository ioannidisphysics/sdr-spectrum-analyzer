"""
Real hardware entry point for the SDR Spectrum Analyzer.

main.py analyses a synthetic signal. This script drives an RTL-SDR and has
two jobs.

Antenna measurement. Find where a dipole is resonant by reading the rise in
the receiver's noise floor when the antenna is connected instead of a 50 ohm
load, with no absolute calibration needed anywhere:

    python main_sdr.py capture --label A_antenna --start 100e6 --stop 190e6
    python main_sdr.py capture --label A_load    --start 100e6 --stop 190e6
    python main_sdr.py compare --antenna A_antenna --reference A_load \
                               --vna ../dipole-sim-vs-measurement/data/A_meas_1.s1p

Wideband survey. Sweep everything the tuner can reach, several times over,
and report what is on the air, how wide it is, how often it was there and
which service it belongs to:

    python main_sdr.py survey --label home --passes 5
    python main_sdr.py report --label home

Capture and analysis are separate commands throughout, so a measurement can
be re-analysed and re-plotted later without the hardware present. Run

    python main_sdr.py selftest

to check the analysis chain against a source of known power spectral density,
which needs no hardware at all.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from src.antenna_noise import (
    antenna_noise_delta,
    find_peak,
    load_s1p_min,
    plot_antenna_noise,
    rise_summary,
)
from src.config import (
    SDR_SAMPLE_RATE_HZ,
    SDR_GAIN_DB,
    SDR_SAMPLES_PER_STEP,
    SDR_NFFT,
    SDR_OVERLAP,
    SDR_WINDOW,
    SDR_USABLE_FRACTION,
    SDR_DC_EXCLUSION_HZ,
    SDR_SETTLE_SAMPLES,
    SDR_SMOOTH_BINS,
    SURVEY_START_HZ,
    SURVEY_STOP_HZ,
    SURVEY_BIN_HZ,
    SURVEY_PASSES,
    SURVEY_SAMPLES_PER_STEP,
    SURVEY_OCCUPANCY_SNR_DB,
    SURVEY_DETECT_SNR_DB,
    SURVEY_FLOOR_WINDOW_HZ,
    SURVEY_FLOOR_PERCENTILE,
    SURVEY_MERGE_GAP_BINS,
    SURVEY_OVERVIEW_ROWS,
    SURVEY_ZOOM_BANDS,
)
from src.raster_check import check_uhf_raster, find_blocks, format_report
from src.sdr_source import open_sdr
from src.sweep import load_sweep, run_sweep, save_sweep
from src.survey import estimate_pass_seconds, load_survey, run_survey, save_survey

SWEEP_DIR = Path("outputs") / "sweeps"
SURVEY_DIR = Path("outputs") / "surveys"
FIGURE_DIR = Path("outputs") / "figures"


# ---------------------------------------------------------------
# Antenna measurement
# ---------------------------------------------------------------

def cmd_capture(args):
    """Sweep the band once and save it."""

    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    out_path = SWEEP_DIR / f"{args.label}.npz"

    if out_path.exists() and not args.overwrite:
        sys.exit(f"{out_path} already exists, pass --overwrite to replace it")

    print(f"Opening RTL-SDR at {args.gain:.1f} dB gain, {args.sample_rate/1e6:.3f} MS/s")
    sdr, settings = open_sdr(args.sample_rate, args.gain, args.ppm)

    if settings["gain_db_applied"] != settings["gain_db_requested"]:
        print(
            f"  gain snapped to the nearest supported value: "
            f"{settings['gain_db_applied']:.1f} dB"
        )

    try:
        print(f"Sweeping {args.start/1e6:.1f} to {args.stop/1e6:.1f} MHz")
        freqs_hz, psd_db, meta = run_sweep(
            sdr,
            f_start_hz=args.start,
            f_stop_hz=args.stop,
            sample_rate_hz=args.sample_rate,
            samples_per_step=args.samples,
            nfft=SDR_NFFT,
            overlap=SDR_OVERLAP,
            window_name=SDR_WINDOW,
            usable_fraction=SDR_USABLE_FRACTION,
            dc_exclusion_hz=SDR_DC_EXCLUSION_HZ,
            settle_samples=SDR_SETTLE_SAMPLES,
        )
    finally:
        sdr.close()

    print(f"\nResolution bandwidth: {meta['rbw_hz']/1e3:.2f} kHz")
    print(f"Welch segments per step: {meta['segments_per_step']}")
    print(f"Points in the stitched trace: {len(freqs_hz)}")

    for warning in meta["warnings"]:
        print(f"WARNING  {warning}")

    save_sweep(out_path, freqs_hz, psd_db, meta, settings)
    print(f"Saved {out_path}")


def cmd_compare(args):
    """Difference two sweeps and locate the antenna's resonance."""

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    f_ant, psd_ant, meta_ant, set_ant = load_sweep(SWEEP_DIR / f"{args.antenna}.npz")
    f_ref, psd_ref, meta_ref, set_ref = load_sweep(SWEEP_DIR / f"{args.reference}.npz")

    _check_pair(meta_ant, set_ant, meta_ref, set_ref)

    freqs_hz, delta_db, delta_smooth_db = antenna_noise_delta(
        f_ant, psd_ant, f_ref, psd_ref, smooth_bins=args.smooth
    )

    f_peak_hz, level_db = find_peak(
        freqs_hz, delta_smooth_db, args.search_lo, args.search_hi
    )
    rise = rise_summary(freqs_hz, delta_smooth_db, args.search_lo, args.search_hi,
                        args.min_rise)

    print(f"Gain (both sweeps):       {set_ant['gain_db_applied']:.1f} dB")
    print(f"Resolution bandwidth:     {meta_ant['rbw_hz']/1e3:.2f} kHz")
    print(f"Smoothing:                {args.smooth} bins")
    print(f"Broadband noise rise:     {rise['median_db']:+.2f} dB median, "
          f"{rise['p90_db']:+.2f} dB at the 90th percentile")
    print(f"Noise rise peak:          {f_peak_hz/1e6:.2f} MHz at {level_db:.2f} dB")

    if not rise["usable"]:
        print(
            f"\nRECEIVER NOISE LIMITED. The antenna raises the floor by only "
            f"{rise['median_db']:+.2f} dB, against the {rise['usable_db']:.1f} dB "
            f"this method needs.\n"
            f"The receiver's own noise is louder than what the antenna delivers, so the "
            f"difference is flat except at transmitters and the peak above is a\n"
            f"broadcast carrier, not a resonance. Raise the tuner gain and repeat, "
            f"stopping at the highest gain that does not report clipping.\n"
        )

    f_vna_hz = None
    if args.vna:
        result = load_s1p_min(args.vna)
        if result is None:
            print("scikit-rf not installed, skipping the VNA comparison")
        else:
            f_vna_hz, s11_db = result
            error_pct = 100 * (f_peak_hz - f_vna_hz) / f_vna_hz
            print(f"VNA |S11| minimum:        {f_vna_hz/1e6:.2f} MHz at {s11_db:.2f} dB")
            if rise["usable"]:
                print(f"Noise peak vs VNA:        {error_pct:+.2f} %")

    save_path = FIGURE_DIR / f"{args.antenna}_vs_{args.reference}.png"
    psd_ref_on_grid = np.interp(freqs_hz, f_ref, psd_ref)

    plot_antenna_noise(
        freqs_hz=freqs_hz,
        psd_ant_db=psd_ant,
        psd_ref_db=psd_ref_on_grid,
        delta_db=delta_db,
        delta_smooth_db=delta_smooth_db,
        save_path=save_path,
        f_peak_hz=f_peak_hz,
        f_vna_hz=f_vna_hz,
        title=f"{args.antenna} minus {args.reference}",
    )
    print(f"Saved {save_path}")


def cmd_raster(args):
    """Check the frequency axis against the UHF television channel plan."""

    f_ant, psd_ant, meta_ant, set_ant = load_sweep(SWEEP_DIR / f"{args.antenna}.npz")
    f_ref, psd_ref, meta_ref, set_ref = load_sweep(SWEEP_DIR / f"{args.reference}.npz")
    _check_pair(meta_ant, set_ant, meta_ref, set_ref)

    freqs_hz, _, delta_smooth_db = antenna_noise_delta(
        f_ant, psd_ant, f_ref, psd_ref, smooth_bins=args.smooth
    )

    blocks = find_blocks(freqs_hz, delta_smooth_db, args.threshold)
    rows, summary = check_uhf_raster(blocks)
    print(format_report(rows, summary))


def _check_pair(meta_a, set_a, meta_b, set_b):
    """Refuse to subtract two sweeps that were not taken the same way."""

    problems = []

    if set_a["gain_db_applied"] != set_b["gain_db_applied"]:
        problems.append(
            f"gain differs: {set_a['gain_db_applied']} vs {set_b['gain_db_applied']} dB"
        )
    if meta_a["sample_rate_hz"] != meta_b["sample_rate_hz"]:
        problems.append("sample rate differs")
    if meta_a["rbw_hz"] != meta_b["rbw_hz"]:
        problems.append("resolution bandwidth differs")
    if abs(meta_a["f_start_hz"] - meta_b["f_start_hz"]) > 1e3:
        problems.append("start frequency differs")
    if abs(meta_a["f_stop_hz"] - meta_b["f_stop_hz"]) > 1e3:
        problems.append("stop frequency differs")

    if problems:
        sys.exit(
            "These two sweeps cannot be compared:\n  "
            + "\n  ".join(problems)
            + "\nThe difference only cancels the receiver if both runs used identical settings."
        )


# ---------------------------------------------------------------
# Wideband survey
# ---------------------------------------------------------------

def cmd_survey(args):
    """Sweep the whole tuning range repeatedly and store the result."""

    SURVEY_DIR.mkdir(parents=True, exist_ok=True)
    out_path = SURVEY_DIR / f"{args.label}.npz"

    if out_path.exists() and not args.overwrite:
        sys.exit(f"{out_path} already exists, pass --overwrite to replace it")

    if args.simulate:
        from src.simulated import open_simulated

        print("Simulated source. Nothing here is a measurement.")
        sdr, settings = open_simulated(args.sample_rate, args.gain, args.ppm, seed=args.seed)
    else:
        print(f"Opening RTL-SDR at {args.gain:.1f} dB gain, {args.sample_rate/1e6:.3f} MS/s")
        sdr, settings = open_sdr(args.sample_rate, args.gain, args.ppm)

        if settings["gain_db_applied"] != settings["gain_db_requested"]:
            print(f"  gain snapped to {settings['gain_db_applied']:.1f} dB")

    seconds, steps = estimate_pass_seconds(
        args.start, args.stop, args.sample_rate, args.samples,
        SDR_SETTLE_SAMPLES, SDR_USABLE_FRACTION,
    )

    print(
        f"{args.start/1e6:.1f} to {args.stop/1e6:.1f} MHz, {steps} steps per pass, "
        f"{args.passes} passes"
    )
    print(
        f"Estimated {seconds/60:.1f} min per pass, {args.passes*seconds/60:.1f} min in total"
    )

    try:
        centers_hz, result, meta = run_survey(
            sdr,
            f_start_hz=args.start,
            f_stop_hz=args.stop,
            sample_rate_hz=args.sample_rate,
            passes=args.passes,
            bin_hz=args.bin,
            samples_per_step=args.samples,
            nfft=SDR_NFFT,
            overlap=SDR_OVERLAP,
            window_name=SDR_WINDOW,
            usable_fraction=SDR_USABLE_FRACTION,
            dc_exclusion_hz=SDR_DC_EXCLUSION_HZ,
            settle_samples=SDR_SETTLE_SAMPLES,
            occupancy_snr_db=args.occupancy_snr,
            floor_window_hz=SURVEY_FLOOR_WINDOW_HZ,
            floor_percentile=SURVEY_FLOOR_PERCENTILE,
        )
    finally:
        sdr.close()

    print(f"\nGrid: {len(centers_hz)} bins of {args.bin/1e3:.1f} kHz")
    print(f"Resolution bandwidth inside each step: {meta['rbw_hz']:.0f} Hz")
    print(f"Bins active in at least one pass: {100*np.mean(result['occupancy'] > 0):.2f}%")
    print(f"Took {meta['duration_s']/60:.1f} min")

    if meta["warnings"]:
        print(f"\n{len(meta['warnings'])} clipping warnings, first few:")
        for warning in meta["warnings"][:5]:
            print(f"  {warning}")

    save_survey(out_path, centers_hz, result, meta, settings)
    print(f"Saved {out_path}")
    print(f"Now run:  python main_sdr.py report --label {args.label}")


def cmd_report(args):
    """Detect, identify and plot the signals in a stored survey."""

    from src.report import build_markdown
    from src.signals import (
        detect_signals,
        flag_artifacts,
        local_noise_floor,
        signals_to_csv,
        window_bins_for,
    )
    from src.waterfall import busiest_bands, plot_band, plot_overview, plot_waterfall

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    survey_path = SURVEY_DIR / f"{args.label}.npz"
    if not survey_path.exists():
        sys.exit(f"{survey_path} not found. Run the survey command first.")

    centers_hz, result, meta, settings = load_survey(survey_path)
    simulated = settings.get("source") == "simulated"

    maxhold_db = result["maxhold_db"]
    mean_db = result["mean_db"]

    floor_window_bins = window_bins_for(
        args.floor_window * 1e6, meta["bin_hz"], len(centers_hz)
    )
    floor_db = local_noise_floor(
        maxhold_db,
        window_bins=floor_window_bins,
        percentile=meta.get("floor_percentile", SURVEY_FLOOR_PERCENTILE),
    )
    print(
        f"Noise floor: {meta.get('floor_percentile', SURVEY_FLOOR_PERCENTILE):.0f}th percentile "
        f"over {floor_window_bins * meta['bin_hz']/1e6:.1f} MHz"
    )

    detections = detect_signals(
        centers_hz,
        maxhold_db,
        floor_db=floor_db,
        snr_db=args.snr,
        merge_gap_bins=args.merge_gap,
        occupancy=result["occupancy"],
    )
    flag_artifacts(detections, bin_hz=meta["bin_hz"])

    print(f"{len(detections)} signals above floor + {args.snr:.0f} dB")
    print(f"{sum(1 for d in detections if d['flags'])} flagged as probable receiver artefacts")
    print()
    print(f"{'Frequency':>13}  {'Width':>9}  {'SNR':>6}  {'Occ':>5}  Band")

    for entry in sorted(detections, key=lambda d: d["peak_db"], reverse=True)[: args.list]:
        occupancy = entry["occupancy_max"]
        occupancy_text = "-" if not np.isfinite(occupancy) else f"{100*occupancy:.0f}%"
        print(
            f"{entry['f_peak_hz']/1e6:10.4f} MHz  {entry['width_hz']/1e3:6.1f} kHz  "
            f"{entry['snr_db']:5.1f}  {occupancy_text:>5}  {entry['band']}"
        )

    # -----------------------------
    # Files
    # -----------------------------

    csv_path = SURVEY_DIR / f"{args.label}_signals.csv"
    signals_to_csv(detections, csv_path)

    figures = []

    overview_path = FIGURE_DIR / f"{args.label}_overview.png"
    plot_overview(
        centers_hz, mean_db, maxhold_db, floor_db, detections, overview_path,
        rows=args.rows,
        title=f"Survey {args.label}: {meta['f_start_hz']/1e6:.0f} to {meta['f_stop_hz']/1e6:.0f} MHz",
        subtitle=(
            f"{meta['passes']} passes, {meta['bin_hz']/1e3:.0f} kHz grid, "
            f"gain {settings.get('gain_db_applied', float('nan')):.1f} dB, "
            f"detection threshold floor + {args.snr:.0f} dB"
        ),
        simulated=simulated,
    )
    figures.append(overview_path.as_posix())

    waterfall_path = FIGURE_DIR / f"{args.label}_waterfall.png"
    if plot_waterfall(
        result["waterfall_freqs_hz"], result["waterfall_db"], waterfall_path,
        pass_times=result["pass_times"],
        title=f"Survey {args.label}, {meta['passes']} passes",
        simulated=simulated,
    ):
        figures.append(waterfall_path.as_posix())

    for band in busiest_bands(detections, limit=args.zooms):
        name = band["short"].lower().replace(" ", "_").replace("/", "-")
        path = FIGURE_DIR / f"{args.label}_band_{name}.png"

        margin = 0.02 * (band["f_hi_hz"] - band["f_lo_hz"])
        drawn = plot_band(
            centers_hz, mean_db, maxhold_db, floor_db, result["occupancy"], detections,
            band["f_lo_hz"] - margin, band["f_hi_hz"] + margin, path,
            title=f"{band['short']}: {band['f_lo_hz']/1e6:.1f} to {band['f_hi_hz']/1e6:.1f} MHz, "
                  f"{band['count']} signals",
            snr_db=args.snr,
            simulated=simulated,
        )
        if drawn:
            figures.append(path.as_posix())

    report_path = SURVEY_DIR / f"{args.label}_report.md"
    report_path.write_text(
        build_markdown(
            centers_hz, result, meta, settings, detections, floor_db, args.snr,
            figures=figures, label=args.label,
        ),
        encoding="utf-8",
    )

    print()
    print(f"Saved {csv_path}")
    print(f"Saved {report_path}")
    for path in figures:
        print(f"Saved {path}")


def cmd_selftest(args):
    """Check the analysis chain against a source of known density."""

    from src.simulated import selftest

    print("Comparing the analyser against a simulated source whose power spectral")
    print("density is known by construction. No hardware involved.")
    print()

    worst = 0.0

    for row in selftest(sample_rate_hz=args.sample_rate, nfft=SDR_NFFT):
        worst = max(worst, abs(row["error_db"]))
        print(
            f"  {row['what']:<28} expected {row['expected_db']:8.2f} dBFS/Hz   "
            f"measured {row['measured_db']:8.2f}   error {row['error_db']:+.3f} dB"
        )

    print()

    if worst > args.tolerance:
        sys.exit(f"FAIL: worst error {worst:.3f} dB exceeds the {args.tolerance:.2f} dB tolerance")

    print(f"PASS: worst error {worst:.3f} dB, within {args.tolerance:.2f} dB")


# ---------------------------------------------------------------
# Command line
# ---------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="RTL-SDR capture, wideband survey and analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # -- antenna measurement ------------------------------------------------

    cap = sub.add_parser("capture", help="sweep a band once and save it")
    cap.add_argument("--label", required=True, help="name for the saved sweep")
    cap.add_argument("--start", type=float, required=True, help="start frequency in Hz")
    cap.add_argument("--stop", type=float, required=True, help="stop frequency in Hz")
    cap.add_argument("--gain", type=float, default=SDR_GAIN_DB)
    cap.add_argument("--sample-rate", type=float, default=SDR_SAMPLE_RATE_HZ)
    cap.add_argument("--samples", type=int, default=SDR_SAMPLES_PER_STEP)
    cap.add_argument("--ppm", type=int, default=0)
    cap.add_argument("--overwrite", action="store_true")
    cap.set_defaults(func=cmd_capture)

    cmp_ = sub.add_parser("compare", help="difference an antenna sweep against a load sweep")
    cmp_.add_argument("--antenna", required=True, help="label of the antenna sweep")
    cmp_.add_argument("--reference", required=True, help="label of the 50 ohm load sweep")
    cmp_.add_argument("--vna", help="Touchstone .s1p to compare against")
    cmp_.add_argument("--smooth", type=int, default=SDR_SMOOTH_BINS)
    cmp_.add_argument("--search-lo", type=float, default=None, help="Hz")
    cmp_.add_argument("--search-hi", type=float, default=None, help="Hz")
    cmp_.add_argument("--min-rise", type=float, default=3.0,
                      help="dB of broadband noise rise below which the peak is reported "
                           "as receiver-noise-limited rather than as a result")
    cmp_.set_defaults(func=cmd_compare)

    # -- frequency axis check ----------------------------------------------
    ras = sub.add_parser("raster",
                         help="check the frequency axis against the UHF TV channel plan")
    ras.add_argument("--antenna", required=True)
    ras.add_argument("--reference", required=True)
    ras.add_argument("--smooth", type=int, default=SDR_SMOOTH_BINS)
    ras.add_argument("--threshold", type=float, default=15.0,
                     help="dB above the floor that counts as a multiplex")
    ras.set_defaults(func=cmd_raster)

    # -- wideband survey ----------------------------------------------------

    srv = sub.add_parser("survey", help="sweep the whole range repeatedly")
    srv.add_argument("--label", required=True, help="name for the saved survey")
    srv.add_argument("--start", type=float, default=SURVEY_START_HZ, help="Hz")
    srv.add_argument("--stop", type=float, default=SURVEY_STOP_HZ, help="Hz")
    srv.add_argument("--passes", type=int, default=SURVEY_PASSES)
    srv.add_argument("--bin", type=float, default=SURVEY_BIN_HZ, help="grid cell width in Hz")
    srv.add_argument("--gain", type=float, default=SDR_GAIN_DB)
    srv.add_argument("--sample-rate", type=float, default=SDR_SAMPLE_RATE_HZ)
    srv.add_argument("--samples", type=int, default=SURVEY_SAMPLES_PER_STEP)
    srv.add_argument("--occupancy-snr", type=float, default=SURVEY_OCCUPANCY_SNR_DB,
                     help="dB above the local floor for a bin to count as occupied in a pass")
    srv.add_argument("--ppm", type=int, default=0)
    srv.add_argument("--simulate", action="store_true",
                     help="use the simulated source instead of the radio")
    srv.add_argument("--seed", type=int, default=None, help="seed for the simulated source")
    srv.add_argument("--overwrite", action="store_true")
    srv.set_defaults(func=cmd_survey)

    rep = sub.add_parser("report", help="analyse a stored survey")
    rep.add_argument("--label", required=True)
    rep.add_argument("--snr", type=float, default=SURVEY_DETECT_SNR_DB,
                     help="dB above the local floor for a detection")
    rep.add_argument("--merge-gap", type=int, default=SURVEY_MERGE_GAP_BINS,
                     help="bins; runs closer than this are reported as one signal")
    rep.add_argument("--floor-window", type=float, default=SURVEY_FLOOR_WINDOW_HZ / 1e6,
                     help="width of the sliding noise floor window in MHz")
    rep.add_argument("--rows", type=int, default=SURVEY_OVERVIEW_ROWS,
                     help="panels in the overview figure")
    rep.add_argument("--zooms", type=int, default=SURVEY_ZOOM_BANDS,
                     help="bands given their own figure")
    rep.add_argument("--list", type=int, default=25, help="rows printed to the terminal")
    rep.set_defaults(func=cmd_report)

    # -- self test ----------------------------------------------------------

    test = sub.add_parser("selftest", help="check the analyser against a known source")
    test.add_argument("--sample-rate", type=float, default=SDR_SAMPLE_RATE_HZ)
    test.add_argument("--tolerance", type=float, default=0.5, help="dB")
    test.set_defaults(func=cmd_selftest)

    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.func(arguments)
