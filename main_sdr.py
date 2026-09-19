"""
Real hardware entry point for the SDR Spectrum Analyzer.

main.py analyses a synthetic signal. This script captures from an RTL-SDR and
uses the analyser to measure something unknown: where a dipole is resonant,
read from the rise in the receiver's noise floor when the antenna is connected
instead of a 50 ohm load.

Two subcommands, deliberately separate so the analysis can be repeated and
re-plotted later without the hardware present.

    python main_sdr.py capture --label A_antenna --start 100e6 --stop 190e6
    python main_sdr.py capture --label A_load    --start 100e6 --stop 190e6
    python main_sdr.py compare --antenna A_antenna --reference A_load \
                               --vna ../dipole-sim-vs-measurement/data/A_meas_1.s1p

Both captures of a pair must use the same gain, sample rate and band. The
script records the settings in each .npz and refuses to compare mismatched
pairs.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from src.antenna_noise import (
    antenna_noise_delta,
    find_peak,
    load_s1p_min,
    plot_antenna_noise,
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
)
from src.sdr_source import open_sdr
from src.sweep import load_sweep, run_sweep, save_sweep

SWEEP_DIR = Path("outputs") / "sweeps"
FIGURE_DIR = Path("outputs") / "figures"


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

    print(f"Gain (both sweeps):       {set_ant['gain_db_applied']:.1f} dB")
    print(f"Resolution bandwidth:     {meta_ant['rbw_hz']/1e3:.2f} kHz")
    print(f"Smoothing:                {args.smooth} bins")
    print(f"Noise rise peak:          {f_peak_hz/1e6:.2f} MHz at {level_db:.2f} dB")
    print(f"Median rise across band:  {np.median(delta_smooth_db):.2f} dB")

    f_vna_hz = None
    if args.vna:
        result = load_s1p_min(args.vna)
        if result is None:
            print("scikit-rf not installed, skipping the VNA comparison")
        else:
            f_vna_hz, s11_db = result
            error_pct = 100 * (f_peak_hz - f_vna_hz) / f_vna_hz
            print(f"VNA |S11| minimum:        {f_vna_hz/1e6:.2f} MHz at {s11_db:.2f} dB")
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


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = parser.add_subparsers(dest="command", required=True)

    cap = sub.add_parser("capture", help="sweep a band and save it")
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
    cmp_.set_defaults(func=cmd_compare)

    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.func(arguments)
