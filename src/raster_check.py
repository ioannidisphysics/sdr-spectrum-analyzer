"""
Verify the swept spectrum against a known broadcast raster.

A spectrum analyser built out of an RTL-SDR and a Welch estimator has two
things that can be silently wrong: the frequency axis, and the width it
assigns to a signal. Both are easy to check without any test equipment,
because terrestrial television is a comb of channels on an exact grid.

In the ITU Region 1 UHF plan a DVB-T channel of index n is centred at

    f_c = 306 + 8n   MHz,     n = 21 .. 60

and the 8 MHz slot carries an occupied bandwidth of 7.61 MHz. Finding those
blocks in a sweep and comparing their measured centres and widths against the
plan is an absolute check: the numbers come from the broadcaster's licence,
not from anything in this repository.

The blocks are found in the antenna-minus-load difference rather than in the
antenna sweep, because the difference already has the receiver's own shape
divided out and a DVB-T multiplex stands 20 dB above the floor in it.
"""

import numpy as np

UHF_OFFSET_MHZ = 306.0
UHF_SPACING_MHZ = 8.0
DVBT_OCCUPIED_MHZ = 7.61


def find_blocks(freqs_hz, delta_smooth_db, threshold_db=15.0, min_width_hz=6e6):
    """
    Contiguous runs that stand threshold_db above the noise floor.

    Parameters
    ----------
    freqs_hz, delta_smooth_db : np.ndarray
        Smoothed antenna-minus-load difference.

    threshold_db : float
        A DVB-T multiplex is 20 to 30 dB up. 15 dB keeps the weaker ones
        without picking up narrowband carriers.

    min_width_hz : float
        Blocks narrower than this are partial: the multiplex ran off the end
        of the sweep, or two channels merged. They are dropped rather than
        reported, because their centre is not their centre.

    Returns
    -------
    list of (f_lo_hz, f_hi_hz)
    """

    above = delta_smooth_db > threshold_db
    runs, start = [], None
    for i, v in enumerate(above):
        if v and start is None:
            start = i
        elif not v and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(above) - 1))

    return [(freqs_hz[a], freqs_hz[b]) for a, b in runs
            if freqs_hz[b] - freqs_hz[a] >= min_width_hz]


def check_uhf_raster(blocks, width_tolerance_hz=3e5):
    """
    Match each block to the nearest UHF channel and report the error.

    Blocks whose width is not within width_tolerance_hz of the DVB-T occupied
    bandwidth are marked but not used for the frequency statistics: a wrong
    width means the block is not one clean multiplex, so its centre carries no
    information about the frequency axis.

    Returns
    -------
    rows : list of dict
        channel, raster_hz, centre_hz, width_hz, error_hz, ppm, clean

    summary : dict
        n_clean, mean_ppm, spread_ppm, mean_width_hz  (NaN when nothing is clean)
    """

    rows = []
    for lo, hi in blocks:
        centre = 0.5 * (lo + hi)
        width = hi - lo
        n = int(round((centre / 1e6 - UHF_OFFSET_MHZ) / UHF_SPACING_MHZ))
        raster = (UHF_OFFSET_MHZ + UHF_SPACING_MHZ * n) * 1e6
        clean = abs(width - DVBT_OCCUPIED_MHZ * 1e6) <= width_tolerance_hz
        rows.append(dict(channel=n, raster_hz=raster, centre_hz=centre, width_hz=width,
                         error_hz=centre - raster, ppm=1e6 * (centre - raster) / raster,
                         clean=clean))

    good = [r for r in rows if r["clean"]]
    if good:
        ppm = np.array([r["ppm"] for r in good])
        summary = dict(n_clean=len(good), mean_ppm=float(ppm.mean()),
                       spread_ppm=float(ppm.std(ddof=1)) if len(good) > 1 else 0.0,
                       mean_width_hz=float(np.mean([r["width_hz"] for r in good])))
    else:
        summary = dict(n_clean=0, mean_ppm=np.nan, spread_ppm=np.nan,
                       mean_width_hz=np.nan)

    return rows, summary


def format_report(rows, summary):
    """Plain text report, one line per block."""

    out = [f"{'channel':>8} {'raster MHz':>11} {'measured MHz':>13} "
           f"{'width MHz':>10} {'error kHz':>10} {'ppm':>8}  flag",
           "-" * 72]
    for r in rows:
        out.append(f"{r['channel']:8d} {r['raster_hz']/1e6:11.0f} {r['centre_hz']/1e6:13.3f} "
                   f"{r['width_hz']/1e6:10.2f} {r['error_hz']/1e3:10.1f} {r['ppm']:8.1f}"
                   f"  {'' if r['clean'] else 'partial block, not counted'}")

    if summary["n_clean"]:
        out += ["",
                f"{summary['n_clean']} clean multiplexes",
                f"mean occupied bandwidth {summary['mean_width_hz']/1e6:.2f} MHz "
                f"against the DVB-T figure of {DVBT_OCCUPIED_MHZ:.2f} MHz",
                f"frequency error {summary['mean_ppm']:+.1f} ppm, "
                f"spread {summary['spread_ppm']:.1f} ppm"]
    else:
        out += ["", "no clean multiplex in this sweep, nothing to check the axis against"]

    return "\n".join(out)
