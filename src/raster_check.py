"""
Verify the swept spectrum against a known broadcast raster.

A spectrum analyser built out of an RTL-SDR and a Welch estimator has two
things that can be silently wrong: the frequency axis, and the width it
assigns to a signal. Both are easy to check without any test equipment,
because terrestrial television is a comb of channels on an exact grid.

In the ITU Region 1 UHF plan a DVB-T channel of index n is centred at

    f_c = 306 + 8n   MHz

and the 8 MHz slot carries an occupied bandwidth of 7.61 MHz. Finding those
blocks in a sweep and comparing their measured centres and widths against the
plan is an absolute check: the numbers come from the broadcaster's licence,
not from anything in this repository.

**The comb stops at channel 48.** The band above 694 MHz was taken away from
television in the second digital dividend and sold to the mobile operators;
in Greece it was auctioned in December 2020. 694 to 790 MHz is now 3GPP band
n28, whose downlink occupies 758 to 788 MHz in carriers that are typically
10 MHz wide and therefore about 9 MHz of occupied bandwidth. Those carriers
sit close enough to the old channel positions that extending the comb to
channel 60 labels them as channels 57 and 58 and then reports them as
television that is 1 to 3 MHz off frequency.

Blocks outside 470 to 694 MHz are therefore identified through the band plan
in `bands.py` rather than forced onto the television grid, and they take no
part in the frequency statistics. A mobile carrier is not a worse reference
than a television channel in principle — LTE centres sit on a 100 kHz raster
— but which carriers an operator has lit is not known in advance, while the
television comb is published.

The blocks are found in the antenna-minus-load difference rather than in the
antenna sweep, because the difference already has the receiver's own shape
divided out and a DVB-T multiplex stands 20 dB above the floor in it.
"""

import numpy as np

from . import bands

UHF_OFFSET_MHZ = 306.0
UHF_SPACING_MHZ = 8.0
DVBT_OCCUPIED_MHZ = 7.61

# Lowest and highest UHF television channel still carrying DVB-T in Region 1,
# and the band edges they imply.
UHF_CHANNEL_MIN = 21
UHF_CHANNEL_MAX = 48
UHF_BAND_HZ = (470e6, 694e6)


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
    Match each block to the UHF television plan and report the error.

    A block counts towards the frequency statistics only if it is inside the
    television band, lands on a channel between UHF_CHANNEL_MIN and
    UHF_CHANNEL_MAX, and has an occupied bandwidth within width_tolerance_hz
    of the DVB-T figure. Anything else is reported with the reason it was not
    counted: a wrong width means the block is not one clean multiplex, and a
    centre outside the band means it is not television at all.

    Returns
    -------
    rows : list of dict
        channel, raster_hz, centre_hz, width_hz, error_hz, ppm, clean, note.
        channel, raster_hz, error_hz and ppm are None for a block that is not
        on the television grid.

    summary : dict
        n_clean, n_other, mean_ppm, spread_ppm, mean_width_hz
        (NaN when nothing is clean)
    """

    rows = []

    for lo, hi in blocks:
        centre = 0.5 * (lo + hi)
        width = hi - lo

        if not UHF_BAND_HZ[0] <= centre <= UHF_BAND_HZ[1]:
            entry = bands.lookup(centre)
            where = f"{entry['short']}, {entry['description']}" if entry else "outside the band plan"
            rows.append(dict(channel=None, raster_hz=None, centre_hz=centre, width_hz=width,
                             error_hz=None, ppm=None, clean=False,
                             note=f"not television: {where}"))
            continue

        n = int(round((centre / 1e6 - UHF_OFFSET_MHZ) / UHF_SPACING_MHZ))
        raster = (UHF_OFFSET_MHZ + UHF_SPACING_MHZ * n) * 1e6

        if not UHF_CHANNEL_MIN <= n <= UHF_CHANNEL_MAX:
            rows.append(dict(channel=n, raster_hz=raster, centre_hz=centre, width_hz=width,
                             error_hz=centre - raster, ppm=1e6 * (centre - raster) / raster,
                             clean=False,
                             note=f"channel {n} is outside the {UHF_CHANNEL_MIN}-{UHF_CHANNEL_MAX} plan"))
            continue

        clean = abs(width - DVBT_OCCUPIED_MHZ * 1e6) <= width_tolerance_hz
        rows.append(dict(channel=n, raster_hz=raster, centre_hz=centre, width_hz=width,
                         error_hz=centre - raster, ppm=1e6 * (centre - raster) / raster,
                         clean=clean,
                         note="" if clean else "width is not one multiplex, not counted"))

    good = [r for r in rows if r["clean"]]

    if good:
        ppm = np.array([r["ppm"] for r in good])
        summary = dict(n_clean=len(good), n_other=len(rows) - len(good),
                       mean_ppm=float(ppm.mean()),
                       spread_ppm=float(ppm.std(ddof=1)) if len(good) > 1 else 0.0,
                       mean_width_hz=float(np.mean([r["width_hz"] for r in good])))
    else:
        summary = dict(n_clean=0, n_other=len(rows), mean_ppm=np.nan,
                       spread_ppm=np.nan, mean_width_hz=np.nan)

    return rows, summary


def format_report(rows, summary):
    """
    Plain text report, one line per block.

    The last paragraph is the part that matters and the one easiest to get
    wrong. A receiver whose reference oscillator is off produces the *same*
    fractional error on every channel, so it shows up as a mean offset with a
    small spread. An edge-finding limit produces an error that changes block
    to block, so it shows up as a spread with a mean near zero. Quoting the
    mean alone, and calling it the crystal error, confuses the two.
    """

    out = [f"{'channel':>8} {'raster MHz':>11} {'measured MHz':>13} "
           f"{'width MHz':>10} {'error kHz':>10} {'ppm':>8}  flag",
           "-" * 78]

    for r in rows:
        channel = f"{r['channel']:8d}" if r["channel"] is not None else f"{'-':>8}"
        raster = f"{r['raster_hz']/1e6:11.0f}" if r["raster_hz"] is not None else f"{'-':>11}"
        error = f"{r['error_hz']/1e3:10.1f}" if r["error_hz"] is not None else f"{'-':>10}"
        ppm = f"{r['ppm']:8.1f}" if r["ppm"] is not None else f"{'-':>8}"
        out.append(f"{channel} {raster} {r['centre_hz']/1e6:13.3f} "
                   f"{r['width_hz']/1e6:10.2f} {error} {ppm}  {r['note']}")

    if not summary["n_clean"]:
        out += ["", "no clean multiplex in this sweep, nothing to check the axis against"]
        return "\n".join(out)

    mean, spread = summary["mean_ppm"], summary["spread_ppm"]

    out += ["",
            f"{summary['n_clean']} clean multiplexes, {summary['n_other']} blocks not counted",
            f"mean occupied bandwidth {summary['mean_width_hz']/1e6:.2f} MHz "
            f"against the DVB-T figure of {DVBT_OCCUPIED_MHZ:.2f} MHz",
            f"frequency error {mean:+.1f} ppm, spread {spread:.1f} ppm"]

    if summary["n_clean"] >= 3:
        if spread > 0 and abs(mean) < spread:
            out.append("the spread is larger than the mean, so this is dominated by how well "
                       "the block edges can be located, not by the receiver's reference")
        elif spread > 0:
            out.append("the mean is larger than the spread, which is the signature of a "
                       f"reference oscillator error; try --ppm {mean:+.0f}")

    return "\n".join(out)
