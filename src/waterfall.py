"""
Plots for a wideband survey.

Three figures, each answering a different question.

The overview answers what is on the air. It is split into several stacked
panels because 1.7 GHz on one axis is about two megahertz per pixel, which
hides everything narrower than a DAB block. Behind the traces the band plan
is shaded, so a detection is read against the service that is supposed to be
there rather than against a bare frequency.

The waterfall answers what changes. Frequency runs across, passes run down,
and colour is the largest power seen in each cell during that pass. A carrier
that is always on draws a continuous vertical line; a handset or a TETRA
channel draws dashes; a broadcast transmitter and a spur look identical here,
which is why the spur has to be identified some other way.

The band plot answers what is inside one service, with the occupancy of every
bin underneath it on its own axis rather than squeezed onto a second scale.

Levels are dBFS per Hz throughout. The receiver has no absolute calibration,
so these numbers compare within a survey and against another survey taken at
the same gain, and are not dBm.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from . import bands

# ---------------------------------------------------------------
# Palette
# ---------------------------------------------------------------
# Three categorical slots, taken in fixed order, plus recessive ink for
# everything that is context rather than data. Validated for colour vision
# deficiency as a set; the detection colour is the one that sits below 3:1
# against the surface, which is why every detection also carries a text label
# rather than relying on its colour.

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d9d8d3"

COLOR_MAXHOLD = "#2a78d6"   # slot 1, blue
COLOR_MEAN = "#eb6834"      # slot 2, orange
COLOR_MARK = "#1baf7a"      # slot 3, aqua
COLOR_REFERENCE = "#8a8983"

# Sequential one hue ramp for the waterfall: light is near the noise floor,
# dark is a strong signal.
WATERFALL_CMAP = LinearSegmentedColormap.from_list(
    "survey_blue",
    ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
)

BAND_SHADE = "#eceae4"


def _style(ax):
    """Recessive axes: the trace is the figure, the frame is not."""

    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.9)
    ax.set_axisbelow(True)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    ax.tick_params(colors=INK_MUTED, labelsize=8)


def _shade_bands(ax, f_lo_hz, f_hi_hz, label_fraction=0.02):
    """
    Shade the band plan behind the trace and label the wider entries.

    Only entries covering more than label_fraction of the visible width get a
    name, otherwise the labels overlap into a grey smear. The narrow ones are
    still shaded, and the detection list names them individually.
    """

    span = f_hi_hz - f_lo_hz
    shown = bands.bands_in_range(f_lo_hz, f_hi_hz)

    for index, entry in enumerate(shown):
        lo = max(entry["f_lo_hz"], f_lo_hz)
        hi = min(entry["f_hi_hz"], f_hi_hz)

        if hi <= lo:
            continue

        ax.axvspan(lo / 1e6, hi / 1e6, color=BAND_SHADE, alpha=0.55 if index % 2 else 0.85, lw=0)

        if (hi - lo) / span >= label_fraction:
            # Against the left edge of the band rather than its centre, which
            # is where the tallest signal in a band tends to be and therefore
            # where its own label wants to go.
            ax.text(
                (lo + 0.012 * span) / 1e6,
                0.97,
                entry["short"],
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=6.5,
                color=INK_MUTED,
                rotation=90,
            )


def _headroom_ylim(ax, trace_db, marker_height=0.66, floor_percentile=0.5, floor_margin=4.0):
    """
    Leave the top third of the panel free for text.

    The band names run vertically along the top of the axes and the detection
    labels sit just above their markers, so an autoscaled panel puts the two
    on top of each other. Placing the tallest marker at a fixed fraction of
    the panel height reserves the rest for them whatever the dynamic range of
    the band happens to be.
    """

    finite = trace_db[np.isfinite(trace_db)]

    if finite.size == 0:
        return

    bottom = float(np.percentile(finite, floor_percentile)) - floor_margin
    top = float(finite.max())
    span = max(top - bottom, 1.0)

    ax.set_ylim(bottom, bottom + span / marker_height)


def _label_detections(ax, detections, f_lo_hz, f_hi_hz, max_labels=8, min_separation=0.05):
    """
    Mark every detection, name only the strongest few.

    Labels are placed strongest first and one is skipped whenever it would
    land within min_separation of the panel width from a label already
    placed, which is what keeps a busy band such as FM from turning into a
    stack of overlapping text. Every detection still gets its marker, and the
    full list is in the CSV.
    """

    inside = [d for d in detections if f_lo_hz <= d["f_peak_hz"] <= f_hi_hz]

    if not inside:
        return

    ax.plot(
        [d["f_peak_hz"] / 1e6 for d in inside],
        [d["peak_db"] for d in inside],
        linestyle="none",
        marker="o",
        markersize=4.5,
        markerfacecolor="none",
        markeredgecolor=COLOR_MARK,
        markeredgewidth=1.2,
        label="detected signal",
        zorder=5,
    )

    span_mhz = (f_hi_hz - f_lo_hz) / 1e6
    placed = []

    for entry in sorted(inside, key=lambda d: d["peak_db"], reverse=True):
        if len(placed) >= max_labels:
            break

        f_mhz = entry["f_peak_hz"] / 1e6
        if any(abs(f_mhz - other) < min_separation * span_mhz for other in placed):
            continue

        placed.append(f_mhz)

        ax.annotate(
            f"{f_mhz:.2f}\n{entry['band']}",
            xy=(f_mhz, entry["peak_db"]),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=6.5,
            color=INK,
            zorder=6,
        )


def plot_overview(
    centers_hz,
    mean_db,
    maxhold_db,
    floor_db,
    detections,
    save_path,
    rows=4,
    title="Wideband survey",
    subtitle=None,
    simulated=False,
):
    """
    The whole surveyed range, split over several panels.

    Parameters
    ----------
    centers_hz : np.ndarray
        Survey grid.

    mean_db, maxhold_db, floor_db : np.ndarray
        Traces on that grid.

    detections : list of dict
        From signals.detect_signals().

    save_path : pathlib.Path or str

    rows : int
        Panels. Each covers an equal slice of the range.

    title, subtitle : str
        Headings.

    simulated : bool
        Marks the figure as coming from the simulated source.
    """

    f_lo, f_hi = float(centers_hz[0]), float(centers_hz[-1])
    edges = np.linspace(f_lo, f_hi, rows + 1)

    header_inches = 0.75
    fig, axes = plt.subplots(
        rows, 1, figsize=(11, 2.35 * rows + header_inches), facecolor=SURFACE
    )
    axes = np.atleast_1d(axes)

    for index, ax in enumerate(axes):
        lo, hi = edges[index], edges[index + 1]
        window = (centers_hz >= lo) & (centers_hz <= hi)

        _style(ax)
        _shade_bands(ax, lo, hi)

        ax.plot(centers_hz[window] / 1e6, mean_db[window], lw=0.55, color=COLOR_MEAN,
                alpha=0.9, label="mean")
        ax.plot(centers_hz[window] / 1e6, maxhold_db[window], lw=0.7, color=COLOR_MAXHOLD,
                label="max hold")
        ax.plot(centers_hz[window] / 1e6, floor_db[window], lw=0.9, color=COLOR_REFERENCE,
                ls="--", label="noise floor")

        _label_detections(ax, detections, lo, hi)

        ax.set_xlim(lo / 1e6, hi / 1e6)
        ax.set_ylabel("dBFS/Hz", fontsize=8, color=INK_MUTED)

        _headroom_ylim(ax, maxhold_db[window])

        if index == 0:
            handles, labels = ax.get_legend_handles_labels()
            unique = dict(zip(labels, handles))
            ax.legend(unique.values(), unique.keys(), loc="upper right", fontsize=7,
                      framealpha=0.9, edgecolor=GRID)

    axes[-1].set_xlabel("Frequency (MHz)", fontsize=9, color=INK_MUTED)

    height_inches = fig.get_size_inches()[1]
    reserved = 1.0 - header_inches / height_inches

    fig.tight_layout(rect=(0, 0, 1, reserved))

    heading = f"SIMULATED — {title}" if simulated else title
    fig.text(0.012, 1.0 - 0.28 * header_inches / height_inches, heading,
             fontsize=13, color=INK, ha="left", va="center")

    if subtitle:
        fig.text(0.012, 1.0 - 0.62 * header_inches / height_inches, subtitle,
                 fontsize=8, color=INK_MUTED, ha="left", va="center")
    fig.savefig(save_path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def plot_waterfall(
    waterfall_freqs_hz,
    waterfall_db,
    save_path,
    pass_times=None,
    title="Survey waterfall",
    simulated=False,
    clip_percentiles=(5, 99.5),
):
    """
    Passes against frequency, coloured by the peak power in each cell.

    Parameters
    ----------
    waterfall_freqs_hz : np.ndarray
        Decimated frequency axis from the survey.

    waterfall_db : np.ndarray
        One row per pass.

    save_path : pathlib.Path or str

    pass_times : np.ndarray or None
        Unix time at the start of each pass; used to label the axis in
        minutes from the start of the survey.

    title : str

    simulated : bool

    clip_percentiles : tuple
        Colour scale limits, as percentiles of the data. Clipping the top
        keeps one strong broadcast carrier from compressing everything else
        into the lightest step.
    """

    if waterfall_db.ndim != 2 or waterfall_db.shape[0] < 2:
        return False

    finite = waterfall_db[np.isfinite(waterfall_db)]
    if finite.size == 0:
        return False

    vmin, vmax = np.percentile(finite, clip_percentiles)

    if pass_times is not None and len(pass_times) == waterfall_db.shape[0]:
        # pass_times holds the start of each pass, so the bottom of the image
        # is one pass beyond the last start, not at it.
        minutes = (np.asarray(pass_times) - pass_times[0]) / 60.0
        step = float(np.median(np.diff(minutes))) if len(minutes) > 1 else 1.0
        y_hi = float(minutes[-1]) + step
        y_label = "Minutes from the start of the survey"
    else:
        y_hi = float(waterfall_db.shape[0])
        y_label = "Pass"

    fig, ax = plt.subplots(figsize=(11, 5.0), facecolor=SURFACE)
    _style(ax)
    ax.grid(False)

    image = ax.imshow(
        waterfall_db,
        aspect="auto",
        origin="upper",
        interpolation="nearest",
        cmap=WATERFALL_CMAP,
        vmin=vmin,
        vmax=vmax,
        extent=[
            waterfall_freqs_hz[0] / 1e6,
            waterfall_freqs_hz[-1] / 1e6,
            y_hi,
            0.0,
        ],
    )

    bar = fig.colorbar(image, ax=ax, pad=0.015)
    bar.set_label("Peak power in the cell (dBFS/Hz)", fontsize=8, color=INK_MUTED)
    bar.ax.tick_params(labelsize=7, colors=INK_MUTED)
    bar.outline.set_edgecolor(GRID)

    ax.set_xlabel("Frequency (MHz)", fontsize=9, color=INK_MUTED)
    ax.set_ylabel(y_label, fontsize=9, color=INK_MUTED)
    ax.set_title(f"SIMULATED — {title}" if simulated else title, fontsize=12, color=INK, loc="left")

    fig.tight_layout()
    fig.savefig(save_path, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    return True


def plot_band(
    centers_hz,
    mean_db,
    maxhold_db,
    floor_db,
    occupancy,
    detections,
    f_lo_hz,
    f_hi_hz,
    save_path,
    title=None,
    snr_db=None,
    simulated=False,
):
    """
    One band, with occupancy on its own panel underneath.

    The occupancy is kept on a separate axis rather than a second y scale on
    the same panel. Two scales sharing one frame invite the reader to compare
    the heights of two curves that have nothing to do with each other.

    Parameters
    ----------
    centers_hz, mean_db, maxhold_db, floor_db, occupancy : np.ndarray
        Survey grid and traces.

    detections : list of dict

    f_lo_hz, f_hi_hz : float
        Range to show.

    save_path : pathlib.Path or str

    title : str or None

    snr_db : float or None
        Detection threshold, drawn above the floor if given.

    simulated : bool
    """

    window = (centers_hz >= f_lo_hz) & (centers_hz <= f_hi_hz)

    if not window.any():
        return False

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10, 6.0), height_ratios=(3, 1), sharex=True, facecolor=SURFACE
    )

    _style(ax_top)
    _style(ax_bot)
    _shade_bands(ax_top, f_lo_hz, f_hi_hz, label_fraction=0.06)

    f_mhz = centers_hz[window] / 1e6

    ax_top.plot(f_mhz, mean_db[window], lw=0.7, color=COLOR_MEAN, alpha=0.9, label="mean")
    ax_top.plot(f_mhz, maxhold_db[window], lw=1.0, color=COLOR_MAXHOLD, label="max hold")
    ax_top.plot(f_mhz, floor_db[window], lw=1.0, color=COLOR_REFERENCE, ls="--", label="noise floor")

    if snr_db is not None:
        ax_top.plot(f_mhz, floor_db[window] + snr_db, lw=0.9, color=COLOR_REFERENCE, ls=":",
                    label=f"threshold, floor + {snr_db:.0f} dB")

    _label_detections(ax_top, detections, f_lo_hz, f_hi_hz, max_labels=12,
                      min_separation=0.035)

    _headroom_ylim(ax_top, maxhold_db[window], marker_height=0.62)

    ax_top.set_ylabel("dBFS/Hz", fontsize=9, color=INK_MUTED)
    ax_top.legend(loc="upper right", fontsize=7, framealpha=0.9, edgecolor=GRID)

    heading = title or f"{f_lo_hz/1e6:.1f} to {f_hi_hz/1e6:.1f} MHz"
    ax_top.set_title(f"SIMULATED — {heading}" if simulated else heading,
                     fontsize=12, color=INK, loc="left")

    ax_bot.fill_between(f_mhz, 0, 100 * occupancy[window], step="mid",
                        color=COLOR_MAXHOLD, alpha=0.35, lw=0)
    ax_bot.plot(f_mhz, 100 * occupancy[window], drawstyle="steps-mid", lw=0.9, color=COLOR_MAXHOLD)
    ax_bot.set_ylim(0, 105)
    ax_bot.set_ylabel("Occupancy (%)", fontsize=9, color=INK_MUTED)
    ax_bot.set_xlabel("Frequency (MHz)", fontsize=9, color=INK_MUTED)
    ax_bot.set_xlim(f_lo_hz / 1e6, f_hi_hz / 1e6)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    return True


def busiest_bands(detections, limit=6, min_count=1):
    """
    Bands worth a zoom plot, strongest first.

    Ranked by the strongest detection in the band rather than by how many
    detections it holds, because a single loud transmitter is more worth
    looking at than a dozen marginal ones.

    Returns
    -------
    list of dict
        Keys short, f_lo_hz, f_hi_hz, count, peak_db.
    """

    groups = {}

    for entry in detections:
        band = bands.lookup(entry["f_peak_hz"])
        if band is None:
            continue

        key = (band["short"], band["f_lo_hz"], band["f_hi_hz"])
        groups.setdefault(key, []).append(entry)

    rows = [
        {
            "short": short,
            "f_lo_hz": lo,
            "f_hi_hz": hi,
            "count": len(members),
            "peak_db": max(m["peak_db"] for m in members),
        }
        for (short, lo, hi), members in groups.items()
        if len(members) >= min_count
    ]

    rows.sort(key=lambda row: row["peak_db"], reverse=True)

    return rows[:limit]
