"""
Antenna match measured from the receiver's noise floor.

The idea
--------
A receiver's noise floor has two contributions: its own internal noise, and
external noise picked up by whatever is on its input. Replace the antenna with
a 50 ohm load and only the internal noise remains, because a matched resistor
at room temperature delivers far less power than the sky, the ground and the
city do through a resonant antenna.

Sweep the band twice, once with the antenna and once with the load, and
subtract. The difference is large where the antenna is well matched and
collapses to zero where it is not. The peak of the difference should therefore
line up with the minimum of |S11| measured on the VNA.

Why the subtraction matters
---------------------------
The RTL-SDR reports dBFS, not dBm. There is no absolute calibration, and the
gain of the tuner, the conversion loss and the ADC scaling are all unknown. But
they are the *same* unknown in both sweeps, as long as the gain is fixed and
the AGC is off, so they cancel in the difference. The measurement is designed
around the missing calibration rather than in spite of it.

What it does not give
---------------------
Antenna efficiency, absolute gain, or anything in dBm. The external noise level
also varies with frequency and with what is transmitting nearby, so the
difference trace is not a clean image of 1 - |S11|^2. It locates the resonance;
it does not replace the VNA.
"""

import numpy as np

from .psd import moving_average_db


def antenna_noise_delta(freqs_ant, psd_ant, freqs_ref, psd_ref, smooth_bins=201):
    """
    Difference between an antenna sweep and a reference load sweep.

    Parameters
    ----------
    freqs_ant, psd_ant : np.ndarray
        Sweep with the antenna connected.

    freqs_ref, psd_ref : np.ndarray
        Sweep with a 50 ohm load connected, same settings.

    smooth_bins : int
        Moving average length. The difference of two noisy traces is noisier
        than either, so this is usually a few hundred bins. It must stay well
        below the width of the antenna's match, which for these dipoles is 8 to
        12 percent of the centre frequency.

    Returns
    -------
    freqs_hz : np.ndarray
        Frequency axis of the antenna sweep.

    delta_db : np.ndarray
        Raw difference in dB.

    delta_smooth_db : np.ndarray
        Smoothed difference in dB.
    """

    psd_ref_on_ant_grid = np.interp(freqs_ant, freqs_ref, psd_ref)
    delta_db = psd_ant - psd_ref_on_ant_grid

    return freqs_ant, delta_db, moving_average_db(delta_db, smooth_bins)


def rise_summary(freqs_hz, delta_smooth_db, search_lo_hz=None, search_hi_hz=None,
                 usable_db=3.0):
    """
    Is there enough noise rise for the peak to mean anything?

    The method only works while external noise dominates the receiver's own.
    When it does not, the difference of the two sweeps is flat at 0 dB except
    at transmitters, and argmax returns the strongest broadcast carrier in the
    search range rather than the antenna's resonance. That failure is silent:
    the number looks like a measurement.

    A robust measure of the broadband rise is the median of the smoothed
    difference, which ignores the transmitters. Below usable_db the trace is
    receiver-noise-limited and the peak should not be quoted.

    Returns
    -------
    dict with median_db, p90_db, usable (bool) and usable_db.
    """

    mask = np.ones(len(freqs_hz), dtype=bool)
    if search_lo_hz is not None:
        mask &= freqs_hz >= search_lo_hz
    if search_hi_hz is not None:
        mask &= freqs_hz <= search_hi_hz

    d = delta_smooth_db[mask]
    median = float(np.median(d))

    return dict(median_db=median, p90_db=float(np.percentile(d, 90)),
                usable=median >= usable_db, usable_db=usable_db)


def find_peak(freqs_hz, delta_smooth_db, search_lo_hz=None, search_hi_hz=None):
    """
    Frequency of maximum noise rise, with a parabolic refinement.

    Uses the same three point interpolation as the dipole project's
    analyze_dipole.py, so the two measurements are read the same way.

    Parameters
    ----------
    freqs_hz, delta_smooth_db : np.ndarray
        Output of antenna_noise_delta().

    search_lo_hz, search_hi_hz : float or None
        Restrict the search. Without a restriction the peak may land on a
        strong broadcast transmitter rather than on the antenna's resonance.

    Returns
    -------
    f_peak_hz : float

    level_db : float
    """

    mask = np.ones(len(freqs_hz), dtype=bool)
    if search_lo_hz is not None:
        mask &= freqs_hz >= search_lo_hz
    if search_hi_hz is not None:
        mask &= freqs_hz <= search_hi_hz

    if not mask.any():
        raise ValueError("search range contains no points")

    f = freqs_hz[mask]
    d = delta_smooth_db[mask]
    i = int(np.argmax(d))

    if 0 < i < len(f) - 1:
        y0, y1, y2 = d[i - 1:i + 2]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            return f[i] + 0.5 * (y0 - y2) / denom * (f[i + 1] - f[i]), float(d[i])

    return float(f[i]), float(d[i])


def load_s1p_min(path):
    """
    Frequency of minimum |S11| in a Touchstone file, for comparison.

    Returns (frequency in Hz, level in dB), or None if scikit-rf is missing.
    """

    try:
        import skrf as rf
    except ImportError:
        return None

    network = rf.Network(path)
    f = network.f
    s_db = 20 * np.log10(np.abs(network.s[:, 0, 0]))
    i = int(np.argmin(s_db))

    if 0 < i < len(f) - 1:
        y0, y1, y2 = s_db[i - 1:i + 2]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            return f[i] + 0.5 * (y0 - y2) / denom * (f[i + 1] - f[i]), float(s_db[i])

    return float(f[i]), float(s_db[i])


# -----------------------------
# Plot
# -----------------------------

def plot_antenna_noise(
    freqs_hz,
    psd_ant_db,
    psd_ref_db,
    delta_db,
    delta_smooth_db,
    save_path,
    f_peak_hz=None,
    f_vna_hz=None,
    title=None,
):
    """
    Two panels: the two noise floors, and their difference.

    Parameters
    ----------
    freqs_hz : np.ndarray
        Frequency axis.

    psd_ant_db, psd_ref_db : np.ndarray
        The two sweeps, on the same axis.

    delta_db, delta_smooth_db : np.ndarray
        Raw and smoothed difference.

    save_path : pathlib.Path or str

    f_peak_hz : float or None
        Marked on the lower panel.

    f_vna_hz : float or None
        |S11| minimum from the VNA, marked for comparison.

    title : str or None
    """

    import matplotlib.pyplot as plt

    f_mhz = freqs_hz / 1e6
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)

    ax_top.plot(f_mhz, psd_ant_db, lw=0.6, alpha=0.85, label="antenna")
    ax_top.plot(f_mhz, psd_ref_db, lw=0.6, alpha=0.85, label="50 ohm load")
    ax_top.set_ylabel("PSD (dBFS/Hz)")
    ax_top.grid(alpha=0.3)
    ax_top.legend(fontsize=8)
    if title:
        ax_top.set_title(title)

    ax_bot.plot(f_mhz, delta_db, lw=0.5, alpha=0.3, color="C0", label="difference")
    ax_bot.plot(f_mhz, delta_smooth_db, lw=2, color="C0", label="smoothed")
    ax_bot.axhline(0, color="grey", ls=":", lw=1)

    if f_peak_hz is not None:
        ax_bot.axvline(
            f_peak_hz / 1e6, color="C0", ls="--", lw=1.2,
            label=f"noise peak {f_peak_hz/1e6:.1f} MHz",
        )
    if f_vna_hz is not None:
        ax_bot.axvline(
            f_vna_hz / 1e6, color="C3", ls="--", lw=1.2,
            label=f"VNA |S11| min {f_vna_hz/1e6:.1f} MHz",
        )

    ax_bot.set_xlabel("Frequency (MHz)")
    ax_bot.set_ylabel("Antenna minus load (dB)")
    ax_bot.grid(alpha=0.3)
    ax_bot.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)
