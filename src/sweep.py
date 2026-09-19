"""
Stitched wideband sweep for the SDR Spectrum Analyzer.

The RTL-SDR sees about 2.4 MHz at a time. A sweep from 100 MHz to 1 GHz is
therefore built by retuning across the band, keeping the usable middle of each
capture and concatenating the pieces.

Two parts of each capture are thrown away:

    the centre, where the DC offset and 1/f noise of the direct conversion
    receiver sit
    the edges, where the analogue anti-alias filter is already rolling off

Keeping the edges would put a periodic scallop into the stitched trace with the
period of the tuning step, which is easy to mistake for structure in the signal.
"""

import numpy as np

from .psd import welch_psd
from .sdr_source import capture_at, check_clipping


def plan_sweep(f_start_hz, f_stop_hz, sample_rate_hz, usable_fraction=0.8):
    """
    Choose the tuning steps for a stitched sweep.

    Parameters
    ----------
    f_start_hz, f_stop_hz : float
        Requested band.

    sample_rate_hz : float
        Sampling rate, equal to the captured bandwidth.

    usable_fraction : float
        Fraction of each capture kept, centred. 0.8 of 2.4 MHz leaves 1.92 MHz
        and stays clear of the filter skirts.

    Returns
    -------
    centers_hz : np.ndarray
        Centre frequency of each step.

    usable_bw_hz : float
        Width kept from each step.
    """

    if f_stop_hz <= f_start_hz:
        raise ValueError("f_stop must be above f_start")

    usable_bw_hz = sample_rate_hz * usable_fraction
    n_steps = int(np.ceil((f_stop_hz - f_start_hz) / usable_bw_hz))
    centers_hz = f_start_hz + usable_bw_hz * (np.arange(n_steps) + 0.5)

    return centers_hz, usable_bw_hz


def run_sweep(
    sdr,
    f_start_hz,
    f_stop_hz,
    sample_rate_hz,
    samples_per_step,
    nfft=4096,
    overlap=0.5,
    window_name="hann",
    usable_fraction=0.8,
    dc_exclusion_hz=40e3,
    settle_samples=65536,
    progress=True,
):
    """
    Sweep the band and return one stitched power spectral density trace.

    Parameters
    ----------
    sdr : RtlSdr
        Device from sdr_source.open_sdr().

    f_start_hz, f_stop_hz : float
        Band to cover.

    sample_rate_hz : float
        Must match what the device was opened with.

    samples_per_step : int
        Samples captured at each centre frequency. More samples means more
        Welch segments and a smoother noise floor: the standard deviation of
        the estimate falls as one over the square root of the segment count.

    nfft, overlap, window_name :
        Passed to welch_psd.

    usable_fraction, dc_exclusion_hz :
        Trimming, see the module docstring.

    settle_samples : int
        Discarded after each retune.

    progress : bool
        Print a line per step.

    Returns
    -------
    freqs_hz : np.ndarray
        Absolute frequency axis, ascending.

    psd_db : np.ndarray
        Stitched power spectral density in dBFS per Hz.

    meta : dict
        Sweep settings, resolution bandwidth and any clipping warnings.
    """

    centers_hz, usable_bw_hz = plan_sweep(
        f_start_hz, f_stop_hz, sample_rate_hz, usable_fraction
    )

    all_freqs = []
    all_psd = []
    warnings = []
    rbw_hz = None
    segments = None

    for i, fc in enumerate(centers_hz, 1):
        x = capture_at(sdr, fc, samples_per_step, settle_samples=settle_samples)

        clipped, fraction = check_clipping(x)
        if clipped:
            warnings.append(
                f"{fc/1e6:.1f} MHz: {100*fraction:.2f}% of samples at full scale, reduce gain"
            )

        f_rel, psd_db_step, info = welch_psd(
            x, sample_rate_hz, nfft=nfft, overlap=overlap, window_name=window_name
        )
        rbw_hz = info["rbw_hz"]
        segments = info["segments"]

        keep = (np.abs(f_rel) <= usable_bw_hz / 2) & (np.abs(f_rel) >= dc_exclusion_hz / 2)

        all_freqs.append(fc + f_rel[keep])
        all_psd.append(psd_db_step[keep])

        if progress:
            print(
                f"  step {i:3d}/{len(centers_hz)}  {fc/1e6:8.2f} MHz  "
                f"{segments} segments{'  CLIPPING' if clipped else ''}"
            )

    freqs_hz = np.concatenate(all_freqs)
    psd_db = np.concatenate(all_psd)

    order = np.argsort(freqs_hz)
    freqs_hz = freqs_hz[order]
    psd_db = psd_db[order]

    meta = {
        "f_start_hz": float(f_start_hz),
        "f_stop_hz": float(f_stop_hz),
        "sample_rate_hz": float(sample_rate_hz),
        "samples_per_step": int(samples_per_step),
        "n_steps": int(len(centers_hz)),
        "usable_bw_hz": float(usable_bw_hz),
        "dc_exclusion_hz": float(dc_exclusion_hz),
        "nfft": int(nfft),
        "overlap": float(overlap),
        "window": window_name,
        "rbw_hz": float(rbw_hz),
        "segments_per_step": int(segments),
        "warnings": warnings,
    }

    return freqs_hz, psd_db, meta


def save_sweep(path, freqs_hz, psd_db, meta, settings):
    """
    Write a sweep to a compressed .npz so it can be re-analysed without the
    hardware connected.
    """

    import json

    np.savez_compressed(
        path,
        freqs_hz=freqs_hz,
        psd_db=psd_db,
        meta_json=json.dumps(meta),
        settings_json=json.dumps(settings),
    )


def load_sweep(path):
    """
    Read a sweep written by save_sweep().

    Returns
    -------
    freqs_hz, psd_db, meta, settings
    """

    import json

    data = np.load(path, allow_pickle=False)

    return (
        data["freqs_hz"],
        data["psd_db"],
        json.loads(str(data["meta_json"])),
        json.loads(str(data["settings_json"])),
    )
