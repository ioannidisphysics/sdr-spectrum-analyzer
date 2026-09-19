"""
RTL-SDR capture for the SDR Spectrum Analyzer.

Provides real IQ samples with the same shape of contract as
synthetic.generate_synthetic_iq(), so the rest of the toolkit does not care
where the samples came from.

Hardware notes
--------------
The RTL-SDR Blog V4 uses the R828D tuner and covers roughly 24 MHz to 1766 MHz.
librtlsdr must be installed and reachable:

    pip install pyrtlsdr

On Windows the driver has to be replaced with WinUSB using Zadig, and
rtlsdr.dll must be on the PATH. The RTL-SDR Blog release bundle contains both.

Two settings matter more than anything else for a measurement that will be
compared against another measurement:

    the gain must be manual and identical between runs
    the AGC of the RTL2832 must be off

Otherwise the receiver silently rescales itself between the antenna run and the
reference run, and the difference of the two says nothing about the antenna.
"""

import numpy as np

V4_FREQ_MIN_HZ = 24e6
V4_FREQ_MAX_HZ = 1766e6


def open_sdr(sample_rate_hz, gain_db, freq_correction_ppm=0):
    """
    Open the first RTL-SDR and put it in a fixed, repeatable state.

    Parameters
    ----------
    sample_rate_hz : float
        Sampling rate. 2.4e6 is the highest rate the RTL-SDR streams reliably
        on most machines; 2.048e6 is the safe fallback if samples are dropped.

    gain_db : float
        Tuner gain in dB. Must be one of the discrete values the tuner
        supports; the closest valid value is selected and reported.

    freq_correction_ppm : int
        Crystal error correction. Leave at 0 unless it has been calibrated
        against a known transmitter.

    Returns
    -------
    sdr : RtlSdr
        Open device, caller is responsible for closing it.

    settings : dict
        What was actually applied, for the run metadata.
    """

    from rtlsdr import RtlSdr  # imported here so the module loads without hardware

    sdr = RtlSdr()
    sdr.sample_rate = sample_rate_hz

    # RTL2832 digital AGC off. Older pyrtlsdr builds do not expose this.
    try:
        sdr.set_agc_mode(False)
    except Exception:
        pass

    # Manual tuner gain. Assigning a number already selects manual mode in
    # pyrtlsdr, but being explicit costs nothing and documents the intent.
    try:
        sdr.set_manual_gain_enabled(True)
    except Exception:
        pass

    valid = list(getattr(sdr, "valid_gains_db", []) or [])
    applied_gain = min(valid, key=lambda g: abs(g - gain_db)) if valid else gain_db
    sdr.gain = applied_gain

    if freq_correction_ppm:
        sdr.freq_correction = int(freq_correction_ppm)

    settings = {
        "sample_rate_hz": float(sdr.sample_rate),
        "gain_db_requested": float(gain_db),
        "gain_db_applied": float(applied_gain),
        "valid_gains_db": valid,
        "agc": False,
        "freq_correction_ppm": int(freq_correction_ppm),
    }

    return sdr, settings


def capture_at(sdr, center_hz, num_samples, settle_samples=65536):
    """
    Retune and read one block of IQ samples.

    The first block after a retune is discarded. The PLL needs time to lock and
    the tuner's own AGC loops need time to settle; samples read immediately
    after setting center_freq contain a transient that shows up as a smear
    across the whole band.

    Parameters
    ----------
    sdr : RtlSdr
        Device from open_sdr().

    center_hz : float
        Centre frequency in Hz.

    num_samples : int
        Samples to keep. Rounded up to a multiple of 512, which is what
        librtlsdr transfers in.

    settle_samples : int
        Samples to read and throw away after retuning.

    Returns
    -------
    np.ndarray
        Complex IQ samples.
    """

    if not V4_FREQ_MIN_HZ <= center_hz <= V4_FREQ_MAX_HZ:
        raise ValueError(
            f"centre {center_hz/1e6:.1f} MHz is outside the V4 range "
            f"{V4_FREQ_MIN_HZ/1e6:.0f}-{V4_FREQ_MAX_HZ/1e6:.0f} MHz"
        )

    sdr.center_freq = center_hz

    if settle_samples > 0:
        sdr.read_samples(_round512(settle_samples))

    return sdr.read_samples(_round512(num_samples))


def check_clipping(x, threshold=0.99, max_fraction=1e-4):
    """
    Warn if the front end is being overdriven.

    An 8 bit receiver clips ungracefully, and a clipped capture produces
    intermodulation products that look exactly like real signals. Rather than
    trusting the level, count how many samples are at the rail.

    Parameters
    ----------
    x : np.ndarray
        Complex IQ samples.

    threshold : float
        Magnitude counted as full scale, per I and Q component.

    max_fraction : float
        Fraction of samples above which the capture is called clipped.

    Returns
    -------
    clipped : bool

    fraction : float
        Fraction of samples at or above the threshold.
    """

    at_rail = (np.abs(x.real) >= threshold) | (np.abs(x.imag) >= threshold)
    fraction = float(np.mean(at_rail))

    return fraction > max_fraction, fraction


def _round512(n):
    """librtlsdr transfers in multiples of 512 samples."""
    return int(np.ceil(n / 512) * 512)
