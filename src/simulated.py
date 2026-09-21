"""
A simulated RTL-SDR, for developing and testing the survey without hardware.

Two reasons this exists rather than being a convenience.

The whole chain can be tested against a known answer. Every emitter below is
defined by its power spectral density in dBFS per Hz, and the samples are
built by giving each FFT bin the amplitude that density implies. Whatever
welch_psd() reports afterwards can be compared against the number that went
in, which turns the window normalisation, the Welch averaging, the sweep
stitching and the survey rebinning into something testable instead of
something believed. selftest() does exactly that.

The detector, the band identification and the artefact flagging can be
checked against ground truth. The scene contains a harmonic of a strong FM
carrier and a reference oscillator spur, both of which the analysis should
name, and several bursty emitters that a single pass will miss and a
multi pass survey should find.

It is not a propagation model and not a substitute for a measurement. Nothing
produced with it belongs in a results table. Every file written from this
source carries source = simulated in its metadata, and every plot made from
one is titled SIMULATED.
"""

import numpy as np

from .bands import RX_MAX_HZ, RX_MIN_HZ

MHZ = 1e6

# (f_center_hz, bandwidth_hz, psd_dbfs_per_hz, duty_cycle, kind, description)
#
# kind is the shape of the emission, which matters because the analysis reads
# a different frequency from each:
#
#   flat    power spread evenly across the bandwidth, which is what a digital
#           multiplex looks like. It has no peak of its own, so the largest
#           bin lands anywhere inside it and only the centroid is stable.
#   fm      flat sidebands with the carrier standing above them, broadcast FM
#   am      a carrier with narrow sidebands, aeronautical voice
#
# Levels are loosely in the range an RTL-SDR with a small antenna in a city
# reports, but they are chosen to exercise the analysis, not to imitate any
# particular location.
SCENE = [
    # Broadcast FM, always on, 180 kHz wide
    (88.7 * MHZ, 180e3, -78.0, 1.0, "fm", "FM station"),
    (91.6 * MHZ, 180e3, -84.0, 1.0, "fm", "FM station"),
    (94.7 * MHZ, 180e3, -72.0, 1.0, "fm", "FM station, strong"),
    (98.4 * MHZ, 180e3, -88.0, 1.0, "fm", "FM station, weak"),
    (100.3 * MHZ, 180e3, -70.0, 1.0, "fm", "FM station, strongest"),
    (103.7 * MHZ, 180e3, -80.0, 1.0, "fm", "FM station"),
    (107.1 * MHZ, 180e3, -86.0, 1.0, "fm", "FM station"),

    # Second harmonic of the strongest FM carrier, produced in the receiver
    (200.6 * MHZ, 180e3, -97.0, 1.0, "fm", "harmonic 2 of 100.3 MHz"),

    # Reference oscillator leakage, one bin wide
    (5 * 28.8 * MHZ, 2e3, -80.0, 1.0, "flat", "reference spur, 5 x 28.8 MHz"),
    (12 * 28.8 * MHZ, 2e3, -84.0, 1.0, "flat", "reference spur, 12 x 28.8 MHz"),

    # Airband, AM carrier with narrow sidebands, keyed occasionally
    (124.15 * MHZ, 8e3, -92.0, 0.30, "am", "airband voice"),

    # Weather satellite, only when a pass is overhead
    (137.62 * MHZ, 34e3, -92.0, 0.10, "flat", "NOAA APT"),

    # Marine and AIS
    (156.8 * MHZ, 16e3, -94.0, 0.08, "am", "marine VHF channel 16"),
    (161.975 * MHZ, 25e3, -90.0, 0.55, "flat", "AIS"),

    # DAB blocks, 1.536 MHz of flat noise-like power
    (178.352 * MHZ, 1.536e6, -96.0, 1.0, "flat", "DAB block 5C"),
    (192.352 * MHZ, 1.536e6, -98.0, 1.0, "flat", "DAB block 9A"),
    (227.360 * MHZ, 1.536e6, -99.0, 1.0, "flat", "DAB block 12D"),

    # TETRA, bursty
    (390.2 * MHZ, 25e3, -88.0, 0.25, "flat", "TETRA emergency services"),
    (412.6 * MHZ, 25e3, -92.0, 0.12, "flat", "TETRA trunked PMR"),

    # Licence free
    (433.92 * MHZ, 60e3, -90.0, 0.06, "flat", "433 MHz remote control"),
    (446.09 * MHZ, 12.5e3, -93.0, 0.04, "flat", "PMR446 handheld"),

    # DVB-T multiplexes, 7.6 MHz wide
    (490.0 * MHZ, 7.6e6, -99.0, 1.0, "flat", "DVB-T channel 23"),
    (530.0 * MHZ, 7.6e6, -97.0, 1.0, "flat", "DVB-T channel 28"),
    (626.0 * MHZ, 7.6e6, -100.0, 1.0, "flat", "DVB-T channel 40"),

    # Cellular and short range
    (806.0 * MHZ, 10e6, -101.0, 1.0, "flat", "LTE800 downlink"),
    (868.3 * MHZ, 100e3, -95.0, 0.05, "flat", "868 MHz meter"),
    (900.4 * MHZ, 200e3, -93.0, 0.20, "flat", "GSM900 handset uplink"),
    (942.6 * MHZ, 200e3, -86.0, 1.0, "flat", "GSM900 base station"),
    (947.2 * MHZ, 200e3, -89.0, 1.0, "flat", "GSM900 base station"),
    (1090.0 * MHZ, 2e6, -96.0, 0.15, "flat", "ADS-B transponder"),
    (1575.42 * MHZ, 2e6, -110.0, 1.0, "flat", "GPS L1, below the noise floor"),
    (1842.0 * MHZ, 15e6, -103.0, 1.0, "flat", "LTE1800 downlink"),
]

# How far the carrier of a shaped emission stands above its own sidebands, and
# how wide that carrier is.
CARRIER_EXCESS_DB = {"fm": 14.0, "am": 20.0}
CARRIER_WIDTH_HZ = {"fm": 4e3, "am": 2e3}


def emission_psd(f_abs, f0, bandwidth, level_db, kind):
    """
    Power spectral density added by one emitter, linear, on a frequency axis.

    Parameters
    ----------
    f_abs : np.ndarray
        Absolute frequency of each FFT bin.

    f0, bandwidth, level_db, kind :
        One row of SCENE.

    Returns
    -------
    np.ndarray
        Density to add, zero outside the emission.
    """

    psd = np.zeros_like(f_abs)

    inside = np.abs(f_abs - f0) <= bandwidth / 2
    psd[inside] += np.power(10.0, level_db / 10.0)

    excess_db = CARRIER_EXCESS_DB.get(kind)
    if excess_db is not None:
        carrier = np.abs(f_abs - f0) <= CARRIER_WIDTH_HZ[kind] / 2
        psd[carrier] += np.power(10.0, (level_db + excess_db) / 10.0)

    return psd


def noise_floor_dbfs_per_hz(f_hz):
    """
    Receiver noise floor against frequency.

    Shaped to resemble what the hardware does: man made noise lifts the low
    VHF end, the tuner's own noise figure rises towards the top of its range,
    and there is a gentle slope in between. The exact shape does not matter,
    only that it is not flat, because a flat floor would let a single fixed
    threshold work and would hide the reason local_noise_floor() exists.
    """

    f = np.asarray(f_hz, dtype=float)
    low_vhf_excess = 10.0 * np.exp(-(f - RX_MIN_HZ) / 120e6)
    tuner_rise = 6.0 * (f / RX_MAX_HZ) ** 2

    return -114.0 + low_vhf_excess + tuner_rise


class SimulatedSdr:
    """
    Stands in for RtlSdr, with the attributes capture_at() and open_sdr() use.

    Parameters
    ----------
    sample_rate_hz : float
        Sampling rate.

    gain_db : float
        Recorded but not modelled. Gain changes the noise figure and the
        overload point of the real receiver, neither of which this source
        pretends to reproduce.

    seed : int or None
        Fixes the bursts and the noise, so a survey can be repeated exactly.
    """

    valid_gains_db = [
        0.0, 0.9, 1.4, 2.7, 3.7, 7.7, 8.7, 12.5, 14.4, 15.7, 16.6, 19.7,
        20.7, 22.9, 25.4, 28.0, 29.7, 32.8, 33.8, 36.4, 37.2, 38.6, 40.2,
        42.1, 43.4, 43.9, 44.5, 48.0, 49.6,
    ]

    def __init__(self, sample_rate_hz=2.4e6, gain_db=30.0, seed=None):
        self.sample_rate = float(sample_rate_hz)
        self.center_freq = 100e6
        self.freq_correction = 0
        self._gain = min(self.valid_gains_db, key=lambda g: abs(g - gain_db))
        self._rng = np.random.default_rng(seed)
        self.closed = False

    # -- interface expected by sdr_source.capture_at ------------------------

    @property
    def gain(self):
        return self._gain

    @gain.setter
    def gain(self, value):
        self._gain = min(self.valid_gains_db, key=lambda g: abs(g - value))

    def set_agc_mode(self, enabled):
        pass

    def set_manual_gain_enabled(self, enabled):
        pass

    def close(self):
        self.closed = True

    def read_samples(self, num_samples):
        """
        Complex IQ for the current centre frequency.

        The samples are built in the frequency domain. Each bin is given the
        amplitude that the designed power spectral density implies, with a
        random complex phase, and the block is transformed back. Doing it this
        way rather than summing modulated carriers means the density that
        comes out of the analyser is the density that went in, which is what
        makes the source useful as a test.
        """

        n = int(num_samples)
        fs = self.sample_rate
        fc = self.center_freq

        # Unshifted FFT frequency axis, matching np.fft.ifft's bin order.
        f_rel = np.fft.fftfreq(n, d=1.0 / fs)
        f_abs = fc + f_rel

        psd = np.power(10.0, noise_floor_dbfs_per_hz(f_abs) / 10.0)

        for f0, bandwidth, level_db, duty, kind, _ in SCENE:
            if abs(f0 - fc) > fs / 2 + bandwidth:
                continue
            if duty < 1.0 and self._rng.random() > duty:
                continue

            psd += emission_psd(f_abs, f0, bandwidth, level_db, kind)

        # E|X_k|^2 = N * fs * S(f_k) reproduces mean|x|^2 = sum S_k df under
        # numpy's FFT convention, so the analyser reads back S directly.
        amplitude = np.sqrt(n * fs * psd)
        noise = (self._rng.standard_normal(n) + 1j * self._rng.standard_normal(n)) / np.sqrt(2.0)

        return np.fft.ifft(amplitude * noise)


def open_simulated(sample_rate_hz, gain_db, freq_correction_ppm=0, seed=None):
    """
    Mirror of sdr_source.open_sdr() so the callers do not branch.

    Returns
    -------
    sdr : SimulatedSdr

    settings : dict
        Same keys as the hardware path, plus source = simulated.
    """

    sdr = SimulatedSdr(sample_rate_hz, gain_db, seed=seed)

    settings = {
        "sample_rate_hz": float(sdr.sample_rate),
        "gain_db_requested": float(gain_db),
        "gain_db_applied": float(sdr.gain),
        "valid_gains_db": list(SimulatedSdr.valid_gains_db),
        "agc": False,
        "freq_correction_ppm": int(freq_correction_ppm),
        "source": "simulated",
        "seed": seed,
    }

    return sdr, settings


def selftest(sample_rate_hz=2.4e6, nfft=4096, num_samples=262144, seed=1):
    """
    Check the analyser against a source whose density is known exactly.

    Measures the noise floor away from any emitter, and the density inside the
    strongest FM carrier, and compares both against the values in SCENE. An
    error of more than a few tenths of a dB means the window normalisation or
    the Welch scaling is wrong.

    Returns
    -------
    list of dict
        One row per check, with the expected and measured density and the
        error in dB.
    """

    from .psd import welch_psd

    sdr = SimulatedSdr(sample_rate_hz, 30.0, seed=seed)
    checks = []

    # -- empty stretch of spectrum, floor only ------------------------------
    fc = 700e6
    sdr.center_freq = fc
    freqs, psd_db, _ = welch_psd(sdr.read_samples(num_samples), sample_rate_hz, nfft=nfft)

    window = np.abs(freqs) < 0.3 * sample_rate_hz
    measured = float(np.mean(np.power(10.0, psd_db[window] / 10.0)))
    expected = float(np.mean(np.power(10.0, noise_floor_dbfs_per_hz(fc + freqs[window]) / 10.0)))

    checks.append(
        {
            "what": "noise floor at 700 MHz",
            "expected_db": 10 * np.log10(expected),
            "measured_db": 10 * np.log10(measured),
            "error_db": 10 * np.log10(measured / expected),
        }
    )

    # -- sidebands of the strongest carrier ---------------------------------
    f0, bandwidth, level_db, _, kind, _ = next(
        row for row in SCENE if row[5] == "FM station, strongest"
    )
    sdr.center_freq = f0 - 400e3
    freqs, psd_db, _ = welch_psd(sdr.read_samples(num_samples), sample_rate_hz, nfft=nfft)

    f_abs = sdr.center_freq + freqs

    # Away from the carrier itself, so the check is on the sideband density.
    inside = (np.abs(f_abs - f0) < bandwidth / 2 * 0.9) & (np.abs(f_abs - f0) > 20e3)
    measured = float(np.mean(np.power(10.0, psd_db[inside] / 10.0)))
    expected = float(
        np.power(10.0, level_db / 10.0)
        + np.mean(np.power(10.0, noise_floor_dbfs_per_hz(f_abs[inside]) / 10.0))
    )

    checks.append(
        {
            "what": f"sidebands at {f0/1e6:.1f} MHz",
            "expected_db": 10 * np.log10(expected),
            "measured_db": 10 * np.log10(measured),
            "error_db": 10 * np.log10(measured / expected),
        }
    )

    return checks
