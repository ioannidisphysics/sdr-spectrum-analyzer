import numpy as np


def generate_synthetic_iq(fs, duration, random_seed=None):
    """
    Generate synthetic complex IQ samples.

    Signals included:
    - Tone at +100 kHz
    - Tone at -220 kHz
    - AM signal centered at +250 kHz
    - Chirp from -350 kHz to -50 kHz
    - Complex Gaussian noise

    Parameters
    ----------
    fs : float
        Sampling frequency in Hz.

    duration : float
        Signal duration in seconds.

    Returns
    -------
    t : np.ndarray
        Time vector.

    x : np.ndarray
        Complex IQ samples.

    metadata : dict
        Information about the generated signals.
    """

    num_samples = int(fs * duration)
    t = np.arange(num_samples) / fs
    rng = np.random.default_rng(random_seed)
    # -----------------------------
    # Signal definitions
    # -----------------------------

    f_signal_1 = 100_000      # +100 kHz
    f_signal_2 = -220_000     # -220 kHz

    f_carrier = 250_000       # AM carrier at +250 kHz
    f_mod = 5_000             # AM modulation frequency = 5 kHz
    mod_index = 0.8           # AM modulation depth

    f_chirp_start = -350_000
    f_chirp_end = -50_000

    # -----------------------------
    # Generate tones
    # -----------------------------

    signal_1 = 1.0 * np.exp(1j * 2 * np.pi * f_signal_1 * t)
    signal_2 = 0.5 * np.exp(1j * 2 * np.pi * f_signal_2 * t)

    # -----------------------------
    # Generate AM signal
    # -----------------------------

    modulation = 1 + mod_index * np.cos(2 * np.pi * f_mod * t)
    signal_3 = 0.35 * modulation * np.exp(1j * 2 * np.pi * f_carrier * t)

    # -----------------------------
    # Generate chirp signal
    # -----------------------------

    chirp_rate = (f_chirp_end - f_chirp_start) / duration

    chirp_phase = 2 * np.pi * (
        f_chirp_start * t + 0.5 * chirp_rate * t**2
    )

    signal_4 = 0.15 * np.exp(1j * chirp_phase)

    # -----------------------------
    # Complex Gaussian noise
    # -----------------------------

    noise_power = 0.01
    noise = np.sqrt(noise_power / 2) * (
            rng.standard_normal(len(t)) + 1j * rng.standard_normal(len(t))

    )

    # -----------------------------
    # Final IQ signal
    # -----------------------------

    x = signal_1 + signal_2 + signal_3 + signal_4 + noise

    metadata = {
        "f_signal_1": f_signal_1,
        "f_signal_2": f_signal_2,
        "f_carrier": f_carrier,
        "f_mod": f_mod,
        "mod_index": mod_index,
        "f_chirp_start": f_chirp_start,
        "f_chirp_end": f_chirp_end,
        "noise_power": noise_power,
        "random_seed": random_seed,
    }

    return t, x, metadata