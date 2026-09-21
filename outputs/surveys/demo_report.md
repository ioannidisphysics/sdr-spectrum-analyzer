# Wideband survey: demo

> **Simulated source.** These numbers come from `src/simulated.py`, not from
> the radio. They are here to exercise the analysis, and nothing in this
> report describes the actual radio environment.

Generated 2026-09-21 11:13.

## Run

| Setting | Value |
| --- | --- |
| Range | 24.0 to 1766.0 MHz |
| Passes | 6 |
| Tuning steps per pass | 909 |
| Sample rate | 2.400 MS/s |
| Samples per step | 131072 |
| Welch FFT size | 4096, hann window, 50% overlap |
| Resolution bandwidth | 879 Hz |
| Survey grid | 10.0 kHz, 174200 bins |
| Tuner gain | 29.7 dB, AGC off |
| Occupancy threshold | floor + 6 dB |
| Detection threshold | floor + 8 dB on the max hold |
| Noise floor estimator | 10th percentile over 40 MHz |
| Duration | 2.1 min |

## Result

- Signals detected: **26**
- Spectrum above the detection threshold: **2.39%** of the range
- Bins active in at least one pass: **2.39%**
- Bins active in every pass: **1.98%**
- Strongest signal: **100.3005 MHz**, -54.4 dBFS/Hz, 53.7 dB above the floor (FM)
- Detections flagged as probable receiver artefacts: **3**

## By band

| Band | Service | Signals | Strongest (dBFS/Hz) | Best SNR (dB) | Max occupancy |
| --- | --- | --- | --- | --- | --- |
| FM | FM broadcast, Band II, 100 kHz raster | 7 | -54.4 | 53.7 | 100% |
| Airband | Aeronautical voice, AM, 8.33 kHz channels | 1 | -71.4 | 37.5 | 17% |
| Weather sat | NOAA APT and Meteor-M downlink, Orbcomm | 1 | -90.5 | 18.8 | 17% |
| 2 m | Amateur 2 m | 1 | -80.0 | 29.5 | 83% |
| AIS | AIS 161.975 and 162.025 MHz | 1 | -88.9 | 21.0 | 100% |
| DAB | Band III, DAB+ blocks 5A to 12D, 1.536 MHz wide | 4 | -81.5 | 29.2 | 100% |
| UHF mil | Military aeronautical | 1 | -83.5 | 28.6 | 83% |
| TETRA BOS | TETRA for emergency services, 25 kHz carriers | 1 | -86.6 | 25.6 | 67% |
| TETRA | TETRA and trunked PMR | 1 | -91.2 | 21.1 | 17% |
| DVB-T | UHF television, channels 21 to 48, 8 MHz raster | 3 | -94.9 | 17.4 | 100% |
| LTE800 DL | 3GPP n20 downlink | 1 | -98.6 | 13.0 | 100% |
| GSM900 UL | GSM, UMTS and LTE900 uplink, handsets transmit here | 1 | -91.3 | 20.1 | 17% |
| GSM900 DL | GSM, UMTS and LTE900 downlink, base stations | 2 | -84.3 | 26.9 | 100% |
| ADS-B | Aircraft transponders and ADS-B at 1090 MHz | 1 | -94.2 | 16.5 | 17% |

## Strongest 26 signals

| Frequency (MHz) | Width (kHz) | OBW 99% (kHz) | Peak (dBFS/Hz) | SNR (dB) | Occupancy | Band | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 100.3005 | 200.0 | 178.8 | -54.4 | 53.7 | 100% | FM | +0.5 kHz from the 100.3000 channel |
| 94.7001 | 200.0 | 178.1 | -56.9 | 51.0 | 100% | FM | +0.2 kHz from the 94.7000 channel |
| 88.7001 | 200.0 | 180.7 | -62.8 | 44.8 | 100% | FM | -0.0 kHz from the 88.7000 channel |
| 103.7001 | 200.0 | 179.5 | -64.9 | 43.3 | 100% | FM | +0.1 kHz from the 103.7000 channel |
| 91.6000 | 200.0 | 180.8 | -69.1 | 38.7 | 100% | FM |  |
| 107.1000 | 200.0 | 179.3 | -70.8 | 37.6 | 100% | FM | +0.0 kHz from the 107.1000 channel |
| 124.1550 | 20.0 | 9.9 | -71.4 | 37.5 | 17% | Airband |  |
| 98.3998 | 200.0 | 179.5 | -73.0 | 35.1 | 100% | FM |  |
| 144.0050 | 20.0 | 9.9 | -80.0 | 29.5 | 83% | 2 m | possible reference spur, 5 x 28.8 MHz |
| 200.6003 | 200.0 | 178.4 | -81.5 | 29.2 | 100% | DAB | possible harmonic 2 of 100.300 MHz |
| 345.5950 | 20.0 | 9.9 | -83.5 | 28.6 | 83% | UHF mil | possible reference spur, 12 x 28.8 MHz |
| 942.6573 | 220.0 | 207.7 | -84.3 | 26.9 | 100% | GSM900 DL |  |
| 390.1957 | 40.0 | 29.8 | -86.6 | 25.6 | 67% | TETRA BOS |  |
| 947.1529 | 220.0 | 208.0 | -87.4 | 23.8 | 100% | GSM900 DL |  |
| 161.9650 | 30.0 | 19.9 | -88.9 | 21.0 | 100% | AIS |  |
| 137.6170 | 40.0 | 29.8 | -90.5 | 18.8 | 17% | Weather sat |  |
| 412.6150 | 40.0 | 29.8 | -91.2 | 21.1 | 17% | TETRA |  |
| 900.3332 | 190.0 | 179.0 | -91.3 | 20.1 | 17% | GSM900 UL |  |
| 178.0237 | 1550.0 | 1525.6 | -93.8 | 16.5 | 100% | DAB |  |
| 1090.3469 | 1720.0 | 1701.0 | -94.2 | 16.5 | 17% | ADS-B |  |
| 526.8251 | 7620.0 | 7531.7 | -94.9 | 17.4 | 100% | DVB-T |  |
| 192.2953 | 1550.0 | 1524.8 | -95.8 | 14.8 | 100% | DAB |  |
| 487.6152 | 7620.0 | 7533.5 | -96.7 | 15.6 | 100% | DVB-T |  |
| 227.5940 | 1540.0 | 1521.8 | -96.9 | 14.2 | 100% | DAB |  |
| 629.2818 | 7620.0 | 7532.8 | -97.7 | 14.5 | 100% | DVB-T |  |
| 809.2837 | 10020.0 | 9904.5 | -98.6 | 13.0 | 100% | LTE800 DL |  |

## Probable receiver artefacts

These are hypotheses, not conclusions. A reference spur stays when the
antenna is disconnected and a real signal does not. A harmonic falls by
about 20 to 30 dB when the gain is reduced by 10 dB, a real signal by 10.

| Frequency (MHz) | Peak (dBFS/Hz) | Explanation |
| --- | --- | --- |
| 144.0050 | -80.0 | possible reference spur, 5 x 28.8 MHz |
| 200.6003 | -81.5 | possible harmonic 2 of 100.300 MHz |
| 345.5950 | -83.5 | possible reference spur, 12 x 28.8 MHz |

## Frequency axis check

FM broadcast carriers in Region 1 sit on odd multiples of 100 kHz, which
makes the band a free frequency reference. The offset between the detected
carriers and their channels is the crystal error of the receiver plus what
the sweep stitching contributes.

- Carriers matched: **5**
- Median offset: **+0.13 kHz**, which at 100 MHz is **+1.3 ppm**
- Spread: 0.20 kHz standard deviation

Passing `--ppm +1` to the capture removes most of this.

## Figures

![outputs/figures/demo_overview.png](outputs/figures/demo_overview.png)

![outputs/figures/demo_waterfall.png](outputs/figures/demo_waterfall.png)

![outputs/figures/demo_band_fm.png](outputs/figures/demo_band_fm.png)

![outputs/figures/demo_band_airband.png](outputs/figures/demo_band_airband.png)

![outputs/figures/demo_band_2_m.png](outputs/figures/demo_band_2_m.png)

![outputs/figures/demo_band_dab.png](outputs/figures/demo_band_dab.png)

![outputs/figures/demo_band_uhf_mil.png](outputs/figures/demo_band_uhf_mil.png)

![outputs/figures/demo_band_gsm900_dl.png](outputs/figures/demo_band_gsm900_dl.png)

## Limits

- Levels are dBFS per Hz. The receiver has no absolute calibration, so they
  compare within a survey and against another survey at the same gain, and
  are not dBm. Converting them would need a source of known power.
- The trace is the product of the antenna, the cable and the receiver as much
  as of what is transmitting. A band that looks empty may be a band the
  antenna cannot hear.
- Occupancy is sampled, not monitored. Each bin is observed for a few tens of
  milliseconds per pass, so a channel that transmits rarely can be missed
  entirely: with 6 passes the survey saw each bin 6 times.
- An eight bit receiver overloaded by a strong transmitter invents harmonics.
  The flags above are the first check, reducing the gain is the second.
- The tuner reaches 1766 MHz. Wi-Fi, Bluetooth and the
  2.45 GHz patch antenna are above that and need a downconverter.
