"""
Frequency allocation table for ITU Region 1, as it applies in Greece.

A wideband sweep produces a list of frequencies where something was found.
On its own that list says very little. The same list read against the band
plan says which service each detection belongs to, and that is what turns a
trace into a measurement someone can check: an FM carrier has to land on an
odd multiple of 100 kHz inside 87.5 to 108 MHz, a DAB block has to be 1.536
MHz wide and centred on one of the Band III block frequencies, and a GSM
downlink carrier has to sit in 925 to 960 MHz and not in the uplink half.

Where a detection contradicts the table, the detection is the thing that
needs explaining. Most of the time the explanation is in signals.py: a
harmonic of something strong, or a spur from the 28.8 MHz reference.

Entries are (f_lo_hz, f_hi_hz, short_name, description, category). They may
overlap, and the lookup returns the narrowest entry that contains the
frequency, so 1090 MHz reads as ADS-B rather than as the aeronautical band
that contains it.

Sources for the boundaries are the ITU Radio Regulations Region 1 table and
the Greek national frequency allocation plan (EETT). The edges of the mobile
bands follow 3GPP band numbering, which is what the operators actually use.
"""

import bisect

MHZ = 1e6

# Categories drive the colour of the shaded regions in the overview plot.
CATEGORIES = (
    "broadcast",
    "mobile",
    "aeronautical",
    "maritime",
    "amateur",
    "satellite",
    "land_mobile",
    "srd",
    "radar",
    "passive",
)

# (f_lo_hz, f_hi_hz, short_name, description, category)
BAND_PLAN = [
    (26.960, 27.410, "CB", "Citizens band, 40 channels, AM and SSB", "land_mobile"),
    (28.000, 29.700, "10 m", "Amateur 10 m", "amateur"),
    (30.000, 47.000, "VHF low", "Land mobile, low band", "land_mobile"),
    (47.000, 68.000, "Band I", "Former VHF TV Band I, now land mobile", "land_mobile"),
    (68.000, 74.800, "VHF LM", "Land mobile", "land_mobile"),
    (74.800, 75.200, "Marker", "Aeronautical marker beacons, 75 MHz", "aeronautical"),
    (75.200, 87.500, "VHF LM", "Land mobile and aeronautical", "land_mobile"),
    (87.500, 108.000, "FM", "FM broadcast, Band II, 100 kHz raster", "broadcast"),
    (108.000, 117.975, "VOR/ILS", "Aeronautical navigation, VOR and ILS localiser", "aeronautical"),
    (117.975, 137.000, "Airband", "Aeronautical voice, AM, 8.33 kHz channels", "aeronautical"),
    (137.000, 138.000, "Weather sat", "NOAA APT and Meteor-M downlink, Orbcomm", "satellite"),
    (138.000, 144.000, "VHF mil", "Military and land mobile", "land_mobile"),
    (144.000, 146.000, "2 m", "Amateur 2 m", "amateur"),
    (146.000, 156.000, "VHF LM", "Land mobile", "land_mobile"),
    (156.000, 162.050, "Marine", "Maritime VHF, 25 kHz channels", "maritime"),
    (161.950, 162.050, "AIS", "AIS 161.975 and 162.025 MHz", "maritime"),
    (162.050, 174.000, "VHF LM", "Land mobile and paging", "land_mobile"),
    (174.000, 230.000, "DAB", "Band III, DAB+ blocks 5A to 12D, 1.536 MHz wide", "broadcast"),
    (225.000, 328.600, "UHF mil", "Military aeronautical voice, AM", "aeronautical"),
    (242.975, 243.025, "Guard", "Military distress guard frequency", "aeronautical"),
    (328.600, 335.400, "ILS GS", "ILS glideslope", "aeronautical"),
    (335.400, 380.000, "UHF mil", "Military aeronautical", "aeronautical"),
    (380.000, 400.000, "TETRA BOS", "TETRA for emergency services, 25 kHz carriers", "land_mobile"),
    (400.150, 406.000, "Met aids", "Radiosondes and meteorological satellites", "satellite"),
    (406.000, 406.100, "SAR", "COSPAS-SARSAT emergency beacons", "satellite"),
    (406.100, 410.000, "UHF LM", "Land mobile", "land_mobile"),
    (410.000, 430.000, "TETRA", "TETRA and trunked PMR", "land_mobile"),
    (430.000, 440.000, "70 cm", "Amateur 70 cm, includes repeaters at 433 and 438", "amateur"),
    (433.050, 434.790, "ISM 433", "433 MHz ISM, remotes, sensors, weather stations", "srd"),
    (440.000, 446.000, "UHF LM", "Land mobile", "land_mobile"),
    (446.000, 446.200, "PMR446", "Licence free handhelds, 16 channels", "srd"),
    (446.200, 470.000, "UHF LM", "Land mobile and paging", "land_mobile"),
    (470.000, 694.000, "DVB-T", "UHF television, channels 21 to 48, 8 MHz raster", "broadcast"),
    (694.000, 790.000, "5G 700", "3GPP n28, uplink 703-733, downlink 758-788", "mobile"),
    (791.000, 821.000, "LTE800 DL", "3GPP n20 downlink", "mobile"),
    (832.000, 862.000, "LTE800 UL", "3GPP n20 uplink", "mobile"),
    (863.000, 870.000, "SRD 868", "868 MHz short range devices, LoRa, Sigfox, meters", "srd"),
    (870.000, 880.000, "UHF LM", "Land mobile", "land_mobile"),
    (880.000, 915.000, "GSM900 UL", "GSM, UMTS and LTE900 uplink, handsets transmit here", "mobile"),
    (915.000, 925.000, "Guard", "Guard band between the GSM900 halves", "mobile"),
    (925.000, 960.000, "GSM900 DL", "GSM, UMTS and LTE900 downlink, base stations", "mobile"),
    (960.000, 1164.000, "DME/TACAN", "Aeronautical ranging, 1 MHz channels", "aeronautical"),
    (1087.000, 1093.000, "ADS-B", "Aircraft transponders and ADS-B at 1090 MHz", "aeronautical"),
    (1164.000, 1215.000, "GNSS L5", "Galileo E5 and GPS L5 at 1176.45 MHz", "satellite"),
    (1215.000, 1240.000, "GNSS L2", "GPS L2 at 1227.6 MHz, GLONASS L2, radiolocation", "satellite"),
    (1240.000, 1300.000, "23 cm", "Amateur 23 cm, shared with radiolocation", "amateur"),
    (1300.000, 1350.000, "Air radar", "Aeronautical radionavigation radar", "radar"),
    (1350.000, 1400.000, "Radar", "Radiolocation and fixed links", "radar"),
    (1400.000, 1427.000, "Passive", "Radio astronomy and passive sensing, no transmission allowed", "passive"),
    (1427.000, 1452.000, "SRD/LM", "Land mobile and short range devices", "srd"),
    (1452.000, 1492.000, "SDL", "Supplemental downlink, former L band DAB", "mobile"),
    (1492.000, 1518.000, "LM", "Land mobile and fixed", "land_mobile"),
    (1518.000, 1559.000, "MSS DL", "Mobile satellite downlink, Inmarsat", "satellite"),
    (1559.000, 1610.000, "GNSS L1", "GPS L1 1575.42, Galileo E1, GLONASS 1598-1606", "satellite"),
    (1610.000, 1626.500, "Iridium", "Iridium satellite links", "satellite"),
    (1626.500, 1660.500, "MSS UL", "Mobile satellite uplink, Inmarsat", "satellite"),
    (1660.500, 1670.000, "Passive", "Radio astronomy", "passive"),
    (1670.000, 1710.000, "Met sat", "Meteorological satellite downlink, NOAA HRPT", "satellite"),
    (1710.000, 1785.000, "LTE1800 UL", "3GPP n3 uplink, GSM1800 handsets", "mobile"),
    (1805.000, 1880.000, "LTE1800 DL", "3GPP n3 downlink, GSM1800 base stations", "mobile"),
    (1880.000, 1900.000, "DECT", "Cordless telephones, 10 carriers 1.728 MHz apart", "srd"),
    (1900.000, 1920.000, "TDD", "Unpaired TDD spectrum", "mobile"),
    (1920.000, 1980.000, "UMTS2100 UL", "3GPP n1 uplink", "mobile"),
    (2110.000, 2170.000, "UMTS2100 DL", "3GPP n1 downlink", "mobile"),
    (2400.000, 2483.500, "ISM 2.4", "Wi-Fi, Bluetooth, microwave ovens, the patch antenna band", "srd"),
]

# Store in Hz and sort by lower edge.
BAND_PLAN = sorted(
    ((lo * MHZ, hi * MHZ, short, text, cat) for lo, hi, short, text, cat in BAND_PLAN),
    key=lambda row: (row[0], row[1]),
)

_LOWER_EDGES = [row[0] for row in BAND_PLAN]

# Tuning range of the RTL-SDR Blog V4, repeated here so callers can mark the
# part of the band plan the hardware cannot reach without a converter.
RX_MIN_HZ = 24e6
RX_MAX_HZ = 1766e6


def lookup(f_hz):
    """
    Band plan entry containing a frequency.

    Overlapping entries are resolved in favour of the narrowest one, so a
    detection at 1090 MHz is reported as ADS-B rather than as the 960 to 1164
    MHz aeronautical block that also contains it.

    Parameters
    ----------
    f_hz : float
        Frequency in Hz.

    Returns
    -------
    dict or None
        Keys short, description, category, f_lo_hz, f_hi_hz. None if the
        frequency falls outside every entry in the table.
    """

    # Every entry whose lower edge is at or below f_hz is a candidate; walking
    # back from the insertion point is enough because the list is sorted by
    # lower edge and no entry in this table is wider than about 220 MHz.
    index = bisect.bisect_right(_LOWER_EDGES, f_hz)
    best = None

    for lo, hi, short, text, cat in reversed(BAND_PLAN[:index]):
        if hi < f_hz:
            continue
        if best is None or (hi - lo) < (best[1] - best[0]):
            best = (lo, hi, short, text, cat)

    if best is None:
        return None

    lo, hi, short, text, cat = best

    return {
        "short": short,
        "description": text,
        "category": cat,
        "f_lo_hz": lo,
        "f_hi_hz": hi,
    }


def label(f_hz):
    """Short name of the band containing f_hz, or 'unallocated'."""

    entry = lookup(f_hz)

    return entry["short"] if entry else "unallocated"


def bands_in_range(f_start_hz, f_stop_hz, categories=None):
    """
    Every band plan entry overlapping a range.

    Used by the plots to shade the services behind the trace and by the report
    to group detections by band.

    Parameters
    ----------
    f_start_hz, f_stop_hz : float
        Range of interest.

    categories : iterable of str or None
        Restrict to these categories.

    Returns
    -------
    list of dict
        Same keys as lookup(), in ascending frequency order.
    """

    out = []

    for lo, hi, short, text, cat in BAND_PLAN:
        if hi < f_start_hz or lo > f_stop_hz:
            continue
        if categories is not None and cat not in categories:
            continue
        out.append(
            {
                "short": short,
                "description": text,
                "category": cat,
                "f_lo_hz": lo,
                "f_hi_hz": hi,
            }
        )

    return out


def fm_channel(f_hz, tolerance_hz=50e3):
    """
    Nearest FM broadcast channel, as a check on the frequency axis.

    FM carriers in Region 1 sit on odd multiples of 100 kHz, which makes the
    band a free frequency reference: the offset between a detected carrier and
    the nearest channel is the crystal error of the receiver plus whatever the
    sweep stitching contributes. A dongle without a TCXO is typically tens of
    ppm off, which at 100 MHz is a few kHz.

    Parameters
    ----------
    f_hz : float
        Detected carrier frequency.

    tolerance_hz : float
        Largest offset still counted as a match.

    Returns
    -------
    dict or None
        Keys channel_hz and offset_hz, or None if the frequency is outside
        the FM band or further than the tolerance from any channel.
    """

    if not 87.5e6 <= f_hz <= 108e6:
        return None

    channel_hz = round(f_hz / 100e3) * 100e3

    # Region 1 uses the odd multiples of 100 kHz, 87.6, 87.8 and so on.
    if int(round(channel_hz / 100e3)) % 2 == 0:
        lower, upper = channel_hz - 100e3, channel_hz + 100e3
        channel_hz = lower if abs(f_hz - lower) < abs(f_hz - upper) else upper

    offset_hz = f_hz - channel_hz

    if abs(offset_hz) > tolerance_hz:
        return None

    return {"channel_hz": channel_hz, "offset_hz": offset_hz}
