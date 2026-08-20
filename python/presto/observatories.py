"""Observatory names and the codes that PRESTO's barycentering uses.

The two-character ITOA codes below are what ``presto.presto.barycenter()``
(and the C code behind it) expects.  They are the same codes that TEMPO's
``obsys.dat`` uses, and the coordinates that go with them live in PRESTO's
``src/observatories.c``.

This table mirrors ``telescope_to_tempocode()`` in PRESTO's
``src/misc_utils.c``, which is not currently wrapped for python.  Keep the
two in sync.  Note that these are *not* the one-character TEMPO ids or the
tempo2 site names that :mod:`presto.polycos` deals with.
"""

from __future__ import annotations

# Telescope name (lowercased, as it appears in a .inf file) to the
# two-character ITOA observatory code, with the familiar name that goes
# with it.
# fmt: off
telescope_to_itoa = {
    "gbt":        ("GB", "GBT"),
    "arecibo":    ("AO", "Arecibo"),
    "vla":        ("VL", "VLA"),
    "parkes":     ("PK", "Parkes"),
    "jodrell":    ("JB", "Jodrell Bank"),
    "gb43m":      ("G1", "GB43m"),
    "gb 140ft":   ("G1", "GB43m"),
    "nrao20":     ("G1", "GB43m"),
    "nancay":     ("NC", "Nancay"),
    "effelsberg": ("EF", "Effelsberg"),
    "srt":        ("SR", "Sardinia Radio Telescope"),
    "fast":       ("FA", "FAST"),
    "wsrt":       ("WT", "WSRT"),
    "gmrt":       ("GM", "GMRT"),
    "chime":      ("CH", "CHIME"),
    "lofar":      ("LF", "LOFAR"),
    "lwa":        ("LW", "LWA1"),
    "mwa":        ("MW", "MWA"),
    "meerkat":    ("MK", "MeerKAT"),
    "ata":        ("AT", "ATA"),
    "k7":         ("K7", "KAT-7"),
    "geocenter":  ("0 ", "Geocenter"),
}
# fmt: on


def telescope_to_tempocode(telescope: str) -> tuple[str, str]:
    """Return the ITOA observatory code and familiar name of a telescope.

    Parameters
    ----------
    telescope : str
        The telescope name, as it appears in a ``.inf`` file or a raw data
        header (e.g. "Parkes", "GBT", "Effelsberg").  The comparison is
        case-insensitive and ignores surrounding whitespace.

    Returns
    -------
    obscode : str
        The two-character ITOA code that ``barycenter()`` expects.
    name : str
        The familiar name of the observatory.

    Raises
    ------
    KeyError
        If the telescope is not one that PRESTO knows about.

    Examples
    --------
    >>> telescope_to_tempocode("Parkes")
    ('PK', 'Parkes')
    """
    try:
        return telescope_to_itoa[telescope.strip().lower()]
    except KeyError:
        raise KeyError(
            f"'{telescope}' is not an observatory that PRESTO knows about.  "
            f"The known names are: {', '.join(sorted(telescope_to_itoa))}"
        ) from None


def obscode(telescope: str) -> str:
    """Return just the two-character ITOA code of a telescope.

    Parameters
    ----------
    telescope : str
        The telescope name (see :func:`telescope_to_tempocode`).

    Returns
    -------
    str
        The two-character ITOA code that ``barycenter()`` expects.

    Raises
    ------
    KeyError
        If the telescope is not one that PRESTO knows about.
    """
    return telescope_to_tempocode(telescope)[0]
