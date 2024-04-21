[![PINT on ASCL](https://img.shields.io/badge/ascl-1107.017-blue.svg)](https://ascl.net/1107.017)

# [PRESTO](https://github.com/scottransom/presto/)

PRESTO is a large suite of pulsar search and analysis software developed primarily by Scott Ransom mostly from scratch, and released under the GPL (v2). It was primarily designed to efficiently search for binary millisecond pulsars from long observations of globular clusters (although it has since been used in several surveys with short integrations and to process a lot of X-ray data as well). It is written primarily in ANSI C, with many of the recent routines in Python. According to Steve Eikenberry, PRESTO stands for: PulsaR Exploration and Search TOolkit!

**PRESTO has discovered well over 1000 pulsars, including ~400 recycled and/or binary pulsars!**

## New in Version 5.0.0!:
 * This is a major release since I've moved to a completely different and modern build system: [meson](https://mesonbuild.com/), along with the [meson-python](https://meson-python.readthedocs.io/en/latest/) backend. This was required since *Numpy* has deprecated `numpy.distutils` and this caused python builds to stop working with Python v3.12.
   * See the [INSTALL.md](https://github.com/scottransom/presto/blob/master/INSTALL.md) for updated installation instructions.
   * You will need to install **meson**, **meson-python**, and **ninja**, but that is easily done via `pip`!
   * Python v3.8 or newer is now required.
 * All of the old Spigot-related codes have been removed. If you need to process Spigot data, please use the `classis` branch that is mentioned in the README.md.
 * All of the `slalib` codes (and the python interface to it) have been removed. If you need that stuff, you should transition to [ERFA](https://github.com/liberfa/erfa) and/or [Astropy](https://www.astropy.org/).
 * There are two nice new python utilities:
   * `binary_utils.py` reads a parfile of a binary pulsar and computes min/max observed barycentric spin periods or velocities as either a function of the orbit (default), or for a prescribed duration of time, and optionally plots those. It also shows basic information about the binary.
   * `compare_periods.py` compares candidate spin periods and their integer and fractional harmonics with one or more parfiles to a prescribed fractional tolerance. 

For information on older versions, please see the [CHANGELOG.md](https://github.com/scottransom/presto/blob/master/CHANGELOG.md).

## About PRESTO:
PRESTO is written with portability, ease-of-use, and memory efficiency in mind, it can currently handle raw data from the following pulsar machines or formats:

 * PSRFITS search-format data (as from GUPPI at the GBT, PUPPI and the Mock Spectrometers at Arecibo, and much new and archived data from Parkes)
 * 1-, 2-, 4-, 8-, and 32-bit (float) filterbank format from SIGPROC
 * A time series composed of single precision (i.e. 4-byte) floating point data (with a text ".inf" file describing it)
 * Photon arrival times (or events) in ASCII or double-precision binary formats

Notice that the following formats which *used* to be supported are not:

 * Wideband Arecibo Pulsar Processor (WAPP) at Arecibo
 * The Parkes and Jodrell Bank 1-bit filterbank formats
 * SPIGOT at the GBT
 * Berkeley-Caltech Pulsar Machine (BCPM) at the GBT

If you need to process them, you can either checkout the "classic" branch of PRESTO (see below), which is not being actively developed. Or you can use DSPSR to convert those formats into SIGPROC filterbank or (even better) PSRFITS search format. You can grab DSPSR [here](http://dspsr.sourceforge.net).  If you *really* need to get one of these machines working in modern PRESTO, let me know and we can probably make it happen.

The software is composed of numerous routines designed to handle three main areas of pulsar analysis:

1. Data Preparation: Interference detection (`rfifind`) and removal (`zapbirds`), de-dispersion (`prepdata`, `prepsubband`, and `mpiprepsubband`), barycentering (via TEMPO).
2. Searching: Fourier-domain acceleration and jerk (`accelsearch`), single-pulse (`single_pulse_search.py`), and phase-modulation or sideband searches (`search_bin`).
3. Folding: Candidate optimization (`prepfold`) and Time-of-Arrival (TOA) generation (`get_TOAs.py`).
4. Misc: Data exploration (`readfile`, `exploredat`, `explorefft`), de-dispersion planning (`DDplan.py`), date conversion (`mjd2cal`, `cal2mjd`), tons of python pulsar/astro libraries, average pulse creation, flux density estimation, and more...
5. Post Single Pulse Searching Tools: Grouping algorithm (`rrattrap.py`), Production and of single pulse diagnostic plots (`make_spd.py`, `plot_spd.py`, and `waterfaller.py`).

Many additional utilities are provided for various tasks that are often required when working with pulsar data such as time conversions, Fourier transforms, time series and FFT exploration, byte-swapping, etc.

**References**: The Fourier-Domain acceleration search technique that PRESTO uses in the routine `accelsearch` is described in [Ransom, Eikenberry, and Middleditch (2002)](https://ui.adsabs.harvard.edu/abs/2002AJ....124.1788R/abstract), the "jerk" search capability is described in [Andersen & Ransom (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...863L..13A/abstract), and the phase-modulation search technique used by `search_bin` is described in [Ransom, Cordes, and Eikenberry (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..911R/abstract). Some other basic information about PRESTO can be found in my [thesis](http://www.cv.nrao.edu/~sransom/ransom_thesis_2001.pdf).

**Support/Docs**:  I may eventually get around to finishing the documentation for PRESTO (or not), but until then you should know that each routine returns its basic usage when you call it with no arguments. I am also willing to provide limited support via email (see below). And make sure to check out the `FAQ.md`!

**Tutorial**: There is a tutorial in the "docs" directory which walks you through all the main steps of finding pulsars using PRESTO.

## Getting it: 
The PRESTO source code is released under the GPL and can be browsed or gotten from here in many different ways (including zipped or tar'd or via git). If you are too lazy to read how to get it but have git on your system do:

    git clone git://github.com/scottransom/presto.git

To update it on a regular basis do

    cd $PRESTO
    git pull

and then re-build things in $PRESTO.

For more detailed installation instructions, see [INSTALL.md](https://github.com/scottransom/presto/blob/master/INSTALL.md).

If you want the "classic" branch, do the following:

    git clone git://github.com/scottransom/presto.git
    cd presto
    git checkout -b classic origin/classic

then build as per the (old) INSTALL file.

### Development:
If you plan to tweak the code, I highly suggest that you use git and clone the directory (or fork it using an account on github).  And if you want to contribute your changes back, please give me a "pull request"!

**Code contributions and/or patches to fix bugs are most welcome!**

### Final Thoughts:
Please let me know if you decide to use PRESTO for any "real" searches, especially if you find pulsars with it!

And if you find anything with it, it would be great if you would cite either my thesis or whichever of the three papers listed above is appropriate.

Also note that many people are now also citing software using the ASCL, in addition to the relevant papers: [PRESTO is there!](https://www.ascl.net/1107.017).

Thanks!

### Acknowledgements:
Big thanks go to Steve Eikenberry for his help developing the algorithms, Dunc Lorimer and David Kaplan for help with (retired) code to process BCPM, SCAMP, and Spigot data, among other things, Jason Hessels and Patrick Lazarus for many contributions to the Python routines, and (alphabetical): Bridget Andersen, Anne Archibald, Cees Bassa, Matteo Bachetti, Slavko Bogdanov, Fernando Camilo, Shami Chatterjee, Kathryn Crowter, Paul Demorest, Paulo Freire, Nate Garver-Daniels, Chen Karako, Mike Keith, Maggie Livingstone, Ryan Lynch, Erik Madsen, Bradley Meyers, Gijs Molenaar, Timothy Olszanski, Chitrang Patel, Paul Ray, Alessandro Ridolfi, Paul Scholz, Maciej Serylak, Ingrid Stairs, Kevin Stovall, Nick Swainston, and Joeri van Leeuwen for many comments, suggestions and patches!

Scott Ransom <sransom@nrao.edu>
