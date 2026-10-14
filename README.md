[![PRESTO on ASCL](https://img.shields.io/badge/ascl-1107.017-blue.svg)](https://ascl.net/1107.017)
[![conda-forge version](https://img.shields.io/conda/vn/conda-forge/presto-pulsar.svg)](https://anaconda.org/conda-forge/presto-pulsar)
[![conda-forge downloads](https://img.shields.io/conda/dn/conda-forge/presto-pulsar.svg)](https://anaconda.org/conda-forge/presto-pulsar)

# [PRESTO](https://github.com/scottransom/presto/)

PRESTO is a large suite of pulsar search and analysis software developed primarily by Scott Ransom mostly from scratch, and released under the GPL (v2). It was primarily designed to efficiently search for binary millisecond pulsars from long observations of globular clusters (although it has since been used in several surveys with short integrations and to process a lot of X-ray data as well). It is written primarily in ANSI C, with many of the recent routines in Python. According to Steve Eikenberry, PRESTO stands for: PulsaR Exploration and Search TOolkit!

**PRESTO has discovered well over 1500 pulsars, including ~500 recycled and/or binary pulsars!**

## Version 6.1.1:
 * Bug fix (#231): `prepfold -accelcand` with a nonzero `-start` now propagates jerk
   (`-wmax`) candidates correctly to the new fold start. Before, `f` and `f-dot` could be
   badly off, often enough to wash out the pulse entirely. Candidates with `w = 0` were not
   affected.
 * Bug fix: signed 8-bit SIGPROC filterbank data is now read correctly on aarch64 (64-bit
   ARM Linux and Apple Silicon). x86 was never affected.

## Version 6.1.0:
 * **Zapped Fourier bins are now filled with noise rather than a constant amplitude**
   (`zapbirds`, `simple_zapbirds.py`, `accelsearch -zaplist`). This changes search results on
   any zapped data; pass `-constamp` for the old behavior. Several other `zapbirds` fixes too.
 * `accelsearch` bug fix (present since 2001): the f-fdot candidate optimizer sized its
   interpolation kernel from a scaled `z`, so it recovered too little power at large `|z|`.
   Happily this did **not** cost anyone candidates &mdash; the sigmas and summed powers
   `accelsearch` reports come from the search itself and are unchanged &mdash; but
   `presto.maximize_rz()`, `explorefft` and `quicklook` did under-report power above
   `|z|` ~ 150. See [CHANGELOG.md](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
   for the details.
 * `single_pulse_search.py` no longer leaves up to 7999 samples at the end of each time series
   unsearched.
 * SIGPROC filterbank bug fix: declinations between 0 and -1 degrees lost their minus sign.
 * Red noise removal is now a shared `libpresto` routine, with new `rednoise -inplace` and
   `realfft -rednoise` options so you can whiten without a separate pass over the file.
 * New `presto.observatories` module: a single home for the telescope-name to observatory-code
   table that used to be copied in three places.
 * `ACCEL_sift.py` and `sifting.DM_re` now find DMs anywhere in a filename, and
   `simple_zapbirds.py` parses its command line with `argparse`.
 * Removed the `presto.singlepulse` subpackage &mdash; it turned out never to have been
   installed &mdash; along with `make_spd.py`, `plot_spd.py` and `rrattrap.py`.
 * `README.md` and `INSTALL.md` now recommend installing from conda-forge.

## Version 6.0.1:
 * `get_TOAs.py` now supports **CHIME** (telescope code `y`, observatory name `chime`).
 * The conda-forge recipe moved to the v1 (rattler-build) `recipe.yaml` format, and now uses
   `astropy-base` and `gsl >=2.8` (the latter so the package is co-installable with tempo2
   while conda-forge's gsl 2.8 migration is in progress).

## Version 6.0.0:
 * **PRESTO is now on conda-forge!** Install the whole suite (C tools + Python package) with
   `conda install -c conda-forge presto-pulsar` (or `pixi add presto-pulsar`).  (The conda
   package is `presto-pulsar` because the plain `presto` name is taken on bioconda by an
   unrelated tool; the Python import is still `presto`.)
 * The `PRESTO` environment variable is **no longer required** &mdash; runtime data (the pulsar
   catalog, etc.) and man pages now install into the prefix and are found automatically.
 * **TEMPO is no longer needed at all:** barycentering is done in-process with ERFA, and
   polyco generation now uses tempo2.
 * The compiled code is now **Fortran-free**: FFTFIT was rewritten in pure Python (NumPy/SciPy,
   with a new `classic`/`aarchiba` selector in `get_TOAs.py`), and the last Fortran least-squares
   solver was replaced by GSL.
 * Vendored Numerical-Recipes/SLALIB/RANDLIB code replaced by GSL and ERFA; many Python modules
   modernized (NumPy docstrings, type hints, NumPy 2.5 fixes); man pages fleshed out.
 * The build uses meson (from v5); `pixi run build` drives a fully self-contained build.

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

1. Data Preparation: Interference detection (`rfifind`) and removal (`zapbirds` and `pfdzap.py`), de-dispersion (`prepdata`, `prepsubband`, and `mpiprepsubband`), barycentering (via TEMPO).
2. Searching: Fourier-domain acceleration and jerk (`accelsearch`), single-pulse (`single_pulse_search.py`), and phase-modulation or sideband searches (`search_bin`).
3. Folding: Candidate optimization (`prepfold` and `fourier_fold.py`) and Time-of-Arrival (TOA) generation (`get_TOAs.py`).
4. Misc: Data exploration (`readfile`, `exploredat`, `explorefft`), de-dispersion planning (`DDplan.py`), date conversion (`mjd2cal`, `cal2mjd`), tons of python pulsar/astro libraries, average pulse creation and flux density estimation (`sum_profiles.py`), and more...
5. Post Single Pulse Searching Tools: Single pulse plotting (`waterfaller.py`).

Many additional utilities are provided for various tasks that are often required when working with pulsar data such as time conversions, Fourier transforms, time series and FFT exploration, byte-swapping, etc.

**References**: The Fourier-Domain acceleration search technique that PRESTO uses in the routine `accelsearch` is described in [Ransom, Eikenberry, and Middleditch (2002)](https://ui.adsabs.harvard.edu/abs/2002AJ....124.1788R/abstract), the "jerk" search capability is described in [Andersen & Ransom (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...863L..13A/abstract), and the phase-modulation search technique used by `search_bin` is described in [Ransom, Cordes, and Eikenberry (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..911R/abstract). Some other basic information about PRESTO can be found in my [thesis](http://www.cv.nrao.edu/~sransom/ransom_thesis_2001.pdf).

**Support/Docs**:  I may eventually get around to finishing the documentation for PRESTO (but probably not), but until then you should know that each routine returns its basic usage when you call it with no arguments. I am also willing to provide limited support via email (see below). And make sure to check out the [FAQ.md](https://github.com/scottransom/presto/blob/master/FAQ.md)!

**Tutorial**: There is a tutorial in the "docs" directory which walks you through all the main steps of finding pulsars using PRESTO.

## Installing it:
**The recommended way to install PRESTO is from [conda-forge](https://conda-forge.org)**, which
gives you the whole suite &mdash; `libpresto`, all of the command-line tools, the man pages, and
the `presto` Python package &mdash; with nothing to compile and no environment variables to set:

    conda install -c conda-forge presto-pulsar

or, if you use [pixi](https://pixi.sh):

    pixi add presto-pulsar

The conda package is named `presto-pulsar` (the plain `presto` name is taken on bioconda by an
unrelated tool), but the Python import is still `import presto` and all of the tool names are
unchanged. Packages are built for linux-64 and osx-64 (Apple Silicon runs the osx-64 build
under Rosetta 2 for now). tempo2 (for `prepfold
-timing` polycos) comes along automatically, and TEMPO is not needed at all.

Build PRESTO from source instead if you want to hack on the code, need `mpiprepsubband` (MPI),
or want a build tuned to your own CPU. See
[INSTALL.md](https://github.com/scottransom/presto/blob/master/INSTALL.md) for the full
instructions; the short version, if you use pixi, is `pixi run build` in a git checkout.

## Getting the source: 
The PRESTO source code is released under the GPL and can be browsed or gotten from here in many different ways (including zipped or tar'd or via git). If you are too lazy to read how to get it but have git on your system do:

    git clone https://github.com/scottransom/presto.git

To update it on a regular basis do

    cd presto
    git pull

and then re-build.

If you want the "classic" branch, do the following:

    git clone https://github.com/scottransom/presto.git
    cd presto
    git checkout -b classic origin/classic

then build as per the (old) INSTALL file.

### Development:
If you plan to tweak the code, I highly suggest that you use git and clone the directory (or fork it using an account on github).  And if you want to contribute your changes back, please give me a "pull request"!

**Code contributions and/or patches to fix bugs are most welcome!**

### Citing PRESTO in your work:
If you use PRESTO at all, please cite it via the [ASCL](https://www.ascl.net/1107.017), which you can do using the [BibTeX entry from ADS](https://scixplorer.org/abs/2011ascl.soft07017R/exportcitation/bibtex). I also request that you add a link to the [PRESTO github repo](https://github.com/scottransom/presto). Some people are still citing my [thesis](https://scixplorer.org/abs/2001PhDT.......123R/abstract), which is nice, but probably not required anymore!

In addition, if you use the specific capabilities for searching that have been published, please cite the relevant papers:
 * For Fourier-domain acceleration searches: [Ransom, Eikenberry, and Middleditch (2002)](https://ui.adsabs.harvard.edu/abs/2002AJ....124.1788R/abstract)
 * For jerk searches: [Andersen & Ransom (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...863L..13A/abstract)
 * For sideband / phase modulation searches: [Ransom, Cordes, and Eikenberry (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..911R/abstract)

### Final Thoughts:
Please let me know if you find pulsars or RRATs or FRBs with PRESTO!

Thanks!

### Acknowledgements:
Big thanks go to Steve Eikenberry for his help developing the algorithms, Dunc Lorimer and David Kaplan for help with (retired) code to process BCPM, SCAMP, and Spigot data, among other things, Jason Hessels and Patrick Lazarus for many contributions to the early Python routines, Bridget Andersen for crucial contributions to jerk search, and (alphabetical): Anne Archibald, Cees Bassa, Matteo Bachetti, Slavko Bogdanov, Fernando Camilo, Shami Chatterjee, Kathryn Crowter, Paul Demorest, Paulo Freire, Nate Garver-Daniels, Chen Karako, Mike Keith, Maggie Livingstone, Ryan Lynch, Erik Madsen, Bradley Meyers, Gijs Molenaar, Timothy Olszanski, Chitrang Patel, Paul Ray, Alessandro Ridolfi, Paul Scholz, Maciej Serylak, Ingrid Stairs, Kevin Stovall, Nick Swainston, Joeri van Leeuwen, Erum Vohra, and Wenke Xia for many comments, suggestions and patches!

Scott Ransom <sransom@nrao.edu>
