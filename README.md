# PRESTO

http://www.cv.nrao.edu/~sransom/presto/

PRESTO is a large suite of pulsar search and analysis software
developed by Scott Ransom mostly from scratch.  It was primarily
designed to efficiently search for binary millisecond pulsars from
long observations of globular clusters (although it has since been
used in several surveys with short integrations and to process a lot
of X-ray data as well).  It is written primarily in ANSI C, with many
of the recent routines in Python.  According to Steve Eikenberry,
PRESTO stands for: **PulsaR Exploration and Search TOolkit**!

**To date, PRESTO has discovered over 300 pulsars, including
more than 150 recycled pulsars, most of which are in binaries!**

## New in Version 2:
 * WAPP, BCPM, Spigot, and 1-bit analog filterbank data are deprecated! 
   (see below)
 * Dramatically improved internal handling (giving better dynamic
   range and RFI removal) of PSRFITS and SIGPROC filterbank data
 * Massive speed-ups (factors of 2 or more) of `accelsearch` when
   all of the F-Fdot plane can fit into core memory (that can be set
   by changing values in `include/meminfo.h`)
 * Many bug fixes and several new scripts (including new orbit fitters)


## About PRESTO:
PRESTO is written with portability, ease-of-use, and memory efficiency
in mind, it can currently handle raw data from the following pulsar
machines or formats:

 * PSRFITS search-format data (as from GUPPI at the GBT, PUPPI and
   the Mock Spectrometers at Arecibo, and much new and archived data
   from Parkes)
 * 1-, 2-, 4-, 8-, and 32-bit (float) filterbank format from SIGPROC
 * A time series composed of single precision (i.e. 4-byte) 
   floating point data
 * Photon arrival times (or events) in ASCII or double-precision 
   binary formats

Notice that the following formats which *used* to be supported are not:

 * Wideband Arecibo Pulsar Processor (WAPP) at Arecibo
 * The Parkes and Jodrell Bank 1-bit filterbank formats
 * SPIGOT at the GBT (may it RIP...)
 * Berkeley-Caltech Pulsar Machine (BCPM) at the GBT (may it RIP...)

If you need to process them, you can either checkout the "classic"
branch of PRESTO (see below), which is not being actively developed.
Or you can use `DSPSR` to convert those formats into SIGPROC
filterbank format (and/or maybe someday soon, PSRFITS search format).
You can grab DSPSR [here](http://dspsr.sourceforge.net).  If you
*really* need to get one of these machines working in PRESTO v2, let
me know and we can probably make it happen.  It will take a day or two
of porting for each backend.

The software is composed of numerous routines designed to handle three
main areas of pulsar analysis:

1. Data Preparation: Interference detection (`rfifind`) and removal
   (`zapbirds`) , de-dispersion (`prepdata`, `prepsubband`, and
   `mpiprepsubband`), barycentering (via TEMPO).

2. Searching: Fourier-domain acceleration (`accelsearch`), single-pulse
   (`single_pulse_search.py`), and phase-modulation or sideband searches
   (`search_bin`).

3. Folding: Candidate optimization (`prepfold`) and Time-of-Arrival
   (TOA) generation (`get_TOAs.py`).

4. Misc: Data exploration (`readfile`, `exploredat`, `explorefft`),
   de-dispersion planning (`DDplan.py`), date conversion (`mjd2cal`,
   `cal2mjd`), tons of python pulsar/astro libraries, average pulse
   creation, flux density estimation, and more...

Many additional utilities are provided for various tasks that are
often required when working with pulsar data such as time conversions,
Fourier transforms, time series and FFT exploration, byte-swapping,
etc.

The Fourier-Domain acceleration search technique that PRESTO uses in
the routine accelsearch is described in [Ransom, Eikenberry, and
Middleditch
(2002)](http://adsabs.harvard.edu/abs/2002AJ....124.1788R), and the
phase-modulation search technique used by search_bin is described in
[Ransom, Cordes, and Eikenberry
(2003)](http://adsabs.harvard.edu/abs/2003ApJ...589..911R).  Some
other basic information about PRESTO can be found in my
[thesis](http://www.cv.nrao.edu/~sransom/ransom_thesis_2001.pdf).  I
will eventually get around to finishing the documentation for PRESTO,
but until then you should know that each routine returns its basic
usage when you call it with no arguments.  I am also willing to
provide limited support via email or telephone (434-296-0320).

**Tutorial**: Note that in the "docs" directory there is a now a
tutorial which walks you through all the main steps of finding pulsars
using PRESTO.  This will need some small modifications given that
PRESTO can't currently process one of the example files (BCPM!).

## Getting it: 
The PRESTO source code is released under the GPL and
can be browsed or gotten from here in many different ways
(including zipped or tar'd or via git).  If you are too lazy to
read how to get it but have git on your system do:

    git clone git://github.com/scottransom/presto.git

To update it on a regular basis do

    cd $PRESTO
    git pull

and then re-make things in $PRESTO/src. 

For more detailed installation instructions, see INSTALL.

If you don't want to mess with git (which means that you will need to
re-install a tarball whenever there are updates) you can get it from
the "Download Source" link on the github page.

If you want the "classic" branch, do the following:

    git clone git://github.com/scottransom/presto.git
    cd presto
    git remote add classic origin/classic 
    git checkout -b classic origin/classic

then build as per the (old) INSTALL file.

### Development:
If you plan to tweak the code, I highly suggest that you use git and
clone the directory (or fork it using an account on github).  Read the
following "living document" on how to develop and collaborate in a
relatively sane way using git:
  http://docs.scipy.org/doc/numpy/dev/gitwash/index.html
If you plan on doing any significant development, please let me know
and I'll either add you as a developer, or we can push/pull changes
via git/github (see the "gitwash" document above).

*Code contributions and/or patches to fix bugs are most welcome!*

### Final Thoughts:
Please let me know if you decide to use PRESTO for any "real"
searches.  And if you find anything with it, it would be great if you
would cite either my thesis or whichever of the two papers listed
above is appropriate.  Thanks!

### Acknowledgements:
Big thanks go to Steve Eikenberry for his help developing the
algorithms, Dunc Lorimer for the basic code which was used to process
BCPM and WAPP data, David Kaplan for lots of help with the GBT SPIGOT
code, Jason Hessels for many contributions to the Python routines, and
a bunch of other contributions of various kinds from (alphabetical):
Anne Archibald, Cees Bassa, Slavko Bogdanov, Fernando Camilo, Paul
Demorest, Paulo Freire, Mike Keith, Patrick Lazarus, Maggie
Livingstone, Paul Ray, Paul Scholz, Ingrid Stairs, Kevin Stovall, and
for many comments, suggestions and patches!

Scott Ransom <sransom@nrao.edu>
