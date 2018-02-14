# PRESTO

http://www.cv.nrao.edu/~sransom/presto/

PRESTO is a large suite of pulsar search and analysis software
developed by Scott Ransom mostly from scratch, and released under the
GPL (v2).  It was primarily designed to efficiently search for binary
millisecond pulsars from long observations of globular clusters
(although it has since been used in several surveys with short
integrations and to process a lot of X-ray data as well).  It is
written primarily in ANSI C, with many of the recent routines in
Python.  According to Steve Eikenberry, PRESTO stands for: PulsaR
Exploration and Search TOolkit!

**PRESTO has discovered over 600 pulsars, including more than 230
recycled and/or binary pulsars!**

## New in Version 2.1:
 * `accelsearch` now has a "jerk" search capability (thanks to UVA
   undergrad Bridget Andersen for help with this!).  This makes
   searches take a *lot* longer, but definitely improves sensitivity
   when the observation duration is 5-15% of the duration of the orbital
   period.  Typically -wmax should be set to 3-5x -zmax (and you probably
   never need to set -zmax to anything larger than 300).
 * Ability to ignore bad channels on the command line (-ignorechan)
   (see `rfifind_stats.py` and `weights_to_ignorechan.py`)
 * Lots of new python utilities (such as for handling RFI, showing
   bandpasses, making waterfall plots, ...)
 * New wrappers for the python interface (will make the transition
   to Python 3.X much smoother later this year)
 * Many bug fixes and minor improvements

## About PRESTO:
PRESTO is written with portability, ease-of-use, and memory efficiency
in mind, it can currently handle raw data from the following pulsar
machines or formats:

 * PSRFITS search-format data (as from GUPPI at the GBT, PUPPI and
   the Mock Spectrometers at Arecibo, and much new and archived data
   from Parkes)
 * 1-, 2-, 4-, 8-, and 32-bit (float) filterbank format from SIGPROC
 * A time series composed of single precision (i.e. 4-byte) 
   floating point data (with a text ".inf" file describing it)
 * Photon arrival times (or events) in ASCII or double-precision 
   binary formats

Notice that the following formats which *used* to be supported are not:

 * Wideband Arecibo Pulsar Processor (WAPP) at Arecibo
 * The Parkes and Jodrell Bank 1-bit filterbank formats
 * SPIGOT at the GBT
 * Berkeley-Caltech Pulsar Machine (BCPM) at the GBT

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
5. Post Single Pulse Searching Tools: Grouping algorithm (`rrattrap.py`),
   Production and of single pulse diagnostic plots (`make_spd.py`,
   `plot_spd.py`, and `waterfaller.py`).

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

**Tutorial**: Note that in the "docs" directory there is a tutorial
which walks you through all the main steps of finding pulsars using
PRESTO.

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
clone the directory (or fork it using an account on github).  And if
you want to contribute your changes back, please give me a "pull
request"!

*Code contributions and/or patches to fix bugs are most welcome!*

### Final Thoughts:
Please let me know if you decide to use PRESTO for any "real"
searches.  And if you find anything with it, it would be great if you
would cite either my thesis or whichever of the two papers listed
above is appropriate.  Thanks!

### Acknowledgements:
Big thanks go to Steve Eikenberry for his help developing the
algorithms, Dunc Lorimer and David Kaplan for help with (retired) code
to process BCPM, SCAMP, and Spigot data, Jason Hessels for many
contributions to the Python routines, and (alphabetical): Bridget
Andersen, Anne Archibald, Cees Bassa, Matteo Bachetti, Slavko
Bogdanov, Fernando Camilo, Paul Demorest, Paulo Freire, Chen Karako,
Mike Keith, Patrick Lazarus, Maggie Livingstone, Chitrang Patel, Paul
Ray, Paul Scholz, Ingrid Stairs, Kevin Stovall, Joeri van Leeuwen for
many comments, suggestions and patches!

Scott Ransom <sransom@nrao.edu>
