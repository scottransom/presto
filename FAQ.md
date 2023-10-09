# PRESTO FAQ

## Introduction

I've had many questions over the years about PRESTO. This FAQ is an attempt at
getting some of them -- and their answers -- down on digital paper. If you have
any additions, please let me know!

Scott Ransom <sransom@nrao.edu>

-----------------

## General PRESTO questions

-----------------

### **I've read the tutorial, but there seem to be a lot of other programs in `$PRESTO/bin`. What are some useful ones and what do they do?**

Here are a few. Remember that for most PRESTO routines, if you run the command
with no additional options, it gives you a brief usage statement.

- `cal2mjd` and `mjd2cal`: These compute a UTC MJD from a calendar date and
  time, or vise-versa.

- `psrfits_quick_bandpass.py`: If you have PSRFITs searchmode data, this routine
  will quickly compute (and plot if requested) the average and standard
  deviation of the bandpass of your data. You can ignore the `DAT_SCL` and
  `DAT_OFFS` if you want, as well (i.e. to check bit occupation).

- `rfifind_stats.py`: Do a kind of integrated statistics on `rfifind` result
  files, including creating files that have an ASCII bandpass (".bandpass"), the
  channels recommended to zap (".zapchans"), and recommended weights
  (".weights") for each channel.

- `weights_to_ignorechan.py`: Convert a ".weights" file from `rfifind_stats.py`
  into a compact format that can be passed to PRESTO routines using the
  `-ignorechan` flag, or to PSRCHIVE's `paz` routine.

- `combine_weights.py`: If you use the above routines to get multiple ".weights"
  files per observation (i.e. by running `rfifind` on subsets of a long
  observation), this will combine the files (via logical "or") to make a
  "combined.weights" file. 

- `toas2dat`: If you have events (e.g. photon arrival times), this simple
  routine will bin those into a time series (i.e. a ".dat" file) that can be
  processed by PRESTO.

- `show_pfd`: Do various things to a `prepfold` ".pfd" file, including zap
  interference as a function of time (`-killparts`) or freq (`-killsubs`),
  partially fixing corrupted $\chi^2$ values due to interference (`-fixchi`),
  and make publication-quality phase-vs-time plots (`-justprofs`).

- `pygaussfit.py`: A interactive gaussian fitter that reads in a ".pfd.bestprof"
  file and outputs text that can be stored in a ".gaussians" file, which can be
  used as a template to get TOAs with `get_TOAs.py`.

- `pyplotres.py`: An interactive timing residuals viewer that works with
  original TEMPO.

- `sum_profiles.py`: A routine that will let you correctly sum various ".pfd"
  files, including rotating them so that they are aligned, to make a high
  signal-to-noise profile. Can also be used to give radiometer equation
  estimates of flux density for the same profiles.

- `dat2sdat` and `sdat2dat`: Not used much anymore, but a way to shrink the size
  of ".dat" files by a factor of 2, by saving them as short integers (and
  vise-versa).

- `gotocand.py`: Would need to be modified for your own computer system(s), but
  a way to much more easily fold `accelsearch` candidates, even if they are
  potentially on other nodes (via `mpiprepsubband`). The `-local` option is
  useful if you are just using a single machine.

- `orbellipsefit.py`: Get an initial orbit fit using the ellipse-fitting
  technique as described in [Freire, Kramer, & Lyne
  2001](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..885F/abstract). Note
  that this does not coherently connect observations via the orbital period.

- `fit_circular_orbit.py`: Fit a circular orbit to data using ".bestprof" files
  (or daily ".par" files) as input, where the orbit is coherently connected
  between observations. This has helped to solve *many* binary pulsars!

- `fitorb.py`: Similar to the above routine, but it can fit eccentric orbits if
  needed (once again uses ".bestprof" or daily ".par" files as input).

- `quickffdots.py`: Quickly plot a small portion of the f/fdot plane in an
  ".fft" file (including summing up to 4 harmonics).

- `quick_prune_cands.py`: A quick sifting script for individual accelsearch
  result files. It is quite good at throwing out non-pulsars.

- `pfd_for_timing.py`: Returns true or false for a ".pfd" file to let you know
  if it can be used to get TOAs or not (i.e. was it used to search for the best
  p/pdot or not).

- `makedata`: A crude (but effective!) routine that lets you generate a time
  series with simulated pulsations. This was the very first routine in PRESTO,
  and my first routine in "C" -- the code is particularly gross. :-)

-----------------

### **Does PRESTO use GPUs to speed up any of the routines?**

There are currently two different versions of GPU-accelerated `accelsearch`,
although to my knowledge, neither of them do jerk searches.

I believe that they can both be used as drop-in replacements for the regular
`accelsearch`, but with everything else coming from the current, standard
PRESTO.

From Jintao Luo: https://github.com/jintaoluo/presto_on_gpu

From Chris Laidler: https://github.com/ChrisLaidler/presto

And slightly further afield, there is a new GPU implementation of the Fourier
Domain Acceleration Search here:
https://github.com/AstroAccelerateOrg/astro-accelerate

-----------------

### **I've read the tutorial, but this is just too complicated.  Is there any easier way to run this software and find pulsars?!**

Yes!  Thanks to Alessandro Ridolfi, there is [PULSAR
MINER](http://alex88ridolfi.altervista.org/pagine/pulsar_software_PULSAR_MINER.html),
which uses PRESTO under the hood, but with a mostly automated pipeline that you
set up with a simple configure script.  It can even use Jintao Luo's
GPU-accelerated `accelsearch` by default.

PULSAR MINER has been used to find many of the newly discovered globular cluster
MSPs from MeerKAT.

-----------------

### **Many of these routines are really slow, and they don't seem to get faster using the `-ncpus` option?  Why is that?**

Most of PRESTO was written before multi-core processing was really much of a
thing.  And so it was designed to be mostly serial, with parallelism coming via
independent searches of different DMs at the same time (which is why
`mpiprepsubband` was written: to speed-up de-dispersion and to distribute the
time series among many processors/nodes).

I've added some OpenMP pragmas over the years, and they make some small
improvements in, for example, dedispersion loops.  But the improvement is not
great.  For `prepdata`, `prepsubband`, and `rfifind`, I'd recommend not using
more than 3-4 CPUs currently with `-ncpus`.  And even then, you won't see
anything close to a 3-4x speedup (probably only a couple tens of percent).

`accelsearch` is the exception, and it does fairly well with `-ncpus`, although
the performance got much worse after the addition of the jerk-search code.

In summary, the OpenMP stuff is very much a work in progress, and I'd love to
work with someone on this if there is any interest.  I suspect that significant
speed-ups could be had without a ton of work.  It is on my ToDo list!

-----------------

### **I'm interested in using `mpiprepsubband` to speed-up dedispersion.  How exactly is it used?**

One of the biggest problems with dedispersion is the I/O bottleneck -- you are
taking data in one big raw-data file and producing hundreds or even many
thousands of time series. Most computer systems don't do well with many
parallel output streams, and so it can take a *long* time to write all of those
dedispersed files.

One way to help the output is to split it over many different machines. That
way, each machine gets a fraction of the time series and there is less
bottleneck on each machine. This is what `mpiprepsubband` does. A single CPU
reads the raw data, it distributes it over the network to different nodes, and
each node (using many cpus on each node) then dedisperses a fraction of the DMs.

In order to translate a dedispersion plan (as generated by `DDplan.py` using the
`-s` subband option) into a proper call for `mpiprepsubband`, there are several
things to remember:
 - You will probably not see much of an advantage unless you dedisperse using
   multiple different machines. Many cores on the same machine does not solve
   the I/O problem.
 - You need one CPU on a machine that can see the raw data. And then N *other*
   CPUs on each of M other machines/nodes to do the dedispersion. That means a
   total of M * N + 1 CPUs. And the bigger M is (i.e. number of dedispersing
   nodes), the better your I/O performance will be.
 - A single `mpiprepsubband` call is like a full line of the `DDplan.py` output:
   each of the M * N dedispersing CPUs is equivalent to an independent
   `prepsubband` *call*.
 - The total number of DMs from an `mpiprepsubband` run must be divisible by M *
   N (i.e. each CPU generates the same number of output DMs).
 - `DDplan.py` knows about these constraints and will use them if you use the
   `-p numprocs` flag.

-----------------

### **I have some data in some format XXX. I want to search it with PRESTO. How do I do that?**

If your data have multiple frequency channels, you first need to integrate them
(with de-dispersion, if necessary) in order create a 1-D time series. A very
good way to do that, which will let you use almost all of PRESTO's tools is to
convert the data to the relatively simple SIGPROC filterbank format (there are
some tools in PRESTO to help with that -- see below).

If your data are events, you need to put them into a binned 1-D time series if
you are going to search them with `accelsearch`. There is a PRESTO utility for
that called `toas2dat`. Note that `prepfold` can fold events directly using the
`-events` option, if you specify their format and have a ".inf" file that
describes the observation.

Finally, your 1-D time series need to be saved in binary form as 32-bit floating
point numbers. You should give it a ".dat" suffix. You then need an ASCII ".inf"
file to go with it. You can use the PRESTO tutorial and de-disperse one of the
test datasets to see what a radio ".inf" file looks like and simply edit it for
your data.

Note that you probably need $10^4$ (at least) points in a time series to make it
worthwhile to search with `accelsearch`. And even then, you need to be very
careful about the frequency limits `-flo` and `-fhi` if they are outside the
normal pulsar-search regime!

------------------

### **I have some baseband data (perhaps from a Software-defined Radio or SDR setup), can I get that into a format for PRESTO to process?**

I'd recommend that you convert your data to the simple (except for the
header...ugh) [SIGPROC](http://sigproc.sourceforge.net) filterbank format after
"detecting" and channelizing.

To get data in SIGPROC format, you will need to write a SIGPROC header, and then
simply append the binary-format channelized Stokes I (i.e. total intensity) data
right after that.

The [DSPSR](http://dspsr.sourceforge.net) software suite does have the
capability to do some of these things (i.e. channelizing and detecting baseband
and writing it to SIGPROC or PSRFITS format), but it might not know about the
exact format of data you are using. The routine in DSPSR you would want is
called `digifil`.

Alternatively, there is some python code in PRESTO that you can hack so that you
can write a filterbank header and then stream the data into the file from
`numpy` arrays (for instance). The code is called `sigproc.py` in
`$PRESTO/python/presto/`. And there is some related code in the same directory
called `filterbank.py`.

-----------------

### **Do the p and pdot (or f and fdot) values returned by PRESTO refer to beginning of the dataset or the middle of the dataset?**

Depends on the search tool:

* For `accelsearch`, the f, fdot, (and fdotdot, if requested) returned are
  average values during the observation, so they apply to the midpoint of the
  dataset.

* For `prepfold`, the reference time is always the beginning of the observation.

-----------------

### **Can I get the average barycentric velocity of an observatory during an observation in a nicer way than running (e.g.) `prepdata` on the raw data file?**

Yes, there is a python convenience function available in the 2nd-level `presto`
module (i.e. `from presto import presto`):

    get_baryv(ra, dec, mjd, T, obs="PK"):
      Determine the average barycentric velocity towards 'ra', 'dec'
      during an observation from 'obs'. The RA and DEC are in the
      standard string format (i.e. 'hh:mm:ss.ssss' and 'dd:mm:ss.ssss').
      'T' is in sec and 'mjd' is (of course) in MJD. The obs variable
      is the standard two character string from TEMPO:  PK, GB, AO, GM, JB, ...

Note:  there is also a general-purpose barycentering code you can run called
`bary`, where you feed it topocentric UTC MJDs in ASCII from STDIN (and give it
a single '.' to tell it you are done).

-----------------

## RFI-related questions

-----------------

### **Are zapping with `.birds` files and masking using `rfifind` the only type of RFI removal in PRESTO?**

RFI is handled multiple ways in PRESTO, most of which are optional. There is a
pretty good discussion about them in [Lazarus et al
2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...812...81L/abstract).

The quick summary is that RFI is handled in at least **six** different ways:

1. Any processing of a raw data file (i.e. SIGPROC filterbank or PSRFITS
    search-mode) by default has clipping turned on, which will zap zero-DM
    short-duration pulses or data dropouts. It is controlled with the flag
    `-clip` which sets the threshold time-domain S/N to clip, or turned off
    completely with `-noclip`. **I highly recommend that you not turn this
    off!** It is almost always useful and will not harm (the vast majority of)
    dispersed astrophysical pulses.
2. `rfifind` finds and masks short duration and/or narrow band periodicities or
    time-domain statistical anomalies in the data. There are many options. You
    can also generate channel zaplists using `rfifind_stats.py` and
    `weights_to_ignorechan.py`, which can used handled by the various `prep*`
    routines.
3. Red-noise suppression in the power spectrum of de-dispersed data (i.e.
    ".fft" file) using the `rednoise` command.
4. ".birds"-file zapping of known periodic signals during Fourier searches via
    `zapbirds`, `simple_zapbirds.py`, or `accelsearch` itself (if searching
    ".dat" files directly). This is used to zap broad-band periodic signals such
    as known pulsars, power line modulation, etc, directly in the power spectrum
    (i.e. ".fft" file).
5. `ACCEL_sift.py` can remove anomalous candidates during search sifting (user
    configurable).
6. You can use `-zerodm` processing (see [Eatough et al,
    2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.395..410E/abstract)) after
    running `rfifind` and using the resulting mask files with the `prep*`
    routines. Be careful with this as it definitely removes power from
    mildly-dedispersed pulsar signals! (But it can also be extremely useful.)

-----------------

### **What is the difference between using `-ignorechan` and explicitly including channels that you want to zap in an `rfifind` mask using `-zapchan`? Is one preferred over the other?**

`-ignorechan` is a recent-ish addition to PRESTO, and, in general, gives
better performance for most cases *if* know know for certain what channels
you want to zap for the full observation. You pass the list of bad
channels to `rfifind` and all of the other raw-data routines in PRESTO
(i.e. `prepdata`, `prepsubband`, `prepfold`, and `mpiprepsubband`), and
those routines then *completely ignore* those channels. If those channels
have bad data (i.e. they are at band edges or something like that), this
will allow `rfifind` to do its job on the channels where there is real
signal and real noise.

The difference between using `-ignorechan` and `-zapchan` is that the
former completely ignores those channels in all processing (as if they
were pure zeros), while the latter, when used with an `rfifind` mask, will
effectively replace the channel values with a smooth-ish running median
(as determined from the `rfifind` "stats" file) for that channel over the
observation.

If the data are well behaved, then they should both do about the same
thing. However, if the input data are badly behaved, for example if the
power levels change a lot during the observation, and especially if the
statistics of the channels in question are highly variable, then masking
the data can end up leaving low-frequency artifacts in the resulting time
series or folds. And those often show up in your searches or folds and are not good.

If you don't know what channels might be bad, and which you might want to
use `-ignorechan` on, you could do a quick first pass with `rfifind` on a
portion of the data (for instance) to identify the bad channels, and then
make an `-ignorechan` list from that. That's exactly what `rfifind` +
`rfifind_stats.py` + `weights_to_ignorechan.py` does. And if you are
lucky, you can simply use the resulting `-ignorechan` and its values and
not even use an `rfifind` mask.

The way I use `-zapchan` (which is only rarely) is almost always after I
do a search (or after closely looking at my `rfifind` mask results) where
I see that some periodic RFI is leaking through into certain channels.  So
I change the mask (via `rfifind -nocompute`) to explicitly add those
channels (they may have been partly masked but didn't pass the `-intfrac`
cutoff to have them zapped completely).

-----------------

### **I'm seeing strong 60 (or 50!) Hz signals in my `accelsearch` results which are obviously from the power mains. Why doesn't `rfifind` filter that out?**

**If** the 50Hz signal is strong enough to show up in individual frequency
channels for relatively short periods of time, then rfifind will flag it (and
mask away that channel). However, usually power mains signals are weaker and
show up across the band (i.e. they are too weak to be detected in a single
channel after only integrating for a few seconds). We usually zap those in the
FFT searches using a ".birds" file (with the `zapbirds` or `simple_zapbirds.py`
commands). That zeros out portions of the power spectrum where there are known
RFI signals (and their harmonics!). Run `simple_zapbirds.py` with no command
line options for more information.

This is discussed a little bit in the PRESTO tutorial, and in more detail in
[Lazarus et al
2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...812...81L/abstract).

-----------------

### **My time series has a lot of low-frequency noise in it. Can I remove that in some way?**

Yes, at least to get a pretty good removal of the rednoise, try this:
1. Prepare and de-disperse your data as per usual to get a ".dat"/".inf" file
    pair.
2. Run `realfft` on the ".dat" file to get a ".fft" file.
3. Run `rednoise` on the ".fft" file. That will create a "_red.fft" file that
    will have been de-reddened.
4. If you simply want to search the data, you can search the "_red.fft" file.
5. If you want a de-reddened time series, perhaps for folding with `prepfold`,
    run `realfft` on the "_red.fft" file, which will inverse-FFT it and create a
    new "_red.dat" file that will be the de-reddened time series.

This tends to work quite well!

-----------------

### **How are periodic signals searched for in `rfifind`?**

Each small chunk (in time) of data from each frequency channel is FFT'd,
converted to powers, and then normalized before being searched by the routine
`search_fft()` which is in `minifft.c`. The normalization constant is computed
from the standard deviation of the time series (and the number of points) for
the small chunk of data that is being FFT'd.

`rfifind` then takes the result of that search and checks if any of the powers
in the short FFT (in one interval in one channel) are above the equivalent
gaussian significance `-freqsig` if you search an FFT of that length (including
trials factors).

-----------------

### **Can I get or change the information in the `rfifind` masks in Python?**

Yes. There is a `rfifind` module in PRESTO's Python tools. That at least lets you read and play with mask and statistics values from `rfifind`. It is not currently possible to write "mask" files, but that could change (and would not be difficult). Here is an example of usage:

    In [1]: import presto.rfifind as r
    
    In [2]: r = r.rfifind("GBT_Lband_PSR_rfifind.mask")
    
    In [3]: r.nint
    Out[3]: 37
    
    In [4]: r.read_mask()
    
    In [5]: shape(r.mask_zap_chans_per_int)
    Out[5]: (37,)
    
    In [6]: r.mask_zap_chans_per_int
    Out[6]:
    [array([56, 59, 94], dtype=int32),
     array([18, 32, 53, 94], dtype=int32),
    ...
     array([44, 50, 94], dtype=int32),
     array([94], dtype=int32)]

Note that the first `r.rfifind()` call does a `read_mask()` automatically. So
that is already available in `r.mask_zap_chans_per_int` as soon as you load the
first file.

-----------------

### **Does rfifind use the same frequency channel numbering as the original raw data file?**

So PRESTO *always* processes data so that the lowest frequency channel is
channel 0. If the bandwidth is reversed in a filterbank or PSRFITS file, it will
flip it in memory automatically.

You can tell if this is happening by running `readfile` on the file. If the BW
is lower sideband (i.e. needs flipped) you will see:

           Invert the band? = True

If you see that, you likely will need to change your channel ordering if you
want to zap channels using other software (i.e. PSRCHIVE).

You can view the band and channel mask by using the command `rfifind_stats.py
MYFILE_rfifind.inf`

-----------------

## Python modules questions

-----------------

### **Do you have any Python code to compute [SOMETHING ABOUT PULSARs]?**

Maybe! I highly recommend you browse the code in `$PRESTO/python`, and
especially the highly useful `psr_utils.py`. I normally load that with: `import
presto.psr_utils as pu`.

You can also read in `rfifind` result files using the `presto.rfifind` module
and `prepfold` files using the `presto.prepfold` module. The ATNF pulsar
database is available in `presto.pypsrcat`. And there are tools for reading
filterbank and PSRFITs files into PRESTO as well (i.e. `psrfits.py`,
`filterbank.py`, and `spectra.py`). You can read TEMPO residuals using
`residuals.py`.

-----------------

## Pulsar searching / `accelsearch` questions

### **What are reasonable ranges to use for `-zmax` and/or `-wmax` if I want to find binary millisecond pulsars?**

This was discussed a bit near the end of section 2.2 in [Andersen & Ransom
2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...863L..13A/abstract). I ran a
bunch of simulations on MSPs around stellar-mass-type companions, and the vast
majority can be detected with `-zmax` < 200 or so, and `-wmax` < 600, *if* the
orbital period is longer than about 10 times observation duration.

Note that all acceleration and jerk effects affect higher harmonics of a pulsar
signal more than they do the fundamental. So you can typically drop the number
of harmonics you search (`-numharm`) down to 4 or so if you are looking for
really accelerated systems.

If you are looking for massive companions (i.e. big black holes), then the
`zmax` and `wmax` values can be quite a bit larger (i.e. in the thousands). The
problem there lies in the fact that as you go to very large `zmax` values, it is
likely that the acceleration isn't constant during the observation, and so you
*need* to use some jerk searching to compensate. The same thing happens for the
constant jerk assumption in jerk searches, but PRESTO doesn't have "snap"
searches yet... ;-)

### **What is a *sigma* in `accelsearch` and what does it mean?**

In short, *sigma* is the probability that a given signal might be due to noise,
but expressed in terms of equivalent gaussian sigma, despite what the true
probability distribution for noise values is. In other words, it is shorthand
for a probability.

The way that sigmas are computed in searches is a tricky thing. And even in an
individual `accelsearch` "cands" file, the sigma term means different things
depending where you are looking in the file.

The "summary" candidate list at the top of the file is from the search through
the full f/fdot plane, where the Fourier powers are normalized by block running
medians.

In the bottom, where the harmonics of each candidate are analyzed, each harmonic
is individually maximized in the f/fdot plane, and the average local power level
is measured using powers around the harmonic but not including the 5 closest
Fourier bins on each side.

If the data are statistically very well behaved and if the harmonic is not super
strong, so that the sinc-function sidelobes don't mess up the local power
computation, then the normalized powers measured in both ways should be quite
similar.

But if there is a lot of RFI or red noise or some other reason why the power
spectrum isn't "nice" (as you would get from pure and unchanging gaussian noise
in the time domain), then the two different normalizations might be quite
different.

The normalization used in the search is faster, while the normalization used for
optimization is "more correct" (since it doesn't at all use the frequency in
question, but only nearby frequencies).

But that is just the power. The sigmas as calculated from the powers have
another important difference:
* The significance at the top of the files in the summary candidate list
  **includes** a correction for the number of independent frequencies/fdots
  searched (**for that single run of accelsearch**, i.e. not including others DM
  trials, for instance). That is the so-called "look-elsewhere" effect.
* The significances at the bottom of the candidates file (where you get detailed
  information about each harmonic of each candidate), assume that you are
  searching only a single Fourier frequency (i.e. single trial, so there is no
  trials correction). That means that for a single-harmonic search, even with
  pure gaussian data, you would see different significances top and bottom.

Note, however, that `ACCEL_sift.py` is made for comparing results between
different searches (i.e. over DMs, for instance), and so the sigmas that it
returns are all single-trial, so that they can be compared!

Finally, as to what the sigmas mean, in general, they are the equivalent
gaussian significance (i.e. meaning as if the powers were gaussian distributed,
which they are not, they are $\chi^2$ distributed) that you would see that power
sum (with or without the trials factor corrections as mentioned above). In other
words, PRESTO computes the CDF of the appropriate $\chi^2$ distribution (the
numbers of degrees of freedom are 2 times the numbers of harmonics summed) and
then converts that to a probability. That probability is then quoted in sigma
as if the PDF were gaussian and not $\chi^2$. That's because most people are
used to thinking about "sigma" as being from a gaussian distribution.

-----------------

### **How can one obtain the spectral signal-to-noise (S/N) from a given sigma value for an accelsearch candidate?**

I actually recommend that you don't do that. I think that S/N in the Fourier
domain (and especially in powers) is not a really useful quantity. You can't
directly compare the S/Ns of different numbers of harmonics, for example. But
you *can* compare their equivalent gaussian significances. That is precisely why
I use sigma (i.e. probabilities) and not S/N.

However, if you need to do it for some reason, basically, the power spectrum is
normalized so that the mean noise values have a level of one. But these are
powers. S/N should be in amplitude. So you would take the normalized power
level of a harmonic, take the sqrt of it (to convert it to an amplitude), and
subtract 1 from the result (since the mean is 1). You can sum the S/N of each
harmonic together to get a rough idea of the S/N of the signal as a whole. But
once again, beware that you cannot properly compare the S/N of a sum of 16
harmonic powers with the S/N of a sum of 8 harmonics, for instance. They have
very different probability distributions ($\chi^2$ with 32 and 16 degrees of
freedom, respectively).

-----------------

### **How do I change the candidate threshold in `accelsearch`? And how is it defined?**

If you simply run `accelsearch` you will see all of the default options. One of
those is `-sigma`, with a default value of 2. That number is the equivalent
gaussian significance of a candidate, using the initial median-block
normalization (or normalization via number of photons if using `-photon`), and
including the number of independent trials searched for that time series.

Given that that is already quite low, you probably don't want to go much below
1.5-sigma or so, or you'll start getting a lot of false positives.

-----------------

### **I've restricted the frequency range of my search using `-flo`/`-fhi` or `-rlo`/`-rhi`, but `accelsearch` is still returning candidates outside of that range! Why is that?**

Those options definitely work, however, they correspond to the frequencies of
the *highest harmonic searched*. The subharmonics of those frequencies will
also get searched, and those subharmonics might be below the frequency of
`-flo`. If you are searching only with a single harmonic (`-numharm 1`), then
the frequency ranges work exactly as you would expect.

The reason `accelsearch` does this is that harmonic summing is done based on the
highest harmonic in order to save operations and memory. If we search a region
around frequency *N*, that corresponds to a single-harmonic search. If we add
to that the frequencies around *N/2*, then it is a two-harmonic search. To get
a four harmonic search, we only need to then add the frequencies around *N/4*
and *3N/4*. And all of these subharmonics are smaller "chunks" of the f/fdot
plane that can be interpolated to exactly correspond to the full-resolution
search around frequency *N*.

Note that the frequency ranges also change the number of independent trials
searched, and so they *will* affect the sigmas returned by `accelsearch`.

-----------------

### **Is there any documentation on how to do sideband modulation analyses?**

No. Very simply, you run `search_bin` on ".fft" files, and manually examine the
resulting candidate files. If there are good candidates, you can try to get a
phase-coherent model of the signal using the `bincand` command (which is quite
time-consuming). But both of these codes are really stale. I don't think I've
run them in 8+ years, so there could easily be bugs. It has been on my to-do
list to update for a long time...

For reference, the details as to what is going on in these codes is described in
[Ransom, Cordes, & Eikenberry,
2003](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..911R/abstract).

-----------------

## `prepfold` questions

-----------------

### **What is all of this chi-square and reduced chi-squared stuff and what does it mean?**

`prepfold` uses the reduced chi-squared (or $\chi^2_{red}$) as a way of
determining how significant pulsations are in your "folded" data (meaning when
you take data and co-add it modulo an assumed pulsation period).  The technique
has been around for quite a while and is often known as "epoch folding" [(See
section IIIb in Leahy et. al 1983 for
details).](https://ui.adsabs.harvard.edu/abs/1983ApJ...266..160L/abstract)

Basically, this is like a normal $\chi^2$ goodness-of-fit statistic, but in this
case, the model that we are fitting to the data is a constant one, i.e. no
pulsations.  We compute the $\chi^2$ of the folded profile compared to the
average value of the folded profile.  The resulting statistic (assuming gaussian
random data in the time series) is $\chi^2$ distributed with Nbins-1 Degrees of
Freedom, where Nbins is the number of bins in the pulse profile.

Because of the exact way that `prepfold` folds binned data, the statistics are
slightly more complicated than that, in reality (see below), but the bottom line
is that the larger the $\chi^2_{red}$ is, the more likely it was that noise
fluctuations didn't cause it (i.e. it is due to real pulsations).  There are
better and alternate statistics we could use, but $\chi^2$ is quick and simple
and quite sensitive to narrow pulse profiles, which we often have in radio
pulsar searches.

Note that for the statistics to be reasonably correct (in the face of strong
RFI, for instance), the off-pulse (i.e. away from the periodicity in question)
$\chi^2_{red}$ noise floor should be approximately 1.  That makes sense
because in that part of parameter space there should be no pulsations, and so
the no-pulsations model should fit the data well and give a $\chi^2_{red} \sim 1$.

-----------------

### **In a `prepfold` plot of a search candidate, there is a significance provided by the chi-square. Does that include the number of trials searched?**

No, the $\chi^2$ from `prepfold` is single-trial significance. And it is only
valid if the mean of the off-pulse (i.e. away from the periodicity in question)
$\chi^2_{red}$ noise floor is approximately 1.

-----------------

### **What is a *sigma* in `prepfold` and what does it mean?**

Just as for *sigma* in `accelsearch`, *sigma* in `prepfold` is the probability
that a given pulsed signal might be due to noise, but expressed in terms of
equivalent gaussian sigma, despite that the true probability distribution for
`prepfold` trials is $\chi^2$. In other words, it is shorthand for a
probability.

As described above and below, `prepfold` uses reduced chi-squared to determine
pulsation significance (which, if all due to noise, would be distributed as a
$\chi^2$ distribution with Nbins-1 degrees of freedom), and so the probability
is based on that number.

For real pulsar signals, that probability can be a tiny number, so I convert it
to the equivalent gaussian significance in sigma (i.e. meaning as if the
significance distribution was gaussian distributed rather than $\chi^2$) since
that is nicer numerically and we are often used to thinking of things in terms
of gaussian significance.

-----------------

### **Why is the number of degrees of freedom a fraction and why does it have a subscript "EFF"?**

The folding algorithm in `prepfold` is different than in many other pulsar
folding codes. Instead of assuming that each datapoint is a delta function in
time, and therefore corresponding to an infinitely short portion of a pulsar's
rotation (which can be put all in one pulse profile bin), `prepfold` assumes
that the data point is finite in duration (i.e. integrated in time), which most
data actually are, especially search-mode data which PRESTO is focused on.

`prepfold` knows the instantaneous start and end phases of each input data bin,
and it effectively "drizzles" the content of that data bin uniformly over as
many pulse profile bins as are needed to contain it in phase. If the duration of
the data bins are similar to or longer than the duration of the pulse profile
bins, the "drizzling" causes significant correlations between neighboring pulse
profile bins. And that leads to effectively smoother profiles (i.e. less RMS
from the noise) and fewer effective degrees of freedom if you use $\chi^2$ for
significance, as `prepfold` does.

There is a semi-analytic correction for this effect that has been tested using
large numbers of simulations with fake data. You can see it in the `DOF_corr`
function/method in `src/fold.c` or `python/presto/prepfold.py`, respectfully.
That correction can also be used to correct for the noise level in the pulse
profile, for more accurate flux densities estimates via the radiometer equation,
for example (`sum_profiles.py` uses the correction; newrms = oldrms /
sqrt(DOF_corr)).

The effect and the correction are described in detail in Appendix E of the
recent paper [Bachetti et al.,
2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...33B/abstract).

-----------------

### **What are all of these options in `prepfold` for? How do I know what I should be using?**

- In general, if you are folding candidates from an `accelsearch` run:
  - Use the `-accelcand #` and `-accelfile FILE` options to specify the basic
    p/pdot/pdotdot values of the fold, directly from the search.
  - If you are folding raw data, make sure you also specify the `-dm`!
  - If you find that your candidate is actually a harmonic (or if you want to
    check), you can adjust the folding parameters automatically to account for
    that using `-pfact` or `-ffact`, which are multiplicative factors on the
    period or spin frequency.

- If you are folding topocentric time series data in order to get
  quick-and-dirty TOAs for solving a timing solution, for example, make sure
  that you fold with `-nosearch` and `-n` which is a power-of-two so that you
  will be able to get TOAs using `get_TOAs.py`.

- For folding raw data, in general, I recommend setting `-nsub` to be between
  50-250, with good values for most pulsar instruments being 64 or 128 (since
  most have powers-of-two channels). That is high enough frequency resolution
  to remove narrow band interference, if needed, but coarse enough frequency
  resolution so that you are integrating enough signal to see your pulsar in
  each subband! The default `-nsub` tries to combine roughly 8-10 frequency
  channels into each subband.

- Also for folding raw data, remember that you can use an `rfifind` mask file
  using `-mask` or specify channels to ignore via an `-ignorechan` option. Those
  can *really* help with bad RFI. You can also use `-zerodm` to fold with
  zero-DM subtraction.

- If you have a particularly faint pulsar candidate, the `-fine` option can help
  make sure that you latch on to it, rather than any other nearby noise peaks,
  since it searches a smaller region of p/pd/pdd/DM volume.

- The `-coarse` option doubles the normal size of each dimension of the
  p/pd/pdd/DM volume, and so is useful for searching for brighter pulsars where
  you don't have an exact ephemeris and so need to search a bit more.

- If you fold a slow pulsar (for most pulsar searches, that would be probably
  periods greater than 100 ms, or so), I recommend that you use the `-slow`
  option. That specifies the `-fine` option, which helps prevent you from
  latching on to nearby RFI, and it also gives you 100 profile bins (i.e. `-n
  100`). Since most slow pulsars don't show much acceleration, you should also
  probably use `-nopdsearch`, which will also help avoid latching on to
  interference. If you don't have a lot of dispersion, I also recommend using
  `-nodmsearch`, or you can easily latch onto DM=0 interference.

- If you have a good timing-based parfile (i.e. ephemeris) for your pulsar, you
  should probably fold with `-timing`. That provides multiple advantages:
  - It folds with polycos and doesn't do any searching. That allows you to get
    TOAs from the `.pfd` file using PRESTO's `get_TOAs.py` or PSRCHIVE's `pat`.
  - It automatically sets `-dm` based on the parfile value.
  - It makes sure that the number of profile bins is a power of two (needed for
    PRESTO's implementation of FFTFIT to get TOAs).
  - It sets the number of intervals to 60 by default, which is highly factorable
    (and allows you to generate a variety of different numbers of TOAs using
    `get_TOAs.py`).
  - It uses `-fine`, which makes (IMO) the final plots look nicer.
  - Note that if you want to use the absolute phase information in a parfile
    (which might be able to keep all of your pulse profiles aligned over many
    days), you can use the `-absphase` option. You can't currently get good
    TOAs using data folded that way, though. If you need this functionality,
    please let me know.
  - Folding with `-par` only uses polycos to fold, and doesn't do any of the
    other things. So you probably don't want to use that ever.

- If you are just doing a very rough and non-scientific fold of some data (i.e.
  a test pulsar), you can use the `-psr` option and specify a pulsar name.
  `prepfold` will grab the best information for the pulsar from the ATNF Pulsar
  Database and fold with that (allowing searching). You *cannot* use `.pfd`
  files folded in that way for TOAs.

- For some PSRFITs data, if you fold with `-noscales` and `-nooffsets` you
  effectively get a bandpass-flattened fold. That can often be useful.

- Folding with binary parameters specified is almost never needed. If you have
  a binary, you should be folding with `accelsearch`-determined parameters from
  a search, or with `-timing` if you have a good ephemeris.

- You should only use `-searchpdd` if you are not also searching in DM!

- You normally do not need to specify `-psrfits` or `-filterbank` unless your
  filename is bizarre.

- If you notice after the fold that interference or some other issue has made
  the off-candidate $\chi^2$ values far from 1 (meaning you get bad statistics
  and pulse significance, and sometimes see drops in the accumulated $\chi^2$
  curve), you can attempt a post-facto fix of that using `show_pfd -fixchi
  MYFILE.pfd`

-----------------

### **There is some RFI in my `prepfold` plot, can I get rid of that post-facto?**

You can use the `show_pfd` routine to re-generate the `prepfold` plots and
statistics. That routine has `-killparts` and `-killsubs` options for zapping
individual, or ranges, of intervals in time, or subbands in frequency.

In the `prepfold` plots, the numbering of the intervals and subbands is shown on
the right side of the phase-vs-time and phase-vs-frequency greyscale plots.

As an example, if you want to zap channels 10, 13-17, 20-30, and 94 from you plot,
the command would be:
```
> show_pfd -killsubs 10,13:17,20:30,94 MYFILE.pfd
```
Note that ":" replaces "-" for the range, and there are no spaces! And also
remember that the 1st channel/subband is 0 and not 1!

-----------------

### **How reliable are the uncertainties in the best period and period-derivative(s) as determined by `prepfold`?**

The answer is that they can be pretty good if your data is well-behaved (i.e.
gaussian noise, constant background levels, no interference, etc). But that
isn't like most radio data, obviously.

In general, those errors come from combining the uncertainties for each fourier
harmonic of the signal (via equations derived by John Middleditch, and available
in [Ransom, Eikenberry, &
Middleditch](https://ui.adsabs.harvard.edu/abs/2002AJ....124.1788R/abstract)).
But they seem to underestimate the uncertainty for real signals a bit --
especially really bright ones. It helps if your stepsize in the P/Pdot plane is
very fine (use the `-fine` option, for instance, and use more bins in the pulse
profile).

De-reddening the time series data also helps as that makes the noise more like
gaussian white noise (which the error derivations assumed).

Bottom line: they are likely underestimated by between 20-100% for real data,
although very recently (June 2021) I pushed up a commit related to the
correlation of `prepfold`'s profile bins that likely fixes most of this issue.
I'm planning to test it with fake data in the future.

-----------------

### **Can I work with the folding data cube from `prepfold` (i.e. the `.pfd` file) with Python?**

Yes! The `prepfold.py` module has a bunch of methods to let you do many
interesting things with your ".pfd" files. Here is an example usage:

    In [1]: import presto.prepfold as pp
    
    In [2]: a = pp.pfd("GBT_Lband_PSR_4.62ms_Cand.pfd")
    
    In [3]: a.use_for_timing()  # can we use this pfd file for TOAs?
    Out[3]: False
    
    In [4]: a.dedisperse(16)  # dedisperse at a DM of 16 pc/cm^3
    
    In [5]: a.sumprof  # return the summed profile of the dedispersed data
    Out[5]: 
    array([5656790.31402834, 5652502.51116502, 5654345.94100014,
           5656388.48898718, 5656145.69576171, 5655103.75782315,
           5656093.92149403, 5654931.8004717 , 5654154.6155577 ,
           5655499.99197552, 5658468.12322909, 5658051.62727781, ...

There are lots of plotting and analysis options. Just take a read through the
`python/presto/prepfold.py` file, especially all the docstrings for the `pfd`
class methods and attributes.

Also, note that PSRCHIVE can load in and work with some ".pfd" files!

-----------------

### **In the "Optimizing" portion of a `prepfold` fold/search, it says "`Warning!:  This is 6535218305 trials! This will take forever!`", and it is! How should I do this?**

When you fold raw data for a search candidate, by default it searches over
nearby DMs, periods, and p-dots ("pd"). That is a 3-dimensional search space,
where each dimension is a few times the number of bins in the pulse profile.
That is a lot of trials. If you then add searching over period second
derivatives ("pdd", i.e. `-searchpdd`) that makes things orders of magnitudes
worse since it is now a much bigger 4-D search.

What you should do if you get a good search candidate is to find the DM where
the search signal was maximized and do a search over p, pd, and (if needed) pdd
on the *time series*. Once that is optimized, and the candidate looks good in
the time series fold, *then* fold the raw data using `-nopsearch` and
`-nodsearch` and only allow it to search over DM. You can specify the best
`-pdd` on the command line if it was needed, but don't specify `-searchpdd`.

Note that you can read the ASCII `.pfd.bestprof` file from the time series fold
to get the optimized best values for `-p`, `-pd`, and `-pdd` instead of reading
them off the `prepfold` plot.

-----------------

### **With `pat` in PSRCHIVE, the left edge of the template is always phase=0, what does `get_TOAs.py` use for phase=0?**

By default, `get_TOAs.py` uses the old "Princeton" method of rotating the
template, so that the phase of the fundamental harmonic of the template is set
to be at zero phase. So if you have a perfectly Gaussian template, then TOAs
will be referenced at its peak, which would be phase zero (since the fundamental
harmonic peaks at phase zero). This method lets you automatically do
"pretty-good" alignment of similar templates even if you don't manually
rotate/align them yourself.

You can get the same "don't rotate my template" behavior of `pat` if you use the
`--norotate` (or `-r`) option to `get_TOAs.py`. And the TOAs generated by `pat`
and `get_TOAs.py`, in that case, should be nearly identical.

-----------------

### **I'm getting bogus TOAs from this observation. Why is that?**

It is likely that you folded allowing searching, and `get_TOAs.py` doesn't have
enough information to properly determine the correct spin parameters and timing.

Try running `pfd_for_timing.py` on your ".pfd" file. If it returns True, you
should be able to get good TOAs. If not, you probably need to re-fold using
`-nosearch`.

-----------------

## De-dispersion or `prepdata`/`prepsubband` questions

-----------------

### **For de-dispersion, where in a frequency channel does PRESTO assume that the frequency for that channel corresponds? Bottom, center, or top of the channel?**

PRESTO does the same as SIGPROC and most other pulsar software and assumes that
the frequency corresponds to the center of the channel. In addition, PRESTO
never shifts the highest frequency channel during de-dispersion, but only brings
the lower frequency channels "forward" in time.

Note that this also means that the MJD Epoch mentioned in a `.inf` file
corresponds to the UTC time of the start of the first sample of the center of
the highest frequency channel.

-----------------

### **How does PRESTO decide what block-size or chunk-size to process raw data?**

This has changed several times in the past. Since PRESTO processes things in
terms of blocks (i.e. number of samples in time), it is useful to have that
number be a big enough chunk of data that you can do useful work with it. But
you also want it to be small enough to give you useful granularity in `rfifind`
or in the `prepfold` subints in time.

For PSRFITS data, the minimum block size is the duration of each PSRFITS SUBINT
row. But for SIGPROC filterbank files, it is hardcoded in PRESTO. If you have
a need (and there are some special reasons why you might), it is possible to
change that size. It is in `$PRESTO/src/sigproc_fb.c` in the line:

    s->spectra_per_subint = 2400;        // use this as the blocksize

You then need to re-run "make".

-----------------

### **I'm trying to downsample my data using `prepdata` or `prepsubband` and it is not letting me use certain downsampling factors. Why is that?**

Your downsampling factor must be evenly divisible into the `Spectra per subint`
for your data, which can be seen using the `readfile` command.

That number is the duration of each PSRFITS row, or for SIGPROC filterbank data,
is hard-coded (see above question) at 2400 (which is highly factorable):

    ❯ factor 2400
    2400: 2 2 2 2 2 3 5 5

So you can use `-downsamp` of 2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 25, ...
etc, which is pretty good.

If you use `DDplan.py` to plan your de-dispersion by pointing it at your raw
data file, it will only give you downsample factors that your data supports.

-----------------

### **I'd like to read the ".dat" and/or ".fft" files into Python. What is the easiest way to do that?**

PRESTO's ".dat" files are straight 32-bit floats (i.e. `numpy.float32`, or
typecode "f"), and the ".fft" files are straight 2x32-bit floats treated as
complex numbers (i.e. `numpy.complex64`, or typecode "F"). This means that you
can easily read them using `numpy.fromfile`:

    In [1]: import numpy as np
    
    In [2]: dat = np.fromfile("myfile.dat", dtype=np.float32)
    
    In [3]: fft = np.fromfile("myfile.fft", dtype=np.complex64)

For the FFT data, the zeroth frequency bin (i.e. `fft[0]`) has the real-valued
DC component in the real portion, and the real-valued Nyquist frequency stored
as the imaginary portion! That is a common trick for packing FFT results so that
the number of complex data points is N/2, where N was the length of the input
time series, rather than N/2+1.  The FFT amplitudes are also completely
unnormalized (equivalent to using the "raw" normalization mode with
`explorefft`).

Note that you can read in the information in the related ".inf" files using either:

    In [4]: from presto import presto
    
    In [5]: inf1 = presto.read_inffile("myfile.inf")
    Reading information from "myfile.inf"
    
    In [6]: inf1
    Out[6]: <presto.presto.prestoswig.infodata; proxy of <Swig Object of type 'INFODATA *' at 0x7f8061d1ac00> 

which uses PRESTO's C-code and structures, as wrapped by `swig`, or:

    In [7]: from presto import infodata
    
    In [8]: inf2 = infodata.infodata("myfile.inf")
    
    In [9]: inf2
    Out[9]: <presto.infodata.infodata at 0x7f80c9d2d400>

which is pure Python (and therefore much easier to modify, if needed).

-----------------

Please let me know if you have other things that you think should go in this file!

Scott Ransom <sransom@nrao.edu>
