## Development (unreleased, since v6.0.1):
 * New module **`presto.observatories`**, holding the telescope-name to two-character ITOA
   observatory code table that `barycenter()` needs, plus `telescope_to_tempocode()` and
   `obscode()` lookups.  It mirrors `telescope_to_tempocode()` in `misc_utils.c` (which is not
   wrapped for python) and replaces three separate hand-maintained copies of that mapping: the
   `scopes` dict in `simple_zapbirds.py` and the `if/elif` chain in
   `presto.presto.bary_to_topo()` (a third copy went away with `presto.singlepulse`, below).
   The `if/elif` copies knew about only 4-5 telescopes each and used `EB` for Effelsberg where
   the C code uses `EF`; the shared table covers all 22 names that PRESTO knows about.  Note that
   these are *not* the one-character TEMPO ids or tempo2 site names that `polycos.py` uses.
 * `presto.presto.bary_to_topo()` no longer passes a stale array-length argument to
   `barycenter()` (it still needs a rewrite for its `numpy.linalg.old` import).
 * `INSTALL.md` and the top-level `meson.build` comments now document what the build flags
   actually buy: `buildtype=release` (`-O3`) is no faster than the default `debugoptimized`
   and loses the debug symbols, while `-Dc_args=-march=native` makes `accelsearch` roughly
   1.5-2x faster on a CPU with AVX2/FMA -- for local builds only, since such binaries will
   not run on older machines and their results differ at the roundoff level.
 * **Red noise removal is now a `libpresto` routine** rather than code living inside the
   `rednoise` program.  `deredden()` (in-core) and `deredden_file()` (streaming, and able to
   work in place) are declared in `misc_utils.h`.  The near-duplicate `deredden()` that lived
   in `accel_utils.c` was removed in favor of the shared one; `accelsearch` results on `.dat`
   input are unchanged.
 * `rednoise` gained **`-inplace`**, which whitens the input `.fft` in place and then renames
   it to `<root>_red.fft` (writing the usual `<root>_red.inf`).  That halves both the disk
   space and the I/O the operation needs, at the cost of destroying the original `.fft`.
 * `realfft` gained **`-rednoise`** (plus `-startwidth`, `-endwidth` and `-endfreq`), which
   removes red noise as part of the transform and writes `<root>_red.fft` and
   `<root>_red.inf` directly.  The un-whitened `.fft` is never written, so this saves a
   file-sized write and read compared to running `realfft` and then `rednoise`.  The results
   are bit-identical to running the two programs in sequence.  An input `.inf` file is
   required (the observation length is needed), and the option only applies to forward
   transforms.
 * `rednoise` bug fix: it opened `argv[1]` rather than the parsed filename, so any command
   line that put an option before the filename (e.g. `rednoise -startwidth 8 foo.fft`)
   tried to open the option instead.  Its output files now always land next to the input
   file, rather than in the current directory when the `.inf`'s recorded name has no path.
 * `realfft` now warns that `-outdir` is not honored for out-of-core (very large) transforms,
   which has always been the case.
 * conda-recipe: a round of packaging improvements from the conda-forge staged-recipes review
   (thanks to @eunos-1128, who packaged both of the new dependencies below):
   - **ERFA is no longer vendored.** A `liberfa` feedstock now exists on conda-forge, so the
     recipe just lists `liberfa` in `requirements.host` instead of fetching and building the
     ERFA tarball as a second source.  `--wrap-mode=nodownload` makes a missing `erfa.pc` a
     loud configure-time error rather than a silent fallback to the vendoring path.  PRESTO's
     `subprojects/erfa.wrap` is unchanged and still covers source builds with no system ERFA.
     The package is now plain GPL-2.0-only (no bundled BSD-3 ERFA license).
   - **`osx-arm64` builds are enabled**, now that `tempo2` is available for that platform;
     only Windows is skipped.
   - `meson setup` uses conda-forge's `${MESON_ARGS}`, which supplies the cross file needed
     for cross-compiled targets like `osx-arm64`, and `meson compile` builds in parallel with
     `-j ${CPU_COUNT}`.
   - Dropped the Fortran compiler from the build requirements: PRESTO's compiled code has been
     Fortran-free since v6.0.0, and PGPLOT supplies its own `libgfortran` runtime.
   - Added a `package_contents` test that checks the installed tools, `libpresto`, `presto.pc`,
     and the runtime data files.

## Version 6.0.1:
 * A small release cut mainly so that the first conda-forge package ships the CHIME support in
   `get_TOAs.py`, plus the conda-forge recipe updates requested during the staged-recipes review.
 * `get_TOAs.py` now knows about **CHIME** (telescope code `y` for TEMPO-style TOAs, observatory
   name `chime` for tempo2/Princeton-style output), so TOAs can be generated from CHIME folds
   without hand-editing the telescope tables.  The script was also re-formatted with ruff.
 * Converted the conda-forge recipe to the **v1 (rattler-build / CEP-13) recipe format**, at the
   request of the `conda-forge/staged-recipes` reviewer: `conda-recipe/meta.yaml` is now
   `conda-recipe/recipe.yaml` (`${{ }}` jinja, a `context:` block, `if:`/`then:` selectors, and a
   `tests:` list).  Also switched the `astropy` run dependency to `astropy-base` (PRESTO only uses
   core astropy), as the conda-forge linter suggested.
 * conda-recipe: request `gsl >=2.8` explicitly.  conda-forge's global pin is still gsl 2.7 while
   the gsl 2.8 migration is in flight, but tempo2 has already been rebuilt against 2.8, so a
   build against 2.7 produced a package that could not be co-installed with tempo2.  Remove the
   constraint once the global pin moves.

## Version 6.0.0:
 * **First conda-forge release, and a major dependency cleanup.** PRESTO v6 no longer
   requires the `$PRESTO` environment variable, no longer needs TEMPO at all (barycentering
   is done in-process with ERFA, and polycos come from tempo2), and its compiled code is now
   Fortran-free (FFTFIT is pure Python; the last Fortran least-squares solver moved to GSL).
   Runtime data and man pages install into the prefix, so `man rfifind` and the pulsar
   catalog work with no extra setup. Together these make PRESTO cleanly installable from
   conda-forge (`conda install -c conda-forge presto-pulsar` &mdash; the `presto` name is
   taken on bioconda by an unrelated tool).  See the bullets below for the individual changes.
 * Added release/packaging infrastructure: a `RELEASE.md` checklist, a `--set X.Y.Z` mode for
   `determine_version.py` (to stamp an explicit release version, not just the git-derived dev
   one), and an in-repo conda-forge recipe (`conda-recipe/meta.yaml` + `build.sh`) that builds
   both the C and Python halves and vendors ERFA (which has no conda-forge feedstock).
 * `prepfold` now prints a one-line hint when its optional Ghostscript-based `.png` conversion
   fails (previously it was skipped silently).  On recent Ubuntu the usual cause is the AppArmor
   `gs` profile, which blocks `gs` from reading/writing outside `$HOME`, `/tmp`, and `/var/tmp`
   (so folding in `/dev/shm` or a cluster scratch dir yields the `.ps` but no `.png`).  Added a
   FAQ entry explaining the cause and the fixes.
 * Fixed a `free(): invalid pointer` crash in `rfifind -nocompute` that occurred at exit, right
   after printing "Done."  The per-thread `tmprfi_arr` was `malloc`ed but only initialized to
   `NULL` inside the compute branch, so with `-nocompute` the cleanup loop's null-check ran on
   uninitialized pointers and freed garbage.  It is now `calloc`ed.
 * Fixed spurious "FFTW wisdom not up-to-date" warnings in conda/pixi environments.  The
   installed executables (and `libpresto`) now carry an RPATH to `<prefix>/<libdir>`, so they
   load the same FFTW that ships in the install prefix instead of falling through the loader
   cache to a possibly-different system `libfftw3f`.  Previously `makewisdom` could write wisdom
   with one FFTW build while the tools (and the Python `_presto` extension, loaded via the conda
   `python`'s own RPATH) read it with another; because FFTW's wisdom version cookie is
   per-build, the mismatched reader rejected the file.  The absolute RPATH is rewritten to a
   relocatable `$ORIGIN` form by conda-build on conda-forge.  Also reworded that FFTW-import
   warning in `fftcalls.c`, which previously mislabeled any rejected wisdom as merely "not
   up-to-date".
 * Removed the long-dead `python/binresponses/` directory (`monte_ffdot`, `monte_short`,
   `monte_sideb`, `montebinresp`).  Nothing referenced it, it was never part of the build, and
   the code relied on the Python-2-era `Numeric`, `mpi`, `miscutils`, and `Statistics` modules,
   so it could no longer even import.
 * Fixed the remaining NumPy 2.5 `DeprecationWarning`s about setting `array.shape` in place by
   switching to `np.reshape()` throughout the Python code (`rfifind`, `psr_utils`, `filterbank`,
   `injectpsr`, `mpfit`, the `presto_src` wrapper, `singlepulse/rrattrap`, and the `detrend_dat`,
   `downsample_filterbank`, `pulsestack`, `single_pulse_search`, and `jerk_example` scripts).
 * **Modernized the main `python/presto` modules** (per ROADMAP "Docstrings, typehints, and
   ruff"): added NumPy-style docstrings and type hints and ran `ruff format` on `psr_utils`,
   `events`, `bestprof`, `harmonic_sum`, `infodata`, `filterbank`, `polycos`, `spectra`,
   `sinc_interp`, `psrfits`, `simple_roots`, `prepfold`, `parfile`, `rfifind`, `binary_psr`,
   `cosine_rand`, `sigproc`, and the `presto_src` SWIG-wrapper layer.  Along the way, dropped the
   Python-2 `from builtins import ...` shims and the now-unused **`six`** dependency, renamed
   module-local `import numpy as Num/num` to `np`, and replaced a few `import *` wildcards with
   explicit imports.  Notable bug fixes uncovered during the pass:
   - `psr_utils`, `events`: fixed APIs so the modules run on modern NumPy/SciPy (`np.complex`,
     `np.array(copy=False)`, `scipy.optimize.zeros`, `np.outerproduct`, `np.float`, and a
     missing `scipy.integrate.quad` import — several `events` functions could not run before).
   - `spectra`: the `Spectra` constructor now stores its `dm` argument (was silently discarded).
   - `infodata.to_file()`: now writes the `waveband`/`beam_diam` fields (were lost on a
     read/write round-trip).
   - `parfile.ELL1_check()` and `sinc_interp.kaiser_window()`: fixed undefined-name bugs
     (`Num.sqrt`, a missing `i0` import) that stopped them running.
   - `rfifind.write_weights_and_offsets()`: `invertband=True` is no longer silently ignored.
   - `binary_psr`: `reflex_motion()` referenced an undefined `omega`, and `most_recent_peri()`
     raised on array input.
   - `psrfits`: fixed an invalid-escape `DATE-OBS` regex, a wrong multi-file `TELESCOP` check,
     and a misplaced `np.any()` parenthesis (all warning/display paths).
   - `cosine_rand`: load its shipped JSON table via `importlib.resources`.
 * Fixed the stale `__version__` in `python/presto/__init__.py` (was `'v5.0.2.dev22'`) to match
   the authoritative `5.3.1`, and taught `determine_version.py --write` to keep it in sync
   (adding `truncate()` so a shorter version string can't corrupt the rewritten files).
 * Removed the vendored **Starlink SLALIB** date/time routines (`src/cldj.c`,
   `src/djcl.c`, `include/slalib.h`, `include/slamac.h`) and replaced their use with
   **ERFA**, which is already a hard dependency (barycentering).  `slaCldj` (Gregorian
   calendar &rarr; MJD) is now `eraCal2jd`, reached through a new libpresto helper
   `cal_to_mjd()` in `misc_utils.c` (used by the PSRFITS/SPIGOT/BPP/WAPP readers and
   the `cal2mjd` tool); `slaDjcl` (MJD &rarr; calendar) is now `eraJd2cal` (in
   `mjd_to_datestr()` and the `mjd2cal` tool).  While in `misc_utils.c`,
   `sphere_ang_diff()` now wraps `eraSeps()` (results agree to ~4e-15 rad), and
   `dms2rad()`/`hms2rad()` now wrap `eraAf2a()`/`eraTf2a()` so ERFA owns the tricky
   sign conventions (e.g. negative declinations) — output is bit-identical to the old
   homegrown arithmetic.  No behavior change for any tool.
 * Removed the vendored **RANDLIB** random-number library (`src/randlib.c`,
   `src/com.c`, `include/randlib.h`) and replaced its use in `makedata` with GSL's
   generators (`gsl_rng` seeded from the clock/PID, `gsl_ran_gaussian` for Gaussian
   noise, `gsl_ran_poisson` for Poissonian noise).  GSL is already a hard dependency,
   so no new dependency is introduced.  RANDLIB was only used by `makedata` and the
   unbuilt legacy test `tests/test_ffts.c` (also ported to GSL).  `makedata` output is
   random by construction, so this changes the specific noise realizations but not
   their statistics.
 * Replaced the Numerical-Recipes-derived Quickselect in `src/median.c` with a thin
   wrapper around GSL's `gsl_stats_float_median()` (GSL is already a hard dependency).
   The `median()` API is unchanged, so all callers and the build files are untouched.
   Results are bit-identical for odd-length inputs; for even-length inputs the value
   is now the mean of the two central elements (GSL's interpolation) rather than the
   lower one — a statistically negligible change in the power-normalization contexts
   where `median()` is used.
 * `prepfold`'s `.pfd.ps` plots are now converted to anti-aliased `.png` images
   using **Ghostscript** (`gs`) called directly, replacing latex2html's `pstoimg`.
   latex2html is large and, crucially, not available on conda-forge; `pstoimg` was
   only a wrapper around `gs`, so the output is identical (`-dTextAlphaBits`/
   `-dGraphicsAlphaBits` give the anti-aliasing, `Orientation 3` reproduces the old
   `-flip cw`, and `-dEPSCrop` reproduces `-crop a`).  `prepfold` now skips the
   `.png` silently if `gs` is missing or fails — the all-important `.ps` is always
   written.  The helper scripts (`bin/pfd2png.sh`, the GBNCC/GBT350/PALFA survey
   scripts, `examplescripts/jerk_example.py`), the `INSTALL.md`/Dockerfile
   dependency lists, and `pixi.toml` were switched from `latex2html` to
   `ghostscript` to match.
 * General cleanup: removed many outdated and stray files, including obsolete
   standalone programs (`src/searchmultidms.c`, `src/showmulti.c`,
   `src/showmulti_dm.c`, `src/testtree.c`), a batch of legacy C test programs along
   with the old `tests/Makefile`, stale data tables (`lib/goodfactors.txt`,
   `lib/parkes_birds.txt`), and old documentation (`docs/presto.tex` and its
   figures).
 * `makewisdom` now defaults to a fast `FFTW_MEASURE` pass that skips the two
   largest FFT sizes in each of its size loops (a few seconds rather than
   ~10-20 min); pass `-patient` for the previous full `FFTW_PATIENT` sweep over
   all sizes for higher-quality wisdom.
 * Replaced the Numerical Recipes `amoeba.c` downhill-simplex with GSL's
   Nelder-Mead simplex (`gsl_multimin_fminimizer_nmsimplex2`) in the f-fdot
   optimizers `max_rz_arr()` and `max_rz_arr_harmonics()` (`src/maximize_rz.c`),
   and removed `src/amoeba.c`.  Validated with a before/after golden harness (now
   committed as `tests/test_maximize_rz.py`): the results match the old solver to
   well below the FFT-plan noise floor on 171/172 synthetic cases, and on the
   remaining case GSL finds the true peak where amoeba got stuck at its initial
   simplex vertex (a candidate started exactly on a harmonic peak).  The GSL
   objective also guards against non-finite values so the simplex is steered away
   from out-of-band harmonics rather than aborting.
 * **PRESTO no longer requires the `PRESTO` environment variable to be set** — a big
   change after ~26 years!  Shared runtime data (the pulsar catalog `pulsars.cat`,
   `psr_catalog.txt`, and `aliases.txt`) is now installed to `{prefix}/share/presto`,
   and both the C tools and the Python package locate it with the same search order:
   `$PRESTO/lib` (still honored as an optional override, e.g. for running from a
   source tree) -> the compiled-in `{prefix}/share/presto` -> a path derived from the
   running executable/module (so relocated conda/pixi environments still work).  This
   is implemented via `presto_data_path()` in the new `src/datadir.c` and
   `data_path()` in the new `python/presto/_datadir.py`.  The man pages (`docs/*.1`)
   now install to `{prefix}/share/man/man1`, so `man rfifind` (etc.) works with no
   extra setup (conda activation already puts that on `MANPATH`).
 * `makewisdom` is now installed and writes `fftw_wisdom.txt` directly into the
   PRESTO data directory (`{prefix}/share/presto`, or `$PRESTO/lib` if set) where the
   tools look for it at runtime, so no manual copying is needed.  An optional path
   argument overrides the destination.
 * Removed the long-obsolete `a2x` ASCII-to-PostScript helper (`bin/a2x.sh`,
   `lib/a2x.ps`, `lib/a2xshell.ps`) and the dead code that used it: the unused
   `src/cand_output.c` routines (`file_reg_candidates`/`file_bin_candidates`) and the
   already-broken legacy test `tests/test_apps.c` (plus its `.mak` files).
 * Polyco generation (e.g. `prepfold -timing`) now uses **tempo2** (which is on
   conda-forge) rather than TEMPO, via `tempo2 -tempo1 -polyco`, which writes
   TEMPO1-format polycos.  Together with the ERFA barycentering below, **PRESTO no
   longer needs TEMPO at all**.  Validated against TEMPO-generated polycos: constant
   phase offsets of only 0.001-0.023 microsec, differential drift <= 3e-6 turns,
   |df|/f <= 4e-11, and identical TOAs from `get_TOAs.py`; the comparison can be
   re-run anytime with `check_polycos_t2` (needs both tempo and tempo2).  The old
   generator is retained as `make_polycos_tempo1()`.  This also fixed two latent
   site-code bugs: ATA polycos were generated for SHAO 65m, and KAT-7 polycos for
   FAST (the tempo2 observatory names are unambiguous).
 * Barycentering no longer uses TEMPO!  `barycenter()` is now computed in-process
   with the [ERFA](https://github.com/liberfa/erfa) library (a new build requirement,
   downloaded and built automatically by meson if not installed), using ERFA's
   built-in analytical ephemeris, leap seconds, and a new internal observatory
   table (`src/observatories.c`).  Agreement with the old TEMPO-based code is at
   the ~10 microsec level (absolute; the ephemeris difference) and ~1 microsec
   differentially across an observation, with v/c matching to ~3e-9.  The `ephem`
   argument (e.g. "DE405") is now ignored.  The old implementation is retained as
   `barycenter_tempo()` for validation (see `src/check_bary_erfa.c` and
   `tests/test_barycenter.py`).  TEMPO is still used for polycos.
 * New program `prepfold_multi`: fold many candidates from a single raw-data pass
   (reads/dedisperses once, folds all candidates, OpenMP-parallelized), built on a
   large refactor of the prepfold internals into a shared pipeline
   (`prepfold_pipeline.c`).  Thanks to Ryan Lynch for the initial refactor work!
 * The compiled code is now completely Fortran-free:
   * FFTFIT was rewritten in pure Python (NumPy/SciPy), with two selectable methods:
     the faithful `classic` port (default) and Anne Archibald's independent
     `aarchiba` algorithm (`get_TOAs.py -A/--fftfit {classic,aarchiba}`).  Validated
     against a frozen golden dataset from the original Fortran.
   * The vendored LAPACK `dgels` least-squares solver used by `bary2topo()` was
     replaced with GSL (`gsl_multifit_linear`), and `fminbr.c` with GSL's Brent
     minimizer.
 * Command-line parsing specs (`clig/*.cli`) are now processed with the new `pyclig`
   generator (the old `clig` is unmaintained); removed the copies of old
   clig-generated man pages which lived in the `clig` directory.
 * Filled in `Description`/`Author`/`SeeAlso` content for the `clig/*.cli` files, so
   `man <tool>` is now real documentation instead of pyclig's placeholder text, with
   usage examples for the main routines (`rfifind`, `prepdata`, `prepsubband`,
   `mpiprepsubband`, `accelsearch`, `realfft`, `show_pfd`, `prepfold`). All CLIG-based
   files were regenerated except `prepfold_multi`, which keeps its existing
   hand-written man page. Removed the unused `dftfold` program, and copied
   `prepfold_multi`'s `-candfile` format into a new `FAQ.md` entry for discoverability.
 * Improved accuracy of chi-squared CDF / significance calculations via better
   transitions to asymptotic expansions in `chi2_logp()`.
 * Fixed `stacksearch.py` FFT-stacking bugs (FFTs shorter/longer than the stack) and
   added a simple sifting script for its outputs over multiple DMs.
 * Fixed prepdata/prepsubband output length when using `-start`/`-offset` without
   `-numout`, fixed the Geocenter observatory code for polycos (issue #221), and
   fixed several pre-existing prepfold memory leaks.
 * The prepfold regression-test suite moved to `tests/prepfold/` with a
   tolerance-based `.bestprof` comparator.
 * Many compiler warnings fixed (thanks Erum Vohra!), plus workarounds for spurious
   gcc 8.5 / gcc 14 warnings; updated Dockerfiles (thanks Wenke Xia!); removed the
   long-stale Travis CI setup.
 * Added `CLAUDE.md` and `ROADMAP.md` describing the codebase and development plans.

## Version 5.3.1:
 * Turned on OpenMP linking for rfifind! Oops.

## Version 5.3.0:
 * rfifind now gets speedup using multiple CPUs via OpenMP (Thanks to Erum Vohra!)
 * Speed and parallelization improvements to accelsearch (Thanks to Erum Vohra!)
 * Bug fix in reading of EPS2DOT in readpar()
 * Bug fix to accelsearch when searching datasets longer than 2^32 points that caused segfaults
 * Bug fix to prepdata, prepsubband, and mpiprepsubband when barycentering files longer than 2^31 points that caused segfaults
 * Bug fix in prepfold where wrong fdotdot falues were reported if folding with fdotdot
 * Turned off OpenMP in the dedispersion routines until I figure out how to handle that better (it works in prepfold, though!)

## Version 5.2.0:
 * Now require GSL (Gnu Scientific Library) to be installed. This library is thread-safe and allows us to more easily parallelize routines.
 * Large set of changes that updated K&R-style declarations so that they would compile with GCC v15. These changes were made by Claude Code(!) with the prompting by Paul Ray. Thanks, Paul and Claude!

## Version 5.1.0:
 * Updated ATNF Pulsar Catalog to v2.65
 * Three new and useful python programs / utilities:
   * `stacksearch.py` Read multiple PRESTO-style `*.fft` files and conduct a stack search for periodicities.
   * `fourier_fold.py` Use the complex amplitudes in a PRESTO `.fft` file (or in multiple FFT files) to generate pulse profiles without having to do time-domain folding.
   * `pfdzap.py` Perform simple time- and/or frequency domain zapping of `.pfd` files. Generate zap commands for `show_pfd`, `get_TOAs.py`, and `sum_profiles.py`.
 * Many small bug fixes and tweaks, including more correct handling of DM smearing in `DDplan.py`

## Version 5.0.3:
 * Updated ATNF Pulsar Catalog to v2.51
 * Added an experimental version of `fit_circular_orbit` using sliders in `examplescripts`
 * Added the ability to use `pygaussfit.py` without middle or right mouse buttons
 * Fixed a couple memory issues in `rednoise`, thanks to @bwmeyers
 * Explicitly set the random number seeds (for reproducibility) in `makedata`
 * Several other very minor tweaks and bug fixes

## Version 5.0.2:
 * Updated the C wrappers for PGPLOT for the Numpy 2.0 C API 
 * Python v3.9 or newer is now required.
 * Several minor bug fixes, including to `injectpsr.py`, thanks to @remsforian

## Version 5.0.1:
 * Minor improvements over v5.0.0
 * Some clarifications and improvements to the build process
 * Addition of new recipes for docker / singularity images (thanks Alessandro Ridolfi!)
 * Bugfix in `rednoise` having to do with absolute paths (thanks Alessandro Ridolfi!)
 * Bugfix that checks to see if the maskfile you intend to use matches the properties of the data you are using

## Version 5.0.0:
 * This is a major release since I've moved to a completely different and modern build system: [meson](https://mesonbuild.com/), along with the [meson-python](https://meson-python.readthedocs.io/en/latest/) backend. This was required since *Numpy* has deprecated `numpy.distutils` and this caused python builds to stop working with Python v3.12.
   * See the [INSTALL.md]()https://github.com/scottransom/presto/blob/master/INSTALL.md) for updated installation instructions.
   * You will need to install **meson**, **meson-python**, and **ninja**, but that is easily done via `pip`!
   * Python v3.8 or newer is now required.
 * All of the old Spigot-related codes have been removed. If you need to process Spigot data, please use the `classic` branch that is mentioned in the README.md.
 * All of the `slalib` codes (and the python interface to it) have been removed. If you need that stuff, you should transition to [ERFA](https://github.com/liberfa/erfa) and/or [Astropy](https://www.astropy.org/).
 * There are two nice new python utilities:
   * `binary_utils.py` reads a parfile of a binary pulsar and computes min/max observed barycentric spin periods or velocities as either a function of the orbit (default), or for a prescribed duration of time, and optionally plots those. It also shows basic information about the binary.
   * `compare_periods.py` compares candidate spin periods and their integer and fractional harmonics with one or more parfiles to a prescribed fractional tolerance. 

## Version 4.0:
 * This is a major release since it involves big changes to the Python portions of the codebase:
   * Python v3.7 or newer is now required.
   * A long-standing memory issue was fixed with Anaconda Python (running `python tests/test_presto_python.py` will tell you if you have that issue or not).
   * Swig v4 is used to generate the Python wrappers of the PRESTO C library.
   * Big thanks to **Shami Chatterjee** and **Bradley Meyers** who helped me get to the bottom of this!
 * There is a [FAQ](https://github.com/scottransom/presto/blob/master/FAQ.md) with lots of information!
 * PRESTO has a dockerfile that allows it to build on Docker Hub automatically. Thanks to **Nick Swainston** for this!  (more testing and improvements would be welcome)
 * `simple_zapbirds.py` makes it much easier to manually zap interference from simple searches (no need for copying ".inf" files and running both `makezaplist.py` and `zapbirds`).
 * `realfft` and `zapbirds` can now be called on many files at once on the command line. This benefits HPC systems which often don't like many programs running serially on many small files.
 * A new python interface to the internal `prepfold` folding code (`simplefold`), as well as wrappers of fast `C` implementations of $\chi^2$ and $Z^2_N$ (thanks to **Matteo Bachetti**).
 * Many bug fixes and minor improvements, including one that would cause segfaults with very large dispersion sweeps in `prepdata` and `prepsubband`, and a problem with `prepfold` significance calculations.

## Version 3.0.1:
 * This is a minor release which fixes several issues and adds some minor improvements:
   * Fix of long-standing `rfifind` bug that could cause the program to hang if channels had zero variance
   * Multiple Python3-related bug fixes
   * Added `-debug` flag to `prepfold` to allow debugging of TEMPO calls to make polycos
   * `DDplan.py` can now read observation parameters from filterbank or PSRFITS input files. And you can write a `dedisp_*.py` dedispersion script, based on the plan, using the `-w` option
   * The `rednoise` program now writes a corresponding *_red.inf file
   * Update of the Tutorial document, including a new slide on red noise

## Version 3.0:
 * This major release of PRESTO includes a massive restructuring of python code and capabilities. Things should work with Python versions 2.7 and Python 3.6 and 3.7 at least. The installation of the python code has changed and has become more "pythonic" so that `PYTHONPATH` is not needed, and all of the various modules are now under a top-level "presto" module. For example, to use the psr_utils module you would now do:
   
   `import presto.psr_utils as pu`
   
   rather than

   `import psr_utils as pu`

   All of these changes will likely lead to code breakage and bugs!

   Please check your code and processing carefully and post issues (and hopefully pull requests) if you find them.

   The installation instructions have been updated in the INSTALL file.

   Huge thanks thanks go to **Gijs Molenaar, Matteo Bachetti, and Paul Ray** for the work that they have done helping with this!

 * There is also a new `examplescripts` directory where you will find some example code to do a lot of important things, like
   * Fully dedispersing an observation: `dedisp.py`
   * Fully searching a dedispersed observation: `full_analysis.py`
   * Sifting the results of a full search: `ACCEL_sift.py`
   * Searching short chunks of a long time series: `short_analysis_simple.py`
   * Making a really nice P-Pdot plane: `ppdot_plane_plot.py`
   * and a few others.

## Version 2.2:
 * Version 2.2 was the last version of PRESTO to work with the old-style python interface which requires Python v2.7 or earlier and is "installed" in-place and used via having `$PRESTO/lib/python` in your `PYTHONPATH`. There will probably be occasional bug fixes for v2.2 in the `v2.2maint` branch of PRESTO. You can get it using: `git checkout -b v2.2maint origin/v2.2maint`, and then installing as per the INSTALL file.

## Version 2.1:
 * `accelsearch` now has a "jerk" search capability (thanks to (then) UVA undergrad **Bridget Andersen** for help with this!). This makes searches take a *lot* longer, but definitely improves sensitivity when the observation duration is 5-15% of the duration of the orbital period.  Typically `-wmax` should be set to 3-5x `-zmax` (and you probably never need to set `-zmax` to anything larger than 300).
 * Ability to ignore bad channels on the command line (-ignorechan) (see `rfifind_stats.py` and `weights_to_ignorechan.py`)
