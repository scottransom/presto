# PRESTO Roadmap

Rough notes on significant changes planned for the coming weeks/months. This is a living,
aspirational document — items may change, be reordered, or be dropped. For stable "how to work
in this repo" guidance, see [CLAUDE.md](CLAUDE.md); for tracked work, see GitHub issues.

## Primary goal: get PRESTO into conda-forge

Package all of PRESTO — compiled code **and** Python modules — for conda-forge, so it can be
installed trivially with `conda`/`pixi`. Most of the items below are prerequisites or
simplifications that make clean packaging possible.

**ERFA packaging wrinkle:** conda-forge has **no C-library ERFA package** (verified 2026-07:
no `erfa`/`liberfa` feedstock exists, and `pyerfa` only ships a private Python extension — no
`erfa.h` or `liberfa.so` to link against; Meson WrapDB has no erfa entry either). PRESTO
handles this generally via `subprojects/erfa.wrap`: meson uses system ERFA when present and
otherwise downloads/builds the release tarball automatically (offline: pre-place the tarball
in `subprojects/packagecache/`). Since conda-forge builds have no network access, the PRESTO
recipe should either (a) list the ERFA tarball as an extra source unpacked into
`subprojects/packagecache/` (works today, but statically vendors it), or (b) — cleaner —
submit a trivial `liberfa` feedstock to conda-forge first (autotools or meson build, ~30-line
recipe) and depend on it. Prefer (b).

### Remove environment-variable requirements

Stop requiring `PRESTO` (and other) environment variables to be set. These currently locate
runtime data such as the pulsar database and the FFTW wisdom file, among other things. For a
conda-forge package these paths need to be discoverable without user-set env vars.

**[DONE]** `PRESTO` is no longer required. Shared runtime data (`pulsars.cat`,
`psr_catalog.txt`, `aliases.txt`) now installs to `{prefix}/share/presto`, and both C and
Python resolve it with the same search order — `$PRESTO/lib` (optional override) →
`{prefix}/share/presto` (compiled in via `-DPRESTO_DATADIR` for C, `sys.prefix` for Python)
→ a path derived from the running executable/module (relocation-safe). New `src/datadir.c`
`presto_data_path()` (used by `database.c`/`fftcalls.c`) and `python/presto/_datadir.py`
`data_path()` (used by `pypsrcat.py`) implement this; `check_meson_build.py` now treats
`PRESTO` as optional. FFTW wisdom is machine-specific and no longer shipped (soft-warns if
absent); `makewisdom` is now installed and writes `fftw_wisdom.txt` straight into that same
data directory (via `presto_data_writepath()`), with an optional path-override argument.
Validated: catalog reads (4303 psrs) and `pypsrcat` import/regen work with `PRESTO` unset;
`test_presto_python.py` passes unset.

### Reduce/remove the Fortran dependency

- Remove Fortran from the main compiled binaries, replacing it with GSL code where possible.
- **[DONE]** Convert the original **FFTFIT** Fortran code and its Python module to pure Python
  (NumPy + SciPy). This can be done independently of the rest and should get thorough tests.
  See the dedicated plan below.
- **[DONE]** Replace the vendored LAPACK/BLAS `dgels` least-squares solver (`src/least_squares.f`,
  the last Fortran source) with GSL's `gsl_multifit_linear` in `bary2topo()`. This drops
  `least_squares.f`, `include/f77.h`, and the Fortran language from the meson build, so
  `libpresto` and all compiled tools are now Fortran-free (PGPLOT remains an external Fortran
  dependency). Validated: GSL vs the old Fortran agree to ~4e-14 relative on the raw fit
  outputs, and end-to-end prepfold topocentric timing parameters differ by ~1e-10 to 1e-9 of
  their quoted uncertainties (QR-vs-SVD roundoff, negligible for timing).

#### FFTFIT pure-Python rewrite (first Fortran-removal target)

FFTFIT (Taylor 1992) measures pulsar times-of-arrival by fitting a shift/scale between an
observed profile and a template, entirely in the Fourier domain. Current implementation:
`python/fftfit_src/*.f` (`fftfit`, `cprof`, `fccf`, `brent`, `ffft`) wrapped via f2py; the
public API used by `bin/get_TOAs.py` and `bin/sum_profiles.py` is just `cprof(template)` and
`fftfit(profile, amp, pha)`. Because this feeds high-precision timing, the replacement must be
validated very carefully against the original.

**Validation philosophy:** "match, but aim to improve" — bulk agreement with the Fortran must be
negligible relative to the reported `eshift`; cases where the port is *more* accurate or stable
(low-SNR, the `eshift=999` bailout, near-half-bin shifts) count as improvements, not failures.
Note the Fortran is single-precision (`real*4`) with truncated constants, so a double-precision
port will not be bit-identical by design.

**Oracles:** (1) a frozen golden dataset generated from the current Fortran [primary]; (2) Anne
Archibald's independent pure-Python FFTFIT [independent cross-check — see below]; (3) real-data
TOAs, kept lightweight by reusing the existing prepfold fixtures in `tests/prepfold/goodfolds/*.pfd`
rather than adding new data.

**Two implementations, selectable by method (decided):** ship both behind a method selector.
  * `"classic"` — the faithful NumPy/SciPy port of Taylor's Fortran. **Default**, for
    backward-compatible TOAs / reproducibility. Validated to ~0.28 sigma vs. the Fortran.
  * `"aarchiba"` — Anne Archibald's independent algorithm (cross-spectrum + upsampled `irfft` +
    bounded minimization + least-squares scale/uncertainty), adopted from her unmerged PINT
    branch. Opt-in improvement: handles the symmetric-profile / very-low-harmonic case the
    Fortran bails on, and different-length profiles.

**Adopting Anne's code (`aarchiba/PINT` branch `fftfit_new`, BSD-3-clause → GPL-compatible,
retain attribution):** her implementation appears abandoned upstream (she has left the field; the
PR, nanograv/PINT#777, was never merged). PRESTO can give it a maintained home. Her 844-line
`hypothesis` test suite is excellent and should be adopted/adapted. Note her API returns `shift`
in *turns* plus `scale`/`offset`/`uncertainty`; a thin adapter maps to PRESTO's
`(shift_bins, eshift, snr, esnr, b, errb, ngood)` tuple.

**Uncertainty calibration [decided: ship now, fix later].** Measured 1-sigma coverage (does the
error bar contain the truth ~68% of the time?) shows *both* methods under-cover by ~5-15%, worse
for broad profiles (e.g. kappa=1: aarchiba 0.656, classic 0.625 vs target 0.683; multi-sigma over
N=2000). This is a pre-existing FFTFIT trait inherited faithfully by the classic port, not a
regression; the shift itself is solid for both. Decision: proceed with the rewrite, document the
uncertainties as approximate/slightly optimistic, and treat calibration as a separate follow-up
(see below) rather than a blocker.

**Follow-up (tracked, non-blocking): FFTFIT uncertainty calibration.** Get 1-sigma coverage to
~0.683 across regimes. aarchiba is the better starting point (marginally better calibrated, and
its least-squares/covariance framework is amenable to a principled fix). Candidate fixes: resolve
Anne's `std`-normalization `FIXME`; investigate linearization bias for broad profiles; possibly an
empirical inflation factor validated against the coverage tests.

**Status: DONE** (except the non-blocking uncertainty-calibration follow-up above). The pure-Python
FFTFIT is the installed implementation (`python/presto/fftfit.py` + `fftfit_aarchiba.py`); the
Fortran and its f2py build are removed. `get_TOAs.py` gained `-A/--fftfit {classic,aarchiba}`
(default classic). Regression + `hypothesis` tests live in `python/fftfit_src/` and drive the
installed package; the frozen oracle (`fftfit_reference.npz`) is retained so tests survive the
Fortran's removal. End-to-end verified on two topocentric `goodfolds` folds: classic reproduces the
prior TOAs exactly; aarchiba agrees within ~0.1 sigma.

**Phased approach:**
0. [done] Golden-reference harness: freeze inputs+outputs of the current Fortran over a broad
   synthetic battery (`fftfit_reference.npz`), plus real-data cases from a `goodfolds` fold.
1. [done] Bottom-up "classic" port (`cprof` via NumPy `rfft`; `fftfit` with `scipy.optimize.brentq`
   for `brent.f` and NumPy for `ffft.f`, letting both Fortran files be dropped). Validated vs the
   frozen oracle: shift within 0.28 sigma, bailouts exact.
2. [done] Independent cross-check: vendored Anne's algorithm agrees with the port to <=0.23 bins on
   all well-determined fits (disagreements only at snr~1.5).
3. [done] Adopt Anne's code + test suite as the `"aarchiba"` method; build the method selector and
   the shift/uncertainty adapter. (Uncertainty calibration deferred as a tracked follow-up.)
4. [done] Assemble behind the public API (drop-in for `get_TOAs.py`, `sum_profiles.py`) with
   `"classic"` default and a `--fftfit` selector; validated against real-data TOAs from `goodfolds`.
5. [done] Swap the build: replace the f2py target with a pure-Python install, remove the Fortran,
   update the tests.

### Remove the TEMPO dependency

- **[DONE]** Barycentering: use **ERFA** instead of TEMPO. `barycenter()` in
  `src/barycenter.c` now computes UTC→TT→TDB (ERFA's built-in leap seconds + `eraDtdb`),
  the Roemer delay (`eraEpv00` analytic ephemeris + `eraPvtob`/`eraC2i06a` for the
  observatory, from a new internal ITRF table in `src/observatories.c`), and the solar
  Shapiro delay, all in-process with **no external files**. v/c comes from a central
  difference of the total delay, matching TEMPO's `femit/fobs - 1` definition (positive
  = moving away). Two conventions matter: PRESTO MJDs count elapsed-seconds/86400 from
  UTC midnight (not ERFA's "stretched" leap-day quasi-JD — up to 1 s different on a
  leap-second day), and UT1≈UTC is used (≤1.4 μs). Validated vs TEMPO/DE405 over 14
  observatories × 7 sky positions × 8 epochs (1995–2023, incl. leap-second days):
  absolute agreement ~14 μs max (the eraEpv00-vs-DE405 ephemeris difference), ~1.3 μs
  differential drift over a dense 12-hr series (mostly TEMPO's measured UT1 vs our
  UT1≈UTC), and |Δ(v/c)| ≲ 3e-9 (vs typical |v/c| ~1e-4). The old code is retained as
  `barycenter_tempo()` and `src/check_bary_erfa.c` re-runs the full comparison anytime
  (needs TEMPO); `tests/test_barycenter.py` checks against frozen TEMPO values
  (`tests/bary_tempo_reference.txt`) with no TEMPO needed. The `ephem` argument is now
  ignored.
- **[DONE]** Polyco creation: use **tempo2** (which is on conda-forge and sets `$TEMPO2`
  automatically) instead of TEMPO. `make_polycos()` in `src/polycos.c` and
  `polycos.create_polycos()` now run `tempo2 -tempo1 -f par -polyco "mjd1 mjd2 nspan
  ncoeff maxha sitename freq"`, which writes TEMPO1-format polycos that the existing
  `getpoly()`/`phcalc()` readers consume (one reader fix: the line-2 observatory field is
  now a site name like "gbt", so it is parsed as a string). Sites are passed as tempo2
  observatory *names* — tempo2's tempo1-code alias table is conflicting ('k' is e-Merlin,
  'y' is aliased to Jodrell), and the old TEMPO codes had two latent bugs anyway (ATA 's'
  was SHAO 65m; KAT-7 'k' was FAST). Validated vs TEMPO with `src/check_polycos_t2.c`
  (needs both programs): constant phase offsets 0.001–0.023 μs, differential drift
  ≤ 3e-6 turns, |Δf|/f ≤ 4e-11; all polyco-based prepfold regression folds pass against
  the *unmodified* goodfolds references and `get_TOAs.py` TOAs are identical. The old
  generator is retained as `make_polycos_tempo1()`. **With this plus ERFA barycentering,
  PRESTO no longer needs TEMPO at all.** (tempo2 caveat: no conda-forge osx-arm64 build.)

### Adopt ERFA/SOFA for astronomy/time routines

Replace hand-rolled code with ERFA/SOFA where possible: `mjd2cal`, `cal2jd`, and certain routines
in `misc_utils.c`. The barycenter↔topocenter calculations are done (see above), which also
provided the two prerequisites listed here: the observatory-position database now lives in
`src/observatories.c` (coordinates from TEMPO's obsys.dat, with PINT/tempo2 filling gaps), and
leap seconds come from ERFA's built-in table (`eraDat` — note this means very new leap seconds,
should they ever resume, require an ERFA update).

**[DONE]** The vendored Starlink SLALIB date/time routines (`src/cldj.c`, `src/djcl.c`,
`include/slalib.h`, `include/slamac.h`) are gone. `slaCldj` (Gregorian→MJD) is now `eraCal2jd`
via a new libpresto helper `cal_to_mjd()` in `misc_utils.c` (used by the PSRFITS/SPIGOT/BPP/WAPP
readers and the `cal2mjd` tool); `slaDjcl` (MJD→Gregorian) is now `eraJd2cal` (in `mjd_to_datestr()`
and the `mjd2cal` tool). While there, `sphere_ang_diff()` now wraps `eraSeps()` (agrees to ~4e-15
rad) and `dms2rad()`/`hms2rad()` wrap `eraAf2a()`/`eraTf2a()` so ERFA owns the sign conventions
(bit-identical output). No behavior change for any tool. Still to do: the Python `parfile.py`
optionally imports the old `slalib` package (`sla_ecleq`/`sla_eqecl`/`sla_eqgal`); those could
move to astropy (already a required dep), minding TEMPO's ecliptic-obliquity conventions.

### Replace numerical routines with GSL

Replace the in-tree median and statistics calls with well-tested GSL equivalents.

**[DONE]** `src/median.c` (a Numerical Recipes Quickselect) is now a thin wrapper around
GSL's `gsl_stats_float_median()`. Odd-length results are bit-identical to the old code;
even-length results now interpolate (mean of the two central values).

### Docstrings, typehints, and ruff

Add these things to at least the python routines python/presto, but especially psr_utils.py, 
prepfold.py, psrfits.py, spectra.py, filterbank.py, and sifting.py.

### Add Descriptions etc in .cli files for the man pages

Have Claude fill in Descriptions, Author, and SeeAlso fields in the pyCLIG .cli files to improve
the generated man pages. Read the code, `docs/PRESTO_search_tutorial.odp`, and the `FAQ.md` for 
relevant information. I would like the man files to start being more useful since there really 
isn't great documentation for PRESTO, in general. Include common usage examples for the main 
PRESTO routines (`rfifind`, `[mpi]prep*`, `accelsearch`, `realfft`, and `show_pfd`). Once complete,
re-generate all the CLIG-based files (except prepfold_multi, given the issues!). On a related note,
the relevant information in `prepfold_multi.1` should also go in a new section in the `FAQ.md`.

**[DONE]** Description/Author/SeeAlso added to all 15 remaining `.cli` files (usage examples
included for `rfifind`, `prepdata`, `prepsubband`, `mpiprepsubband`, `accelsearch`, `realfft`,
and `show_pfd`); all regenerated via pyclig except `prepfold_multi_cmd.cli`, which keeps its
existing hand-written content. `prepfold_multi`'s `-candfile` format is now also copied into a
new FAQ.md entry. Man page installation (next section) is still open.

### As an extension to the previous, figure out man page installation

Discuss with Claude ways to do this seamlessly and usefully, especially once we (soon!) 
get PRESTO into conda-forge.

**[DONE]** The `docs/*.1` man pages now install to `{prefix}/share/man/man1` via
`install_man()` in the top-level `meson.build`, so `man rfifind` (etc.) works with no
`$PRESTO` or extra `MANPATH` setup (conda activation already adds `{prefix}/share/man`).

### Clean up old and unused code

Ask Claude to examine files in `src`, `include`, `tests`, and in `python` to find things that
can safely be removed. These should be things that have not been significantly touched in a
a long while, are irrelevant because of other code. I would like to approve a list item by 
item for those things to be removed.

**[DONE]** `dftfold` (outdated, unused, no references elsewhere in the tree) was removed
while doing the `.cli` Description/Author/SeeAlso pass below.

**[DONE]** `a2x` (the old ASCII→PostScript helper) and its files (`bin/a2x.sh`, `lib/a2x.ps`,
`lib/a2xshell.ps`) were removed, along with the dead code that shelled out to it:
`src/cand_output.c` (both `file_reg_candidates`/`file_bin_candidates` were unused — the former
only by the already-broken legacy `tests/test_apps.c`) and its `presto.h` prototypes, plus the
commented-out a2x block in `minifft.c`. The dead `tests/test_apps.c` (missing `fftapps.h`/
`plotting.h`, unbuildable) and its `test_fdot`/`test_plot`/`test_finterp` `.mak` files were
removed too.

### Replace code of questionable origin or which can be replaced with GSL or ERFA (or scipy)

- One or two of the .c files (e.g. `amoeba.c` and `median.c`, plus parts of `misc_utils.c`) 
  originate in Numerical Recipes. It would be great to replace those.
  **[DONE for `amoeba.c`]** replaced by a GSL simplex in `maximize_rz` (commit 8f41787).
  **[DONE for `median.c`]** now wraps GSL's `gsl_stats_float_median()` (see above).
- **[DONE]** `solvopt.c` and `apprgrdn.c` were dead — no callers anywhere and not even in
  the meson build (never compiled) — so they were simply removed.
- **[DONE]** Replaced the randlib stuff (`src/randlib.c`, `src/com.c`, `include/randlib.h`)
  with GSL's `gsl_rng`/`gsl_ran_gaussian`/`gsl_ran_poisson`. It was only used by `makedata`
  (and the unbuilt legacy `tests/test_ffts.c`); those files are now gone.
- **[DONE]** Replaced the vendored SLALIB time/coordinate utilities (`cldj.c`, `djcl.c`, and
  the `dms2rad`/`hms2rad`/`sphere_ang_diff`/`mjd_to_datestr` routines in `misc_utils.c`, plus
  the `mjd2cal` and `cal2mjd` tools) with ERFA calls. See the ERFA section above for details.

### Develop a plan, with scripts, for doing proper tagged releases

This will be important for simplfying updates when PRESTO is on conda-forge

**[DONE]** `RELEASE.md` documents the repeatable checklist (version stamp &rarr; CHANGELOG/README
&rarr; tag &rarr; push &rarr; `gh release` &rarr; hash &rarr; conda-forge). `determine_version.py`
gained a `--set X.Y.Z` mode to stamp an explicit release version (overriding the git-derived dev
version) into the three files it keeps in sync.

### Develop scripts and recipes for getting PRESTO onto conda-forge

When that is complete (and integration into conda-forge is imminent), we will tag v6.

**[IN PROGRESS]** v6.0.0 is stamped/tagged, and a conda-forge recipe is kept in-repo under
`conda-recipe/` (`recipe.yaml` + `build.sh`). It builds both halves (meson for `libpresto`/tools,
then `pip` for the Python package), vendoring the ERFA release tarball into
`subprojects/erfa-2.0.1/` since ERFA has no conda-forge feedstock. First submission targets
`linux-64` + `osx-64`; `osx-arm64` waits on a conda-forge `tempo2`. Remaining: get the
`conda-forge/staged-recipes` PR reviewed and merged (creates `presto-pulsar-feedstock`; the
plain `presto` name is taken on bioconda by an unrelated tool), after which
the autotick bot handles version bumps.

## Nice-to-haves

### Revitalize CI

The old Travis CI setup (`.travis/`, `.travis.yml`) was stale for years and has been removed
(2026-07). Set up modern CI — presumably GitHub Actions — that builds the meson project and
Python package and runs the test suite: `tests/test_barycenter.py`, `tests/test_presto_python.py`
(now TEMPO-free), the fftfit tests in `python/fftfit_src/`, and ideally the prepfold suite
(`tests/prepfold/prepfold_tests.sh`, which downloads its ~100 MB fixture and needs TEMPO for
the polyco-based folds). A linux + macOS matrix would also catch the recurring macOS
build/rpath issues early.

## Other known issues to fix

### `python/presto/filterbank.py` low-bit-depth support

The pure-Python `filterbank.py` reader only supports 8-, 16-, and 32-bit filterbank files
(`check_nbits` rejects everything else). It should be extended to also read **1-, 2-, and
4-bit** data, which the SIGPROC format and PRESTO's C reader (`readfile`, `prep*`) all handle
(e.g. the 4-bit `GBT_Lband_PSR.fil`). This needs sub-byte unpacking in `get_spectra` (samples
are packed multiple-per-byte) plus `get_dtype`/`check_nbits` updates; the framing math
(`bytes_per_spectrum = nchans*nbits // 8`) is already correct for packed low-bit data. Verify
against `readfile` output on a known file.

### `prepfold_multi` generated command-line files

`prepfold_multi` currently requires **manual** edits to its auto-generated
`src/prepfold_multi_cmd.[ch]` after they are generated from `clig/prepfold_multi_cmd.cli` — the
generator (`pyclig`) output is not yet correct/complete on its own. Fix so `prepfold_multi`
regenerates cleanly from its `.cli` like every other tool. Until then, do not blindly overwrite
the hand-edits when regenerating.

## Postponed, potentially indefinitely

### Drop support for retired pulsar backends

(Note: Postponed since all the code now compiles without warnings and will likely not be touched
again. And we occasionally find old files from those machines where it is (very!) useful to be 
able to see what they are via `readfile`.)
Remove all code for old backends (Spigot, BCPM, WAPP, and others). This is believed to affect 
only the `readfile` code path.
