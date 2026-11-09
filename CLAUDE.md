# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What PRESTO is

PRESTO (PulsaR Exploration and Search TOolkit) is a suite of pulsar search and analysis
software. The performance-critical core is ANSI C (Fortran-free since v6.0.0); higher-level
analysis, plotting, and orchestration are Python. The C code builds into a shared library
(`libpresto`) plus a set of standalone command-line executables. The Python package
(`presto`) wraps the same C routines via a SWIG extension and adds many pure-Python modules.

Raw-data support: PSRFITS search-format, SIGPROC filterbank (1/2/4/8/32-bit), single-precision
float time series (with a `.inf` header), and event/photon arrival times.

Planned larger changes (currently: CI, then a C port of CoherentSearch.jl with Highway SIMD
kernels, developed on a feature branch) live in [ROADMAP.md](ROADMAP.md).

## Build and install

The build uses **meson** (v5+ switched away from Makefiles). The two halves — C code and the
Python package — are built and installed separately. The `PRESTO` environment variable is now
**optional** (installed tools find their data in `{prefix}/share/presto`; setting `PRESTO` to
this top-level directory still works as an override for running from a source tree). Set
`PGPLOT_DIR` (and `TEMPO2` if generating polycos;
conda-forge's tempo2 sets it automatically). Note: `$PRESTO/bin` and
`$PRESTO/lib` should **not** be on `PATH`/`LD_LIBRARY_PATH`/`PYTHONPATH` (unlike older versions).

```sh
# 1. C code and libpresto
meson setup build --prefix=$CONDA_PREFIX   # or $VIRTUAL_ENV, or $HOME, or omit for /usr/local
python check_meson_build.py                # sanity-checks env vars vs. the meson config
meson compile -C build
meson install -C build

# 2. Python package and C extension bindings
cd python && pip install .
```

`check_meson_build.py` is the first thing to run when a build/link problem appears; it reports
the effective `{prefix}`, `{bindir}`, `{libdir}` and flags common environment mistakes. Build
logs are in `build/meson-logs/`. Many compiler warnings are expected and benign.

Alternatively, `pixi run build` drives both stages inside a self-contained pixi environment
(toolchain + deps from conda-forge; see `pixi.toml`). It uses its own `build-pixi/` directory
and installs into `.pixi/envs/default`, so it never touches a manually managed `build/`.

Uninstall: `ninja uninstall -C build` and `pip uninstall presto`.

`makewisdom` (installed) generates `fftw_wisdom.txt` for slightly faster FFTs, writing it
straight to where the tools look for it at runtime — `{prefix}/share/presto` (or `$PRESTO/lib`
if `PRESTO` is set). Pass an explicit path argument to override the destination.

## Tests

`pixi run test` builds everything and runs the known-good suite, and GitHub Actions
(`.github/workflows/ci.yml`) runs that same task on every push and PR, on linux-64,
linux-aarch64 and osx-arm64, plus a linux-64 build with `-funsigned-char` (plain `char` is
unsigned on aarch64), which catches signed-char assumptions on x86. The test list is **explicit** in the `test` task in
`pixi.toml`: `tests/test_barycenter.py`, `tests/test_binary_calcs.py`,
`tests/test_maximize_rz.py`, `tests/test_pipeline_synthetic.py`, `tests/test_sifting.py`, the
fftfit suite in
`python/fftfit_src/` (needs `hypothesis`), and `tests/test_presto_python.py`. Add a new test
to that list when you add it.

`tests/test_pipeline_synthetic.py` is the end-to-end test: `tests/synthpulsar.py` writes a
deterministic 64-channel filterbank with a known injected pulsar (23.4567 Hz, 10% duty,
S/N 12, DM 50). The pipeline then runs on it (prepdata/prepsubband → realfft → rednoise →
accelsearch → prepfold on both the `.dat` and the raw file), and each result is checked
against the truth. Thresholds sit below the measured values, which are recorded at each
assert; re-measure before tightening them. (`sifting` replaces the ACCEL file's sigma with a
single-trial one, so `cand.sigma` there is higher than the ACCEL file's column.) The same
file round-trips 8-bit, signed 8-bit, 16-bit and float filterbanks through `prepdata -dm 0`
sample by sample, and checks 1/2/4-bit files channel by channel against both the Python and C
readers (via `prepsubband -sub`). Use `synthpulsar.py`
for any new test that needs raw data.

**Real-data tests** run separately: `pixi run test-realdata` (~7 min; it is just
`test-realdata-search` plus `test-realdata-prepfold`, which CI runs as separate steps so one
failing does not hide the other), which GitHub Actions
(`.github/workflows/realdata.yml`) runs weekly (Mondays) and on demand (`workflow_dispatch`),
not on every push. They use the 173 MB Terzan 5 GBT PSRFITS file
`Ter5_080912_short2bits.fits` (2-bit, 64 channels, 859 s). It is **not in git**: it is an
asset on the repo's `test-data-v1` GitHub release (a pre-release, kept off "Latest"), and
`tests/prepfold/get_test_data.sh` downloads it into `tests/prepfold/` (gitignored) and checks
its SHA-256. The tag is versioned on purpose: to change or add data, upload to a new tag
(`test-data-v2`, ...) and update the script's tag and checksum, never replace an asset in
place, so older commits still fetch the data their references were made from. The task runs:
- `tests/test_realdata_search.py` (pytest; skips if the file is absent): rfifind →
  prepsubband (DMs 234–248, barycentered, masked) → realfft → rednoise → `accelsearch -zmax 50`
  → sifting → prepfold, plus `single_pulse_search.py`. The truth is the par files in
  `tests/prepfold/`: `binary_psr` predicts the Doppler-shifted frequency and z of Ter5A
  (bright, z ≈ −5) and Ter5N (weak, DM 238), and Ter5A's single pulses at DM 242 must be
  phase-locked to its ephemeris (Rayleigh test). As for the synthetic test, thresholds sit
  below measured values recorded at each assert.
- the `tests/prepfold/` shell suites `prepfold_tests.sh` (`.bestprof` regressions against
  `goodfolds/`, `reader_equiv_check`, the cmd-layout check) and `prepfold_multi_tests.sh`
  (multi == many singles). See `tests/prepfold/README.md`.

**`tests/test_presto_python.py` is one of the most important tests**: it exercises the SWIG
wrappings of crucial C routines (FFTs, folding, etc.). It is a plain script with top-level
asserts, not pytest functions; a failure exits 1.

Many other Python scripts in `tests/` (e.g. `dedisp.py`, `test_max_rzw_harmonics.py`,
`tests/python/`) and `examplescripts/` are old demos or stale smoke tests that no longer
work as tests, so CI deliberately does not glob-collect `tests/`. The legacy C test programs
with `.mak` makefiles are not built by meson. `examplescripts/ffdot_example.py` is still a
useful manual check, but it waits for input at a PGPLOT prompt.

## Architecture

### C library and executables (`src/`, `include/`, `lib/`)

`src/meson.build` is the authoritative map of what gets built. It defines:
- `libpresto` — the shared library: FFT/correlation routines, dedispersion, folding,
  barycentering, response functions, candidate handling, FITS/filterbank I/O, etc.
- Standalone executables (e.g. `accelsearch`, `prepfold`, `prepdata`, `prepsubband`,
  `rfifind`, `zapbirds`, `readfile`), most linking against `libpresto`.
- Grouped object sets reused across executables: `INSTRUMENTOBJS` (backend/format readers),
  `PLOT2DOBJS` (PGPLOT helpers).
- `install: false` on utility/experimental executables that are built but not installed.

There is also a **legacy `src/Makefile`** (pre-meson). It is *not* part of the meson/pixi
build, but the maintainer still uses it for quick standalone compilation checks during
development, so **keep it roughly in sync** when you add/remove/rename C source files: mirror
the change in its `PRESTOOBJS` list (which must match `libpresto`'s source list in
`src/meson.build`) and, for a new/removed program, its `BINARIES` list and build rules. It
won't always be perfectly current, but don't let it drift on file additions/removals you make.

Key external dependencies: **FFTW3 single-precision** (`fftw3f` — the `f` matters), GSL, **ERFA**
(barycentering), glib-2.0, cfitsio, PGPLOT/cpgplot, X11, libpng, and optional OpenMP and MPI
(MPI gates `mpiprepsubband`). The **tempo2** executable is called at runtime for polyco
generation (`prepfold -timing`); TEMPO is not used at all anymore.

### Command-line parsing is generated (`clig/` → `src/*_cmd.c`)

Executables that take many options do **not** hand-write their argument parsing. Each has a
`clig/<name>_cmd.cli` spec file that generates `src/<name>_cmd.c` (and the matching header). The
generated files are committed to the repo. **To change a program's command-line options, edit
the `.cli` file and regenerate, do not edit `_cmd.c` by hand.** `.cli` syntax: typed option
declarations (`Int`, `Float`, `Flag`, `String`) with ranges (`-r min max`) and defaults (`-d`).

The generator is **`pyclig`** (in `~/bin`, not vendored — source is on the maintainer's GitHub).
It is brand-new and only used locally for now; it replaces and updates the old, unmaintained
`clig` program.

**Known wart / TODO:** `prepfold_multi` is a brand-new feature whose auto-generated
`prepfold_multi_cmd.[ch]` files currently require *manual* modification after generation — the
generated output is not yet correct/complete on its own. This should be fixed so `prepfold_multi`
regenerates cleanly from its `.cli` like every other tool. Until then, be careful not to blindly
overwrite hand-edits when regenerating `prepfold_multi_cmd.[ch]`.

### Python package (`python/`)

`python/meson.build` orchestrates four subdirectories, each with its own `meson.build`:
- `presto_src/` — builds the `_presto` SWIG C extension module against `libpresto`; exposes it
  as `presto.presto` (import via `from presto import presto`).
- `presto/` — the pure-Python package: analysis and I/O modules such as `infodata`, `sigproc`,
  `psrfits`, `prepfold`, `polycos`, `parfile`, `residuals`, `sifting`, `psr_utils`,
  `binary_psr`, `observatories` (telescope name to observatory code), and plotting via
  `Pgplot`. Data files like `cosine_rand.json` are installed alongside the code.
- `fftfit_src/` — FFTFIT, now pure Python (NumPy/SciPy); holds its tests and a reference
  dataset. The old Fortran/f2py implementation was removed in v6.0.0.
- `ppgplot_src/` — Python PGPLOT bindings.

`wrappers/` contains older SWIG interface files (`.i`) and is largely legacy relative to the
active `presto_src/` binding.

Python entry-point scripts live in `bin/` (e.g. `DDplan.py`, `single_pulse_search.py`,
`get_TOAs.py`, `waterfaller.py`, `simple_zapbirds.py`); `bin/meson.build` installs them.

## Conventions

- Version numbers live in several places and must be kept in sync when releasing: top-level
  `meson.build`, `python/pyproject.toml`, and `python/presto/__init__.py`. See
  `determine_version.py` and `CHANGELOG.md`.
- **Add a `CHANGELOG.md` entry for any non-trivial change** (new/changed/removed features,
  dependency changes, behavior changes, notable cleanups — anything a user or packager might
  care about). Add it as a bullet at the top of the `## Development (unreleased...)` section
  (newest first), in the same style as the existing entries, as part of the same commit. Trivial
  changes (typo fixes, comment tweaks, minor internal refactors) don't need one.
- Most compiled programs print usage when run with no arguments — the fastest way to learn a
  tool's interface (and to confirm a `.cli` edit worked after rebuild).
- macOS builds need workarounds (linker flags, `fix_rpath_duplicates.py` install script); see
  the macOS section of `INSTALL.md`.
- **Git workflow:** Scott is effectively the sole developer, so routine changes — bug fixes,
  refactors, code replacements (e.g. the amoeba→GSL swap), doc updates, etc. — may be committed
  **directly to `master`** without a branch or PR (this overrides the usual "branch off the
  default branch first" default). Reserve feature branches and PRs for very large-scale changes,
  or for work coming from other contributors. Still only commit when asked.
