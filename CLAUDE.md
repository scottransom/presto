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

Planned larger changes (conda-forge packaging, dependency removal, etc.) live in
[ROADMAP.md](ROADMAP.md).

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

There is no single unified test runner. The main smoke tests are:

```sh
python tests/test_presto_python.py        # exercises the SWIG bindings (FFTs, folding, etc.)
python examplescripts/ffdot_example.py
python python/fftfit_src/test_fftfit.py
```

`tests/` also holds pytest-style Python tests (`test_*.py`, e.g. `test_binary_calcs.py`,
`test_max_rzw_harmonics.py`) and legacy C test programs with `.mak` makefiles. There is
currently no CI (the stale Travis setup was removed; see ROADMAP.md).

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
`get_TOAs.py`, `waterfaller.py`, `rrattrap.py`); `bin/meson.build` installs them.

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
