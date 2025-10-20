# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About PRESTO

PRESTO (PulsaR Exploration and Search TOolkit) is a large suite of pulsar search and analysis software written primarily in ANSI C with Python bindings and utilities. It's designed to efficiently search for binary millisecond pulsars and has discovered over 1500 pulsars. The codebase is ~80k lines of C code plus extensive Python modules.

## Build System

PRESTO uses **meson** (not make) as the primary build system as of v5.0:

### Initial Setup
```bash
# Set PRESTO environment variable to repository root
export PRESTO=/path/to/presto

# Install build dependencies
pip install meson meson-python ninja
pip install --upgrade pip

# Configure meson (choose one prefix option)
meson setup build --prefix=$CONDA_PREFIX        # For Conda
meson setup build --prefix=$VIRTUAL_ENV         # For virtualenv
meson setup build --prefix=$HOME                # For home directory
meson setup build                               # For /usr/local (requires sudo)

# Verify configuration before building
python check_meson_build.py
```

### Building and Installing
```bash
# Build C/Fortran code and libpresto shared library
meson compile -C build
meson install -C build

# Install Python modules
cd $PRESTO/python
pip install .
```

### Testing
```bash
# Run basic tests
python tests/test_presto_python.py
python examplescripts/ffdot_example.py
python python/fftfit_src/test_fftfit.py
```

### Cleaning
```bash
# Uninstall everything
cd $PRESTO/build
ninja uninstall
pip uninstall presto

# Clean old builds (if needed)
cd $PRESTO/src
make cleaner
```

### Legacy Makefile
The `src/Makefile` still exists but is **deprecated**. It was used in PRESTO v4 and earlier. Do not use it for building; use meson instead.

## Architecture Overview

### Core Components

**C Library (`src/` → `lib/libpresto.so` or `.dylib`)**
- Core signal processing and pulsar search algorithms (~80k lines)
- Key modules:
  - FFT operations (`fftcalls.c`, `fastffts.c`, `realfft.c`)
  - Acceleration/jerk search (`accelsearch.c`, `accel_utils.c`)
  - Dispersion handling (`dispersion.c`)
  - Binary orbit integration (`orbint.c`)
  - Data I/O for various formats (`backend_common.c`, `psrfits.c`, `sigproc_fb.c`)
  - Folding and period search (`fold.c`, `prepfold.c`)
  - RFI removal (`rfifind.c`, `zapping.c`)
- Headers in `include/` with `presto.h` as the main interface

**Python Bindings (`python/presto_src/`)**
- SWIG-generated bindings (`presto_wrap.c`, `prestoswig.py`) expose C functions to Python
- Built via meson-python (see `python/meson.build` and `python/pyproject.toml`)

**Python Modules (`python/presto/`)**
- High-level pulsar analysis tools
- Key modules:
  - `psrfits.py`: PSRFITS data handling
  - `infodata.py`: Metadata for time series
  - `binary_psr.py`: Binary pulsar calculations
  - `parfile.py`: Parse pulsar parameter files
  - `fftfit.py`: FFT-based pulse profile fitting
  - `singlepulse/`: Single pulse search tools (rrattrap, spcand)

**Executables (`src/` → `bin/`)**
- Command-line tools for pulsar searching and analysis
- Major programs:
  - `prepdata`, `prepsubband`: Data preparation and de-dispersion
  - `accelsearch`: Fourier-domain acceleration/jerk search
  - `prepfold`: Fold candidates to create diagnostic plots (.pfd files)
  - `rfifind`: RFI detection and removal
  - `single_pulse_search.py`: Find single pulses
  - Data exploration: `readfile`, `exploredat`, `explorefft`

**Python Scripts (`bin/`)**
- High-level analysis scripts (`.py` files)
- Important utilities:
  - `DDplan.py`: Plan de-dispersion strategy
  - `fourier_fold.py`: Fold using FFT amplitudes
  - `get_TOAs.py`: Generate timing TOAs
  - `sum_profiles.py`: Stack pulse profiles
  - `stacksearch.py`: Stack search multiple FFTs
  - Single pulse tools: `rrattrap.py`, `make_spd.py`, `waterfaller.py`

### Data Flow Architecture

1. **Raw Data Input** → Various formats supported via `backend_common.c`:
   - PSRFITS search format (modern standard)
   - SIGPROC filterbank (1-bit to 32-bit float)
   - Time series with `.inf` metadata files
   - Photon events (ASCII or binary)

2. **Preparation Pipeline**:
   - `rfifind` → detects RFI, creates `.mask` files
   - `prepdata`/`prepsubband`/`mpiprepsubband` → de-disperse, create `.dat` and `.inf` files
   - `realfft` → FFT time series to `.fft` files

3. **Search Pipeline**:
   - `accelsearch` → search `.fft` for accelerated signals → `.ACCEL` candidates
   - `single_pulse_search.py` → search `.dat` for single pulses
   - `prepfold` → fold candidates → `.pfd` diagnostic files
   - `show_pfd` → view/manipulate `.pfd` files

4. **Analysis/Timing**:
   - `get_TOAs.py` → extract TOAs from `.pfd` files
   - TEMPO integration for timing solutions

### File Format Conventions

- `.inf`: Text metadata for time series (sample rate, telescope, etc.)
- `.dat`: Time series data (float32 by default)
- `.fft`: FFT output (complex64 amplitudes)
- `.pfd`: Prepfold diagnostic data (binary, contains folded profiles)
- `.bestprof`: Text summary of best fold parameters from `prepfold`
- `.ACCEL`: Acceleration search candidates (text)
- `.mask`: RFI mask from `rfifind`

## Dependencies

**Required external libraries:**
- FFTW3 (single-precision): FFT operations
- PGPLOT: Plotting (especially for `.pfd` diagnostic plots)
- CFITSIO: FITS file I/O
- GLIB2: Data structures and utilities
- TEMPO: Pulsar timing (separate installation, set `$TEMPO` env var)

**MacOS-specific notes:**
- Use MacPorts or Homebrew for dependencies
- May need to set `LIBRARY_PATH`, `CC=gcc`, and `CFLAGS` environment variables
- See INSTALL.md "For MacOS Users" section for detailed instructions
- Fortran compiler (gfortran) required

## Development Notes

### C Code Style
- Mostly K&R C style (older code is being converted to ANSI C - see recent commits)
- Heavy use of custom vector operations (`vectors.c`)
- FFTW malloc wrappers for aligned memory
- OpenMP support available (set `USEOPENMP` in old Makefile, configured via meson)

### Modifying C Code
- Headers in `include/`, implementations in `src/`
- After changes: `meson compile -C build && meson install -C build`
- If modifying Python bindings: rebuild and reinstall Python package

### CLIG Command-Line Interface
- Many programs use CLIG for argument parsing (`clig/` directory)
- `.cli` files define command-line interfaces
- Generate `*_cmd.c` and `*_cmd.h` files (already generated, rarely need to regenerate)

### Python Development
- Python code follows typical scientific Python conventions
- Dependencies: numpy, scipy, astropy, matplotlib
- After modifying Python code: `cd $PRESTO/python && pip install .`

### Common Issues
- **Shared library errors**: Add libpresto location to `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (macOS)
- **Compile-time linking errors**: Set `LIBRARY_PATH` to include libpresto location
- **Memory errors in tests**: Try rebuilding without `-DUSE_FFTW_MALLOC` flag
- Environment variables: Use `python check_meson_build.py` to diagnose setup issues

## Working with Git

- Main branch: `master`
- Current branch: Uses descriptive names (e.g., `bugfix/ansi`)
- Follow standard git workflow for contributions
- Classic branch (`origin/classic`) contains old, unsupported code for legacy data formats
