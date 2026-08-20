# New and improved steps to install PRESTO!

## Easiest route: build with [pixi](https://pixi.sh)

If you use pixi, the repo now has a `pixi.toml` workspace that provides the
entire toolchain and all build dependencies from conda-forge and drives the
full build:

    pixi run build     # C code + libpresto + Python package, installed into the pixi env
    pixi run test      # quick smoke tests
    pixi shell         # drop into the environment (prepfold, python, etc. on PATH)

Use `pixi run install-py-editable` after building if you want live-editable
pure-Python modules while debugging. The pixi tasks use their own meson build
directory (`build-pixi/`), so they will not interfere with a manually managed
`build/`. Notes: ERFA is built automatically (it is not on conda-forge), and
tempo2 (used for polyco generation) is included from conda-forge — no TEMPO
installation is needed at all.

## td;dr for experienced PRESTO builders:
With v5, we have switched to building and installing with [meson](https://mesonbuild.com/).

**MacOS users should see the comments at the bottom of this file!**

**If you are interested in using Docker or Singularity containers of PRESTO, see the bottom of this file!**

As always, there are a set of essential packages required to build PRESTO. This command should do it on a Debian/Ubuntu-like system:
`apt install git build-essential libfftw3-bin libfftw3-dev libgsl28 libgsl-dev liberfa-dev pgplot5 libglib2.0-dev libcfitsio-bin libcfitsio-dev libpng-dev ghostscript gfortran tcsh autoconf libx11-dev python3-dev python3-numpy python3-pip`

Setting the `PRESTO` environment variable is now **optional** — the installed tools locate their runtime data (the pulsar catalog, etc.) automatically in `{prefix}/share/presto`. It is still handy as a shorthand for the commands below (which use `$PRESTO`) and works as an override when running from an uninstalled source tree, so you may want to point it at the top-level PRESTO git checkout. Either way, make sure that `$PRESTO/lib` and `$PRESTO/bin` are **not** in your `PATH` or `LD_LIBRARY_PATH` or `PYTHONPATH` environment variables as we have required in the past. It is probably a good idea to clean your earlier compiles, as well.

From your activated Python virtual or [Conda](https://docs.conda.io/) environment, make sure that you have `meson`, `meson-python`, `ninja`, and a recent `pip` installed (also Python >=3.8):

`pip install meson meson-python ninja` or `conda install meson meson-python ninja`

and

`pip install --upgrade pip`.

Now setup the C/Fortran code builds:

`cd $PRESTO`

`meson setup build --prefix=$CONDA_PREFIX` if you use Conda/Mamba/Anaconda

or

`meson setup build --prefix=$VIRTUAL_ENV` if you use a Python virtual environment

or

`meson setup build --prefix=$HOME` if you want things installed in your `$HOME/bin`, `$HOME/lib` directories.

or, if you have `sudo` permissions and want to install to `/usr/local` (or equivalent system directory):

`meson setup build`.

Note that if you don't want to set `--prefix` on the command line, you can also edit the top-level `meson.build` file and add e.g. `'prefix=/home/sransom'` to the `default_options` variable in the project definition.

In order to avoid the possible issues with linking and running, I recommend doing:

`python check_meson_build.py`.

If all looks good, it will let you know. If not, I recommend trying to fix the issues that it caught, and then starting again.

Now do the actual build and install via:

    meson compile -C build
    meson install -C build

You will see a lot of compiler warnings due to my crappy C-coding (I should fix those...), but that's OK. There should be logs in case anything goes wrong in `$PRESTO/build/meson-logs`. You should be able to run `prepfold` (for example), at this time and see the usage information. If that doesn't work, see the **troubleshooting** information below.

Then finally, install the Python codes and bindings via `pip`:

    cd $PRESTO/python
    pip install .

And that should do it! You can quickly test to see if most things are working by doing:

    cd $PRESTO
    python tests/test_presto_python.py
    python examplescripts/ffdot_example.py
    python python/fftfit_src/test_fftfit.py

Another good test is to see if you can run and fit the default profile in `pygaussfit.py`

If you want slightly faster FFT calls, just run the installed `makewisdom`. By default it uses `FFTW_MEASURE` and skips the largest FFT sizes, so it finishes in a few seconds; add `-patient` to instead use `FFTW_PATIENT` over the full set of sizes for higher-quality (but much slower, ~10-20 min) wisdom. It writes `fftw_wisdom.txt` straight into `{prefix}/share/presto` (or `$PRESTO/lib` if you have `PRESTO` set), which is where PRESTO looks for it at runtime — no manual copying needed. You can pass an explicit output path as an argument to override the destination.

**Optimization flags:** the default `buildtype=debugoptimized` (`-O2` plus `-g`) is a good
choice, and switching to `buildtype=release` (`-O3`, no `-g`) is *not* worth it — it makes no
measurable difference to run times, and the debug symbols cost nothing at runtime while making
any crash report far more useful.

What *does* help is letting the compiler use your CPU's vector instructions:

    meson setup build --prefix=$CONDA_PREFIX -Dc_args=-march=native

On a modern x86-64 machine (i.e. one with AVX2 and FMA) that makes `accelsearch` roughly
1.5-2x faster. Three caveats:
  * **Only use this for builds you will run on the machine that compiled them.** A binary built
    with `-march=native` will die with an illegal-instruction error on an older CPU, so don't do
    this for a shared cluster install, a container, or a package you distribute.
  * It only helps the compute-bound tools. I/O-bound ones like `prepsubband` see no change at
    all (their output is bit-identical).
  * It changes results at the roundoff level, since the compiler may then contract multiplies
    and adds into FMA instructions. In an `accelsearch` test this meant identical candidates
    with identical sigmas, and powers/frequencies differing in the last digit or two — well
    inside the quoted uncertainties — but a `-march` build will not reproduce an earlier run
    bit-for-bit.

Note that you can uninstall everything via:

    cd $PRESTO/build
    ninja uninstall
    pip uninstall presto

**Some quick troubleshooting info for the new meson build process:**
  * If you get an error when trying to run one of the compiled executables (like `prepfold`) about shared libraries (e.g. `prepfold: error while loading shared libraries: libpresto.so: cannot open shared object file: No such file or directory`) you will need to add the install location of `libpresto.so` to your `LD_LIBRARY_PATH` environment variable.
    * Note: You can find the install location by doing: `grep libpresto $PRESTO/build/meson-logs/install-log.txt`
  * If you get an error when running `pip install ...` for the Python codes saying `../meson.build:29:15: ERROR: C shared or static library 'presto' not found`, you will need to do something similar to the above and add the install location of `libpresto.so` to `LIBRARY_PATH` (note that that is **not** `LD_LIBRARY_PATH`! The former is for runtime linking and the latter is for linking at compile time. Note that you might need to use `DYLD_LIBRARY_PATH` on macOS!)

## Detailed install and build instructions:

**(Note:  For MacOS users, please see the bottom of the document!)**

1.  **Install [FFTW3](http://www.fftw.org)**

    I highly recommend that you use pre-compiled packages for your OS/distribution! For example, Ubuntu has good FFTW packages: `libfftw3-bin` and `libfftw3-dev`.

    If you compile your own, you need to compile FFTW for **single** precision. For all architectures I recommend the following configuration: `./configure --enable-shared --enable-single`

    If you are on a modern Intel processor and have a recent version of GCC, you can get much improved performance by adding: `--enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-fma`

2.  **Install [PGPLOT](http://www.astro.caltech.edu/~tjp/pgplot/)**

    I highly recommend that you use pre-compiled packages for your OS/distribution! For example, Ubuntu has a good PGPLOT package: `pgplot5`. You will likely need to set the `PGPLOT_DIR` environment variable. On Ubuntu, I have `PGPLOT_DIR=/usr/lib/pgplot5`

    If you want to try to compile your own (good luck!), you need the X-windows and postscript drivers at a minimum.

3.  **Install [tempo2](https://bitbucket.org/psrsoft/tempo2/)** (only needed for polycos)

    tempo2 is used for polyco generation (e.g. `prepfold -timing`); if you don't fold with polycos you can skip it. The easiest install is from conda-forge (`conda`/`pixi install tempo2`), which also sets the required `TEMPO2` runtime environment variable automatically (note: there is no conda-forge build for osx-arm64). If you build from source instead, set `TEMPO2` to the runtime directory yourself.

    Note: as of the post-v5.3.1 development code, **TEMPO is no longer used at all** — barycentering is done in-process via ERFA, and polycos come from tempo2.

4.  **Install [ERFA](https://github.com/liberfa/erfa)** (optional -- can be automatic!)

    ERFA does the barycentering-related astronomical calculations. If it is not already installed, **meson will automatically download and build it** as a subproject (via `subprojects/erfa.wrap`), so you can skip this step if your build machine has network access. To use a system copy instead: on Ubuntu the packages are `liberfa1` and `liberfa-dev`. Note that the `pyerfa` Python package (conda-forge/PyPI/astropy) does *not* provide the C library and headers that PRESTO needs. For a fully offline build without a system ERFA, download the [release tarball](https://github.com/liberfa/erfa/releases) and place it in `subprojects/packagecache/` (keeping the original filename) before running `meson setup`.

5.  **Install [GLIBv2](http://library.gnome.org/devel/glib/)**

    On Linux machines this is almost certainly already on your system (check in `/usr/lib` and `/usr/include/glib*`). Although you may need to install a glib development package in order to have the required include files. On Ubuntu, the package you need is: `libglib2.0-dev`

6.  **Install [GSL](https://www.gnu.org/software/gsl/)**

    I highly recommend that you use pre-compiled packages for your OS/distribution! For example, Ubuntu has good GSL packages: `libgsl28` and `libgsl-dev`.

7.  **Install [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/)**

    I highly recommend using pre-compiled packages, once again (on Ubuntu they are `libcfitsio-bin` and `libcfitsio-dev`), however, this is a very easy install via source.

8.  **(Optional) Set the `PRESTO` environment variable**

    This is no longer required — the installed tools find their runtime data in `{prefix}/share/presto` on their own. Setting `PRESTO` to the top level directory of the PRESTO distribution (i.e. this directory) is still useful as a shorthand for the commands below and as an override for running from an uninstalled source tree. Whether or not you set it, make sure that `$PRESTO/lib` and `$PRESTO/bin` are **not** in your `PATH` or `LD_LIBRARY_PATH` or `PYTHONPATH` environment variables as we have required in the past.

9.  **Activate your Python virtual environment *or* Conda/Mamba/Anaconda environment**

    * That environment should have `numpy` installed at a minimum
    * Make sure that `pip` is recent (`pip install --upgrade pip`)
    * Install the build tools: `pip install meson meson-python ninja` or `conda install meson meson-python ninja`

10.  **Configure the meson build**

    In the top level PRESTO directory (i.e. `cd $PRESTO`), configure `meson` via:

    `meson setup build --prefix=$CONDA_PREFIX` if you use Conda/Mamba/Anaconda
    
    or

    `meson setup build --prefix=$VIRTUAL_ENV` if you use a Python virtual environment

    or

    `meson setup build --prefix=$HOME` if you want things installed in your `$HOME/bin`, `$HOME/lib` directories.

    or, if you have sudo permissions and want to install to `/usr/local` (or equivalent):

    `meson setup build`.

11. **Check your environment variables against the configuration**

    `python check_meson_build.py`

    If everything looks good, it will tell you. Otherwise, try fixing the issues and starting over from step #8.

12. **Build and install all the C/Fortran codes and the PRESTO shared library (e.g. `libpresto.so`)**

    `meson compile -C build`

    `meson install -C build`

    There should be logs in case anything goes wrong in `$PRESTO/build/meson-logs`. Note that all PRESTO compiled binaries will be installed in `{prefix}/{bindir}`, and the PRESTO shared library (likely either `libpresto.so` or `libpreso.dylib`) will be installed in `{prefix}/{libdir}` as defined by `meson`. You can see the values of `{prefix}`, `{bindir}`, and `{libdir}` using the `check_meson_build.py` script from the previous step.

13. **Try running a PRESTO command like `prepfold`**

    You should get the regular usage screen. If you get a shared library error, see the troubleshooting steps above or below.

14. **Compile and install the PRESTO python codes and libraries**

    `cd $PRESTO/python`

    `pip install .`

    If you get a shared library error, see the troubleshooting steps above or below.

15. **Run some basic tests**

    `cd $PRESTO`

    `python tests/test_presto_python.py`

    `python examplescripts/ffdot_example.py`

    `python python/fftfit_src/test_fftfit.py`

    Another good test is to see if you can run and fit the default profile in `pygaussfit.py`

16. **Run `makewisdom` to have (slightly) fast FFTs**

    Just run `makewisdom` (it is installed to `{prefix}/bin`). By default it uses `FFTW_MEASURE` and skips the largest FFT sizes, so it only takes a few seconds. For higher-quality wisdom over the full set of sizes, run `makewisdom -patient`, which uses `FFTW_PATIENT` and takes about 10-20 min, so be patient. Either way it writes `fftw_wisdom.txt` directly into `{prefix}/share/presto` (or `$PRESTO/lib` if you have `PRESTO` set), which is where PRESTO looks for it at runtime, so no manual move is needed. To write it elsewhere, pass an output path as an argument.

17. **Go find pulsars!**
    
    Everything should be ready to go now, and installed (likely) in the same place as the rest of your Python virtual environment and/or Conda/Mamba/Anaconda environment.

    If you want to uninstall everything you can do:

    `cd $PRESTO/build`

    `ninja uninstall`

    `pip uninstall presto` 

Scott Ransom
Updated Oct 2025, for v5.2.0

-----------------------------------------------------------------

## TROUBLESHOOTING
--------------------

Couple quick trouble-shooting tips if you are having problems compiling and running:

1. Environment variables!  Use `python check_meson_build.py`!
   - Is `PRESTO` set to the top-level PRESTO source directory?
   - Is `TEMPO2` set to the tempo2 runtime directory (only needed for polycos; automatic with conda-forge's tempo2)?
   - Is `PGPLOT_DIR` set to the location of the PGPLOT utility files?
     (Note: On Ubuntu, that should be `/usr/lib/pgplot5`)
   - Is `$PRESTO/bin` in your `PATH`? (It should *not* be!)
   - Is `$PRESTO/lib` in your `LD_LIBRARY_PATH`? (It should *not* be!)
   - Is there any `presto` stuff in your `PYTHONPATH`? (There should *not* be!)

2. Have you have installed the relevant `-dev` packages for `glib2`, `FFTW`, `GSL`, and `CFITSIO` if you are using a Debian-based Linux distribution? Here are the required packages on a clean Ubuntu:
   - `git`
   - `build-essential`
   - `libfftw3-bin`
   - `libfftw3-dev`
   - `libgsl28`
   - `libgsl-dev`
   - `pgplot5`
   - `libglib2.0-dev`
   - `libcfitsio-bin`
   - `libcfitsio-dev`
   - `libpng-dev`
   - `ghostscript`
   - `gfortran`
   - `tcsh`
   - `autoconf`
   - `libx11-dev`
   - `python3-dev`
   - `python3-numpy`
   - `python3-pip`
   
   And the following command should get all of them: `apt install git build-essential libfftw3-bin libfftw3-dev libgsl28 libgsl-dev liberfa-dev pgplot5 libglib2.0-dev libcfitsio-bin libcfitsio-dev libpng-dev ghostscript gfortran tcsh autoconf libx11-dev python3-dev python3-numpy python3-pip`

3. After the Python modules are built and installed, and you run `python tests/test_presto_python.py`, if you get a memory error, please contact Scott! I think that these issues are fixed, but if they are not, we will need to change the build process a tiny bit with a special variable define.
   
4. If you are having trouble with PRESTO creating polycos, you can use `prepfold` with the `-debug` option when folding using `-timing`. That will show you the `tempo2` call and keep all of the (usually) temporary output files.

5. If the python build or install is failing, it can be useful to save and inspect the build directory via `pip install --config-settings=builddir=build .`.

## For MacOS Users
-------------------

If you are using **MacOS**, Paul Ray has been running PRESTO a lot and knows several tricks to get it working:

- PRESTO should build almost "out of the box" on a Mac, once you have the external packages installed and after setting a few environment variables.

- For MacPorts, which has worked well, this should install all the important packages. **This assumes you use Python 3.11 and gcc13.  You may use different versions.**:
    ~~~
    % sudo port install mp-gcc13 python311 py311-ipython pip311 virtualenv311 virtualenvwrapper311
    % sudo port select --set gcc mp-gcc13
    % sudo port select --set python python311
    % sudo port select --set virtualenv virtualenv311
    % sudo port select --set  virtualenvwrapper virtualenvwrapper311
    % sudo port install pgplot cfitsio gsl glib2 fftw-3 fftw-3-single
    ~~~
- tempo2 (only needed for polycos) is available from conda-forge for osx-64; on osx-arm64 you currently need to build it from source.
- Before you build, you will likely need to set the following environment variables. You probably do *not* need to have `DYLD_LIBRARY_PATH` set at runtime.
    ~~~
    # These are needed only at *BUILD* time
    # This points to the MacPorts libraries and those installed in your virtualenv, 
    # as well as a manually installed PGPLOT (not needed if PGPLOT is installed by MacPorts)
    export LIBRARY_PATH=/opt/local/lib:/opt/local/lib/libgcc:$VIRTUAL_ENV/lib:$PGPLOT_DIR

    # This prevents using the macOS native "cc" command, in favor of the MacPorts gcc
    export CC=gcc
    # This makes sure the MacPorts includes can be found
    export CFLAGS="-I/opt/local/include -I$PGPLOT_DIR"
    ~~~

- If you are using Homebrew instead, Alessandro Ridolfi has successfully installed PRESTO using that, as well. Please let Alessandro or I know if you are having issues and we can send you some help. 

## For Docker / Singularity Users
----------------------------------

Alessandro Ridolfi wrote recipes to create Docker images which you can use and modify in the `examplescripts` directory.

He made two versions: one with PNG image support (i.e. providing `gs` via the `ghostscript` package, which `prepfold` uses to make the nice anti-aliased .png versions of the .pfd.ps files) and the other one without the PNG support. Note that `prepfold` calls `gs` only if it is present; if it is missing, the .ps plot is still written and no error is reported. (Earlier versions of PRESTO used `pstoimg` from the large `latex2html` package for this; `pstoimg` was itself only a wrapper around Ghostscript, which is not available on conda-forge, so PRESTO now calls `gs` directly.) The files are:
- examplescripts/Dockerfile_presto5_png_ubuntu24.04.txt
- examplescripts/Dockerfile_presto5_ubuntu24.04.txt

He has also made available pre-built Docker and Singularity images which can be grabbed from his DockerHub with the following commands:

- Docker images:
  ~~~
  docker pull alex88ridolfi/presto5:latest
  docker pull alex88ridolfi/presto5:png
  ~~~

- Apptainer / Singularity images:
  ~~~
  singularity pull docker://alex88ridolfi/presto5:latest
  singularity pull docker://alex88ridolfi/presto5:png
  ~~~

These can  be used with commands like:
  ~~~
  singularity shell -B /path/to/mount:/path/to/mount presto5_latest.sif
  singularity shell -B /path/to/mount:/path/to/mount presto5_png.sif
  ~~~

These images could be of help for people who don’t want / cannot install PRESTO by themselves, or for use on HPC systems.

Alessandro is also planning to make images that include PRESTO5+PULSAR_MINER, and PRESTO5+PRESTO_GPU+PULSAR_MINER in the near future.
