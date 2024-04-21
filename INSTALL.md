# New and improved steps to install PRESTO!

## td;dr for experienced PRESTO builders:
With v5, we have switched to building and installing with [meson](https://mesonbuild.com/).

**MacOS users should see the comments at the bottom of this file!**

As always, there are a set of essential packages required to build PRESTO. This command should do it on a Debian/Ubuntu-like system:
`apt install git build-essential libfftw3-bin libfftw3-dev pgplot5 libglib2.0-dev libcfitsio-bin libcfitsio-dev libpng-dev gfortran tcsh autoconf libx11-dev python3-dev python3-numpy python3-pip`

Make sure that your `PRESTO` environment variable points to the top-level PRESTO git checkout. And make sure that `$PRESTO/lib` and `$PRESTO/bin` are **not** in your `PATH` or `LD_LIBRARY_PATH` or `PYTHONPATH` environment variables as we have required in the past. It is probably a good idea to clean your earlier compiles, as well. Just cd into the `src` directory and do a `make cleaner`, and then come back here.

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
    pip install --config-settings=builddir=build .

And that should do it! You can quickly test to see if most things are working by doing:

    cd $PRESTO
    python tests/test_presto_python.py
    python examplescripts/ffdot_example.py
    python python/fftfit_src/test_fftfit.py

Another good test is to see if you can run and fit the default profile in `pygaussfit.py`

If you want to run `makewisdom` for slightly faster FFT calls, it is located in `$PRESTO/build/src`. Just run it from there, and then copy or move the resulting `fftw_wisdom.txt` file to `$PRESTO/lib`.

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

    I highly recommend that you use pre-compiled packages for your OS/distribution! FOr example, Ubuntu has good FFTW packages: `libfftw3-bin` and `libfftw3-dev`.

    If you compile your own, you need to compile FFTW for **single** precision. For all architectures I recommend the following configuration: `./configure --enable-shared --enable-single`

    If you are on a modern Intel processor and have a recent version of GCC, you can get much improved performance by adding: `--enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-fma`

2.  **Install [PGPLOT](http://www.astro.caltech.edu/~tjp/pgplot/)**

    I highly recommend that you use pre-compiled packages for your OS/distribution! FOr example, Ubuntu has a good PGPLOT package: `pgplot5`. You will likely need to set the `PGPLOT_DIR` environment variable. On Ubuntu, I have `PGPLOT_DIR=/usr/lib/pgplot5`

    If you want to try to compile your own (good luck!), you need the X-windows and postscript drivers at a minimum.

3.  **Install [TEMPO](http://tempo.sourceforge.net/)**

    Make sure to set the `TEMPO` environment variable so that it points to the top level of the `TEMPO` code base.

4.  **Install [GLIBv2](http://library.gnome.org/devel/glib/)**

    On Linux machines this is almost certainly already on your system (check in `/usr/lib` and `/usr/include/glib*`). Although you may need to install a glib development package in order to have the required include files. On Ubuntu, the package you need is: `libglib2.0-dev`

5.  **Install [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/)**

    I highly recommend using pre-compiled packages, once again (on Ubuntu they are `libcfitsio-bin` and `libcfitsio-dev`), however, this is a very easy install via source.

6.  **Set the `PRESTO` environment variable**

    It should be set to the top level directory of the PRESTO distribution (i.e. this directory). And make sure that `$PRESTO/lib` and `$PRESTO/bin` are **not** in your `PATH` or `LD_LIBRARY_PATH` or `PYTHONPATH` environment variables as we have required in the past.

7.  **Activate your Python virtual environment *or* Conda/Mamba/Anaconda environment**

    * That environment should have `numpy` installed at a minimum
    * Make sure that `pip` is recent (`pip install --upgrade pip`)
    * Install the build tools: `pip install meson meson-python ninja` or `conda install meson meson-python ninja`

8.  **Configure the meson build**

    In the top level PRESTO directory (i.e. `cd $PRESTO`), configure `meson` via:

    `meson setup build --prefix=$CONDA_PREFIX` if you use Conda/Mamba/Anaconda
    
    or

    `meson setup build --prefix=$VIRTUAL_ENV` if you use a Python virtual environment

    or

    `meson setup build --prefix=$HOME` if you want things installed in your `$HOME/bin`, `$HOME/lib` directories.

    or, if you have sudo permissions and want to install to `/usr/local` (or equivalent):

    `meson setup build`.

9. **Check your environment variables against the configuration**

    `python check_meson_build.py`

    If everything looks good, it will tell you. Otherwise, try fixing the issues and starting over from step #8.

10. **Build and install all the C/Fortran codes and the PRESTO shared library (e.g. `libpresto.so`)**

    `meson compile -C build`

    `meson install -C build`

    There should be logs in case anything goes wrong in `$PRESTO/build/meson-logs`. Note that all PRESTO compiled binaries will be installed in `{prefix}/{bindir}`, and the PRESTO shared library (likely either `libpresto.so` or `libpreso.dylib`) will be installed in `{prefix}/{libdir}` as defined by `meson`. You can see the values of `{prefix}`, `{bindir}`, and `{libdir}` using the `check_meson_build.py` script from the previous step.

11. **Try running a PRESTO command like `prepfold`**

    You should get the regular usage screen. If you get a shared library error, see the troubleshooting steps above or below.

12. **Compile and install the PRESTO python codes and libraries**

    `cd $PRESTO/python`

    `pip install --config-settings=builddir=build .`

    If you get a shared library error, see the troubleshooting steps above or below.

13. **Run some basic tests**

    `cd $PRESTO`

    `python tests/test_presto_python.py`

    `python examplescripts/ffdot_example.py`

    `python python/fftfit_src/test_fftfit.py`

    Another good test is to see if you can run and fit the default profile in `pygaussfit.py`

14. **Run `makewisdom` to have (slightly) fast FFTs**

    Just run `$PRESTO/build/src/makewisdom`. It takes about 10-20 min to run, so be patient. Note that the `fftw_wosdom.txt` file will be located in `$PRESTO/build/src`, so you will need to move it to `$PRESTO/lib` so that PRESTO can find it.

15. **Go find pulsars!**
    
    Everything should be ready to go now, and installed (likely) in the same place as the rest of your Python virtual environment and/or Conda/Mamba/Anaconda environment.

    If you want to uninstall everything you can do:

    `cd $PRESTO/build`

    `ninja uninstall`

    `pip uninstall presto` 

Scott Ransom
Updated April 2024, for v5.0.0

-----------------------------------------------------------------

## TROUBLESHOOTING
--------------------

Couple quick trouble-shooting tips if you are having problems compiling and
running:

1. Environment variables!
   - Is `PRESTO` set to the top-level PRESTO source directory?
   - Is `TEMPO` set to the top-level TEMPO source directory?
   - Is `PGPLOT_DIR` set to the location of the PGPLOT utility files?
     (Note: On Ubuntu, that should be `/usr/lib/pgplot5`)
   - Is `$PRESTO/bin` in your `PATH`? (It should *not* be!)
   - Is `$PRESTO/lib` in your `LD_LIBRARY_PATH`? (It should *not* be!)
   - Is there any `presto` stuff in your `PYTHONPATH`? (There should *not* be!)

2. Have you have installed the relevant `-dev` packages for `glib2`, `FFTW` and `CFITSIO` if you are using a Debian-based Linux distribution? Here are the required packages on a clean Ubuntu:
   - `git`
   - `build-essential`
   - `libfftw3-bin`
   - `libfftw3-dev`
   - `pgplot5`
   - `libglib2.0-dev`
   - `libcfitsio-bin`
   - `libcfitsio-dev`
   - `libpng-dev`
   - `gfortran`
   - `tcsh`
   - `autoconf`
   - `libx11-dev`
   - `python3-dev`
   - `python3-numpy`
   - `python3-pip`
   
   And the following command should get all of them: `apt install git build-essential libfftw3-bin libfftw3-dev pgplot5 libglib2.0-dev libcfitsio-bin libcfitsio-dev libpng-dev gfortran tcsh autoconf libx11-dev python3-dev python3-numpy python3-pip`

3. After the Python modules are built and installed, and you run `python tests/test_presto_python.py`, if you get a memory error, please contact Scott! I think that these issues are fixed, but if they are not, we will need to change the build process a tiny bit with a special variable define.
   
4. If you are having trouble with PRESTO creating polycos, you can use `prepfold` with the `-debug` option when folding using `-timing`. That will show you the `TEMPO` call and keep all of the (usually) temporary output files.

5. If you are using **MacOS**, Paul Ray has been running PRESTO a lot and knows several tricks to get it working:

    - PRESTO should build almost "out of the box" on a Mac, once you have the external packages installed and after setting a few environment variables.
    
    - For MacPorts, which has worked well, this should install all the important packages. **This assumes you use Python 3.11 and gcc13.  You may use different versions.**:

        ~~~
        % sudo port install mp-gcc13 python311 py311-ipython pip311 virtualenv311 virtualenvwrapper311
        % sudo port select --set gcc mp-gcc13
        % sudo port select --set python python311
        % sudo port select --set virtualenv virtualenv311
        % sudo port select --set  virtualenvwrapper virtualenvwrapper311
        % sudo port install pgplot cfitsio glib2 fftw-3 fftw-3-single
        ~~~

    - TEMPO should build easily with gfortran. I did not make any changes to the distro.

    - Before you build, you will likely need to set the following environment variables. You probably do *not* need to have `DYLD_LIBRARY_PATH` set at runtime.

        ~~~
        # These are needed only at *BUILD* time
        # This points to the MacPorts libraries and those installed in your virtualenv
        export LIBRARY_PATH=/opt/local/lib:/opt/local/lib/libgcc
        # This prevents using the macOS native "cc" command, in favor of the MacPorts gcc
        export CC=gcc
        # This makes sure the MacPorts includes can be found
        export CFLAGS="-I/opt/local/include"
        ~~~