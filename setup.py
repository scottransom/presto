from __future__ import print_function
import os
import sys
import numpy
from numpy.distutils.core import Extension, setup


version = "2.2"

define_macros = []
undef_macros = []
extra_compile_args = ["-DUSEFFTW"]
include_dirs = [numpy.get_include()]

# Note:  you might need to add "gfortran" to the following list if 
# you see errors relating to missing "g" functions....
ppgplot_libraries = ["cpgplot", "pgplot", "X11", "png", "m"]
ppgplot_library_dirs = ["/usr/X11R6/lib"]

presto_libraries = ["presto", "fftw3f", "m"]
presto_library_dirs = []

ppgplot_include_dirs = include_dirs
presto_include_dirs = include_dirs
undef_macros.append('USE_NUMARRAY')

if os.name != "posix":
    raise Exception("os not supported")

if "PGPLOT_DIR" in os.environ:
    ppgplot_library_dirs.append(os.environ["PGPLOT_DIR"])
    ppgplot_include_dirs.append(os.environ["PGPLOT_DIR"])
else:
    print("PGPLOT_DIR env var not defined!", file=sys.stderr)
if "PRESTO" in os.environ:
    presto_library_dirs.append(os.path.join(os.environ["PRESTO"], "lib"))
    presto_include_dirs.append(os.path.join(os.environ["PRESTO"], "include"))
else:
    print("PRESTO env var not defined!", file=sys.stderr)
    presto_include_dirs.append(os.path.join(os.path.dirname(__file__),
                                            'include'))

ext_ppgplot = Extension('_ppgplot',
                        ['python/ppgplot_src/_ppgplot.c'],
                        include_dirs=ppgplot_include_dirs,
                        libraries=ppgplot_libraries,
                        library_dirs=ppgplot_library_dirs,
                        define_macros=define_macros,
                        extra_compile_args=extra_compile_args)

ext_presto = Extension('_presto',
                       ['python/presto_src/presto_wrap.c'],
                       include_dirs=presto_include_dirs,
                       libraries=presto_libraries,
                       library_dirs=presto_library_dirs,
                       define_macros=define_macros,
                       extra_compile_args=extra_compile_args)

# not working yet
ext_fftfit = Extension('fftfit', sources=['python/fftfit_src/brent.f',
                                          'python/fftfit_src/bcprof.f',
                                          'python/fftfit_src/bfccf.f',
                                          'python/fftfit_src/bffft.f',
                                          'python/fftfit_src/bfftfit.f',
                                          'python/fftfit_src/bfftfit.pyf'])


scripts = ['bin/' + i for i in os.listdir('bin')]


setup(name="presto",
      version=version,
      install_requires=['numpy', 'future', 'six', 'scipy', 'matplotlib', 'astropy', 'pyslalib'],
      scripts=scripts,
      description="Python interfaces to PGPLOT and PRESTO",
      author="Scott Ransom (ppgplot from Nick Patavlis)",
      author_email="sransom@nrao.edu",
      url="https://github.com/scottransom/presto",
      packages=['presto', 'presto.ppgplot', 'presto.presto', 'presto.singlepulse'],
      package_dir={'presto.ppgplot': 'python/ppgplot_src',
                   'presto.presto': 'python/presto_src',
                   'presto': 'python/presto',
                   'presto.singlepulse': 'python/presto/singlepulse',
                   },
      ext_modules=[ext_ppgplot, ext_presto]) #, ext_fftfit])
