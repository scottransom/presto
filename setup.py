from __future__ import print_function
import os
import sys
import numpy

# setuptools has to be imported before numpy.distutils.core
import setuptools
from numpy.distutils.core import Extension, setup

version = "3.0.1"

define_macros = []
undef_macros = []
extra_compile_args = ["-DUSEFFTW"]
include_dirs = [numpy.get_include()]
# For MacOS with MacPorts use the following
# include_dirs.append("/opt/local/include")

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
    extra_link_args = ["-Wl,-rpath,"+os.path.join(os.environ["PRESTO"], "lib")]
else:
    print("PRESTO env var not defined!", file=sys.stderr)
    presto_include_dirs.append(os.path.join(os.path.dirname(__file__),
                                            'include'))
    extra_link_args = ["-Wl,-rpath,"+os.path.join(os.path.dirname(__file__),
                                                  "lib")]

ext_ppgplot = Extension('_ppgplot',
                        ['python/ppgplot_src/_ppgplot.c'],
                        include_dirs=ppgplot_include_dirs,
                        libraries=ppgplot_libraries,
                        library_dirs=ppgplot_library_dirs,
                        define_macros=define_macros,
                        extra_compile_args=extra_compile_args,
                        extra_link_args=extra_link_args)

ext_presto = Extension('_presto',
                       ['python/presto_src/presto_wrap.c'],
                       include_dirs=presto_include_dirs,
                       libraries=presto_libraries,
                       library_dirs=presto_library_dirs,
                       define_macros=define_macros,
                       extra_compile_args=extra_compile_args,
                       extra_link_args=extra_link_args)

ext_fftfit = Extension('_fftfit', sources=['python/fftfit_src/brent.f',
                                          'python/fftfit_src/cprof.f',
                                          'python/fftfit_src/fccf.f',
                                          'python/fftfit_src/ffft.f',
                                          'python/fftfit_src/fftfit.f',
                                          'python/fftfit_src/_fftfit.pyf'])


scripts = ['bin/' + i for i in os.listdir('bin') if i.endswith('.py') or i.endswith('.sh')]


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
      package_data={'presto': ['*.json']},
      ext_modules=[ext_ppgplot, ext_presto, ext_fftfit])
