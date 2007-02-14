from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc, get_python_lib
import os
import sys
import numpy

###################################################################
# build the extension
#

define_macros = []
undef_macros = []
extra_compile_args = []
include_dirs = []

ppgplot_libraries = ["cpgplot", "pgplot", "X11", "png", "m", "g2c"]
ppgplot_library_dirs = ["/usr/X11R6/lib"]

# Use NumPy instead of Numeric or numarray
make_extension = Extension
#include_dirs.append(numpy.get_numpy_include())
include_dirs.append(numpy.get_include())
ppgplot_include_dirs = include_dirs
presto_include_dirs = include_dirs
undef_macros.append('USE_NUMARRAY')

if os.name == "posix":
    if os.environ.has_key("PGPLOT_DIR"):
        ppgplot_library_dirs.append(os.environ["PGPLOT_DIR"])
        ppgplot_include_dirs.append(os.environ["PGPLOT_DIR"])
    else:
        print >>sys.stderr, "PGPLOT_DIR env var not defined!"
else:
    raise Exception, "os not supported"

ext_ppgplot = make_extension('_ppgplot',
                             [os.path.join('ppgplot_src','_ppgplot.c')],
                             include_dirs=ppgplot_include_dirs,
                             libraries=ppgplot_libraries,
                             library_dirs=ppgplot_library_dirs,
                             define_macros=define_macros,
                             extra_compile_args=extra_compile_args)


###################################################################
# the package
#

setup(name="ppgplot",
      version="0.99",
      description="Python interface to PGPLOT",
      author="Scott Ransom (ppgplot from Nick Patavlis)",
      author_email="sransom@nrao.edu",
      packages=['ppgplot'],
      package_dir={'ppgplot':'ppgplot_src'},
      ext_modules=[ext_ppgplot])
