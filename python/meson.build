project('PRESTO_python', 'c', 'fortran',
  version: '5.0.0',
  license: 'GPL-2.0',
  default_options: [
    'buildtype=release',
  ],
)
#
# Other possibilities
#    'buildtype=debugoptimized',
#    'buildtype=debug',

cc = meson.get_compiler('c')

# Compile and prepare all PRESTO code and python packages
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()

# Check the python headers
cc.check_header('Python.h', dependencies: [py3_dep], required: true)
glib = dependency('glib-2.0')
fftw = dependency('fftw3f') # Note the 'f' for the single-precision version!
fits = dependency('cfitsio')
x11 = dependency('x11')
png = dependency('libpng')

libm = cc.find_library('m', required: false)
libpresto = cc.find_library('presto', required: true)
pgplot = cc.find_library('pgplot', required: true)
cpgplot = cc.find_library('cpgplot', required: true)

inc = include_directories('../include')

incdir_numpy = run_command(py3,
  [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy)

# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

subdir('fftfit_src')
subdir('ppgplot_src')
subdir('presto_src')
subdir('presto')
