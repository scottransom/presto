project('PRESTO', 'c', 'fortran',
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

glib = dependency('glib-2.0')
fftw = dependency('fftw3f') # Note the 'f' for the single-precision version!
fits = dependency('cfitsio')
x11 = dependency('x11')
png = dependency('libpng')
omp = dependency('openmp', required: false)
mpi = dependency('mpi', language: 'c', required: false,
  not_found_message: 'MPI not found. Skipping mpiprepsubband build.')

libm = cc.find_library('m', required: false)
pgplot = cc.find_library('pgplot', required: true)
cpgplot = cc.find_library('cpgplot', required: true)

inc = include_directories('include')

# Compile and prepare all PRESTO code and python packages
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()

subdir('src')
subdir('bin')

import('pkgconfig').generate(
  libpresto,
  filebase: 'presto',
  name: 'PRESTO',
  description: 'PulsaR Exploration and Search TOolkit by Scott Ransom',
)
