project('PRESTO', 'c', 'fortran',
  version: '5.0.0',
  license: 'GPL-2.0',
  default_options: [
    'buildtype=release',
  ],
)
# Other possibilities
#    'buildtype=debugoptimized',
#    'buildtype=debug',

cc = meson.get_compiler('c')

glib = dependency('glib-2.0')
fftw = dependency('fftw3f') # Note the 'f' for single-precision version!
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

_global_c_args = cc.get_supported_arguments(
  '-DUSE_FFTW_MALLOC',
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
  '-Wno-unused-result',
  '-Wno-unused-but-set-parameter'
)
add_project_arguments(_global_c_args, language : 'c')

subdir('src')
subdir('bin')

import('pkgconfig').generate(
  libpresto,
  filebase: 'presto',
  name: 'PRESTO',
  description: 'PulsaR Exploration and Search TOolkit by Scott Ransom',
)
