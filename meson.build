project('PRESTO', 'c',
  version: '5.0.0',
  license: 'GPL-2.0',
  default_options: [
    'buildtype=release',
    'prefix=/home/sransom'
  ],
)
# Other buildtype possibilities:
#    'buildtype=debugoptimized',
#    'buildtype=debug',
# Can also set install prefix as option via e.g. 'prefix=/home/sransom'

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

# These next couple things are from the Scipy top-level meson.build
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

if host_machine.system() == 'darwin' and cc.has_link_argument('-Wl,-ld_classic')
  # New linker introduced in macOS 14 not working yet, see gh-19357 and gh-19387
  add_project_link_arguments('-Wl,-ld_classic', language : ['c', 'cpp', 'fortran'])
endif

subdir('src')
subdir('bin')

import('pkgconfig').generate(
  libpresto,
  filebase: 'presto',
  name: 'PRESTO',
  description: 'PulsaR Exploration and Search TOolkit by Scott Ransom',
)
