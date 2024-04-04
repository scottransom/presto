fftfit_module = custom_target('fftfit_module',
  output: 'fftfitmodule.c',
  input: 'fftfit.pyf',
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)

py3.extension_module('fftfit',
  [fftfit_module, fortranobject_c, 'brent.f', 'cprof.f', 'fccf.f', 'ffft.f', 'fftfit.f'],
  include_directories: [inc_np, inc_f2py],
  c_args: [numpy_nodepr_api, '-Wno-unused-variable'],
  dependencies : py3_dep,
  link_language: 'fortran',
  subdir: 'presto',
  install : true,
)
