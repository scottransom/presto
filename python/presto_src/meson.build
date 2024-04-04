# Note: will need to fix Numpy API deprecation soon -SMR
# c_args: [numpy_nodepr_api, '-Wno-unused-variable'],
py3.extension_module('_presto', 'presto_wrap.c',
  include_directories: [inc, inc_np],
  dependencies : [py3_dep, libm, fftw, libpresto],
  install : true,
)

py3.install_sources(
  ['__init__.py', 'prestoswig.py'],
  subdir: 'presto/presto'
)
