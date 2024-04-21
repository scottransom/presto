# Note: will need to fix Numpy API deprecation soon -SMR
# c_args: [numpy_nodepr_api, '-Wno-unused-variable'],
py3.extension_module('_ppgplot', '_ppgplot.c',
  include_directories: inc_np,
  dependencies : [py3_dep, libm, pgplot, cpgplot, x11, png],
  install : true,
)

py3.install_sources(['__init__.py'],
  subdir: 'presto/ppgplot'
)
