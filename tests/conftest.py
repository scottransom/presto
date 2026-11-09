"""pytest configuration shared by the tests in this directory."""

import os

# In the sanitizer CI job, ASan's runtime is preloaded (LD_PRELOAD) because the Python
# interpreter is not instrumented but the SWIG module's libpresto is.  This process
# already has it loaded, so drop it from the environment that child processes inherit.
# The PRESTO executables are instrumented and need no preload, and forcing ASan into
# uninstrumented programs that they start (e.g. tempo2, via prepfold -timing) can make
# them hang.
_preload = os.environ.get("LD_PRELOAD", "")
if "libasan" in _preload:
    _rest = " ".join(p for p in _preload.replace(":", " ").split() if "libasan" not in p)
    if _rest:
        os.environ["LD_PRELOAD"] = _rest
    else:
        del os.environ["LD_PRELOAD"]
