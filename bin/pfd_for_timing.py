#!/usr/bin/env python
from __future__ import print_function
import sys
from presto import prepfold


if len(sys.argv) == 1:
    sys.stderr.write("""usage:  pfd_for_timing.py PFDFILES\n
    This script returns 'true' or 'false' if a .pfd file can be
    used for timing via get_TOAs.py or not.\n""")
    sys.exit(0)

for pfdfile in sys.argv[1:]:
    try:
        pfd = prepfold.pfd(pfdfile)
        if pfd.use_for_timing():
            print("%s: true"%pfdfile)
        else:
            print("%s: false"%pfdfile)
    except:
        sys.stderr.write("Error:  Can't check '%s'\n"%pfdfile)
