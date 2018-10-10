#!/usr/bin/env python
from __future__ import (print_function,division)
import presto.psr_utils as pu
import sys
from presto.infodata import infodata

if len(sys.argv) != 2:
    print("chooseN <file.inf|numpoints>")
    print("    Prints a good value for fast FFTs to be used for -numout in prepdata/prepsubband")
    sys.exit(1)

if sys.argv[1].endswith('.inf'):
    inf = infodata(sys.argv[1])
    n = inf.N
else:
    try:
        n = int(sys.argv[1])
    except:
        print("chooseN <file.inf|numpoints>")
        print("    Prints a good value for fast FFTs to be used for -numout in prepdata/prepsubband")
        sys.exit(2)

print(pu.choose_N(n))
