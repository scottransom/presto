#!/usr/bin/env python
import sys
from presto import sifting


if len(sys.argv) < 2:
    sys.stderr.write("\nusage:  quick_prune_cands.py ACCEL_file_name [sigma]\n\n")
    sys.exit()

if len(sys.argv)==3:
    sifting.sigma_threshold = float(sys.argv[2])

cands = sifting.read_candidates([sys.argv[1]], track=True)
cands.print_cand_summary()
cands.to_file()
