from __future__ import print_function
from Numeric import *
from cvects import *

def carr(n):
    return complex_arr(n)

def multcarr(a, val):
    mult_arr(a, val, len(a))

a = open("cvects.h", "r")
print(get_filelen(a))
a.close()

d = arange(10, typecode='d')
print(d)
dgenrotate_1d(d, len(d), 3.5)
print(d)
print(dstats(d, len(d)))

f = complex_arr(10)
print(f)
mult_arr(f, 10.0, len(f))
print(f)
