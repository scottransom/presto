import numpy as np
import newpresto as presto

print "Testing FFT stuff...",
N = 20
x = np.random.standard_normal(N)
nx = presto.rfft(presto.rfft(x, -1), 1)
assert(np.allclose(x, nx, atol=1e-6))
print "success"

print "Testing FFTW call...",
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.fftwcall(cx, -1)
presto.fftwcall(cx, 1)
assert(np.allclose(cx/N, ncx, atol=1e-6))
print "success"

print "Testing tablesixstepfft call...",
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.tablesixstepfft(cx, -1)
presto.tablesixstepfft(cx, 1)
assert(np.allclose(cx/N, ncx, atol=1e-6))
print "success"

print "Testing reading infiles...",
x = presto.read_inffile("1937_DM71.02_zerodm.inf", verbose=False)
assert(x.telescope=="GBT")
assert(x.mjd_i==55267)
assert(x.dt==8.192e-05)
assert(x.numonoff==1)
assert(x.analyzer=="sransom")
print "success"

x.analyzer="test"
x.name="xxx"
print "Testing writing infiles...",
presto.write_inffile(x, verbose=False)
y = presto.read_inffile("xxx", verbose=False)
assert(y.analyzer=="test")
assert(y.bary==0)
assert(y.numonoff==1)
print "success"

