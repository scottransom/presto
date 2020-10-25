from __future__ import print_function
import numpy as np
from presto import presto
import os
from os import path

here = path.dirname(__file__)

print("Testing FFT stuff...", end=' ')
N = 20
x = np.random.standard_normal(N)
nx = presto.rfft(presto.rfft(x, -1), 1)
assert(np.allclose(x, nx, atol=1e-6))
print("success")

print("Testing FFTW call...", end=' ')
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.fftwcall(cx, -1)
presto.fftwcall(cx, 1)
assert(np.allclose(cx/N, ncx, atol=1e-6))
print("success")

print("Testing tablesixstepfft call...", end=' ')
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.tablesixstepfft(cx, -1)
presto.tablesixstepfft(cx, 1)
assert(np.allclose(cx/N, ncx, atol=1e-6))
print("success")

print("Testing reading infiles...", end=' ')
x = presto.read_inffile(path.join(here, "1937_DM71.02_zerodm.inf"), verbose=False)
assert(x.telescope=="GBT")
assert(x.mjd_i==55267)
assert(x.dt==8.192e-05)
assert(x.numonoff==1)
assert(x.analyzer=="sransom")
print("success")

print("Testing writing infiles...", end=' ')
x.analyzer="test"
x.name="xxx"
x.dt=0.125
presto.write_inffile(x, verbose=False)
y = presto.read_inffile("xxx", verbose=False)
assert(y.analyzer=="test")
assert(y.bary==0)
assert(y.numonoff==1)
assert(y.dt==0.125)
os.remove("xxx.inf")
print("success")

print("Testing allocation and freeing of memory...", end=' ')
for ii in range(1024):
    a = presto.gen_fvect(1024 * 32768)
    del a
for ii in range(1024):
    a = presto.gen_cvect(1024 * 16384)
    del a
print("success")

print("Testing psrparams and orbitparams stuff...", end=' ')
psr = presto.psrepoch("J0737-3039A", 56000.0, verbose=False)
assert(round(psr.dm-48.92, 7)==0)
# This needs to change when we start using the actual psrcat.db file
assert(round(psr.orb.p-8834.534998272, 7)==0)
print("success")

print("Testing spectralpower and spectralphase...", end=' ')
a = np.arange(5.0) + complex(0.0, 1.0)
assert(np.allclose(presto.spectralpower(a),
                   np.arange(5.0)**2.0 + 1))
assert(np.allclose(presto.spectralphase(a),
                   np.array([90., 45., 26.56505203, 18.43494797, 14.03624344])))
print("success")

print("Testing vector shifting / rotation...", end=' ')
a = np.arange(4, dtype=np.float32)
presto.frotate(a, 1)
assert(np.allclose(a, np.array([1, 2, 3, 0])))
a = np.arange(4, dtype=np.float64)
presto.drotate(a, 1)
assert(np.allclose(a, np.array([1, 2, 3, 0])))
print("success")

print("Testing orbit integration stuff...", end=' ')
orb = presto.orbitparams()
orb.p = 10000.0
orb.e = 0.1
orb.x = 1.0
orb.t = 1234.0
orb.w = 75.0
orb.pd = orb.wd = 0.0
E0 = presto.keplers_eqn(orb.t+0.0, orb.p, orb.e, 1e-15)
E1 = presto.keplers_eqn(orb.t+100.0, orb.p, orb.e, 1e-15)
E2 = presto.keplers_eqn(orb.t+200.0, orb.p, orb.e, 1e-15)
E3 = presto.keplers_eqn(orb.t+300.0, orb.p, orb.e, 1e-15)
Es = np.asarray([E0, E1, E2, E3])
Es_check = np.asarray([ 0.85050653, 0.9175909,
                        0.9842971, 1.05061346])
assert(np.allclose(Es, Es_check))
Es_new = presto.dorbint(E0, 4, 100.0, orb)
assert(np.allclose(Es_new, Es_check))
presto.E_to_v(Es, orb)
Vs_check = np.asarray([-112.15558594, -122.45299212,
                       -131.9991447, -140.76659065])
assert(np.allclose(Es, Vs_check))
minv, maxv = presto.binary_velocity(300.0, orb)
minv *= presto.SOL/1000.0
maxv *= presto.SOL/1000.0
assert(round(minv-Vs_check.min(), 7)==0)
assert(round(maxv-Vs_check.max(), 7)==0)
print("success")

print("Testing Fourier response generation...", end=' ')
numbetween = 16
z = 5.0
w = 40.0
bins_per_side = max([presto.r_resp_halfwidth(presto.LOWACC), \
                     presto.z_resp_halfwidth(z, presto.LOWACC), \
                     presto.w_resp_halfwidth(z, w, presto.LOWACC)])
nn = numbetween * bins_per_side * 2;
rresp = presto.gen_r_response(0.0, numbetween, nn)
zresp = presto.gen_z_response(0.0, numbetween, nn, z)
wresp = presto.gen_w_response(0.0, numbetween, nn, z, w)
pr = presto.spectralpower(rresp)
pz = presto.spectralpower(zresp)
pw = presto.spectralpower(wresp)
rs = np.arange(float(nn))/numbetween - bins_per_side
if False:
    import matplotlib.pyplot as plt
    plt.plot(rs, pr, 'b-')
    plt.plot(rs, pz, 'g-')
    plt.plot(rs, pw, 'r-')
    plt.show()
assert(rs[nn//2]==0.0)
assert(pr[nn//2]==1.0)
assert(round(pz[nn//2]-0.227675, 6)==0)
assert(round(pw[nn//2]-0.019462, 6)==0)
print("success")

print("Testing angle functions...", end=' ')
dd1 = 15.25
dd2 = presto.dms2rad(*presto.deg2dms(dd1))*presto.RADTODEG
assert(round(dd1-dd2, 12)==0)
dd1 = -0.5
dd2 = presto.dms2rad(*presto.deg2dms(dd1))*presto.RADTODEG
assert(round(dd1-dd2, 12)==0)
hh1 = 12.125
hh2 = presto.hms2rad(*presto.hours2hms(hh1))*presto.RADTODEG/15.0
assert(round(hh1-hh2, 12)==0)
hh1 = -0.5
hh2 = presto.hms2rad(*presto.hours2hms(hh1))*presto.RADTODEG/15.0
assert(round(hh1-hh2, 12)==0)
ang = presto.sphere_ang_diff(10.0*presto.DEGTORAD, 10.0*presto.DEGTORAD,
                             50.0*presto.DEGTORAD, -10.0*presto.DEGTORAD)
assert(round(160334.960*presto.ARCSEC2RAD-ang, 7)==0)
print("success")

# Only run this test if TEMPO is available
envval = os.getenv("TEMPO")
if envval is not None:
    print("Testing get_baryv (barycenter)...", end=' ')
    vavg1 = presto.get_baryv("18:24:32.9520", "-24:52:12.0000",
                             56421.44222222222222, 214.5386496, obs="GB")
    vavg2 = -7.2069293455783169e-05
    assert(round(vavg1-vavg2, 10)==0)
    print("success")
else:
    print("Skipping test of presto.get_baryv() since TEMPO not set.")

print("Testing simple folding code...", end=' ')
prof, phs = presto.fold(np.ones(10000), 0.001, 10, 1)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 200, 1)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 500, 1)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, tlo=-1.12313)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, startphs=0.65765)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, startphs=0.99)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, startphs=0.99, standard=False)
assert(np.allclose(prof, prof.mean()))
assert(np.all(prof>0))
print("success")
