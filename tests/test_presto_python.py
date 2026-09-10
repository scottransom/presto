import numpy as np
import presto.events  # tests whether cosine_rand.json was installed correctly
from presto import presto
import os
from os import path

here = path.dirname(__file__)

print("Testing FFT stuff...", end=" ")
N = 20
x = np.random.standard_normal(N)
nx = presto.rfft(presto.rfft(x, -1), 1)
assert np.allclose(x, nx, atol=1e-6)
print("success")

print("Testing FFTW call...", end=" ")
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.fftwcall(cx, -1)
presto.fftwcall(cx, 1)
assert np.allclose(cx / N, ncx, atol=1e-6)
print("success")

print("Testing tablesixstepfft call...", end=" ")
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.tablesixstepfft(cx, -1)
presto.tablesixstepfft(cx, 1)
assert np.allclose(cx / N, ncx, atol=1e-6)
print("success")

print("Testing reading infiles...", end=" ")
x = presto.read_inffile(path.join(here, "1937_DM71.02_zerodm.inf"), verbose=False)
assert x.telescope == "GBT"
assert x.mjd_i == 55267
assert x.dt == 8.192e-05
assert x.numonoff == 1
assert x.analyzer == "sransom"
print("success")

print("Testing writing infiles...", end=" ")
x.analyzer = "test"
x.name = "xxx"
x.dt = 0.125
presto.write_inffile(x, verbose=False)
y = presto.read_inffile("xxx", verbose=False)
assert y.analyzer == "test"
assert y.bary == 0
assert y.numonoff == 1
assert y.dt == 0.125
os.remove("xxx.inf")
print("success")

print("""Note: If the following fails, edit the top-level meson.build and remove
        '-DUSE_FFTW_MALLOC' from the `_global_c_args` definition, and re-build.""")
print("Testing allocation and freeing of memory...", end=" ")
for ii in range(1024):
    a = presto.gen_fvect(1024 * 32768)
    del a
for ii in range(1024):
    a = presto.gen_cvect(1024 * 16384)
    del a
print("success")

print("Testing psrparams and orbitparams stuff...", end=" ")
psr = presto.psrepoch("J0737-3039A", 56000.0, verbose=False)
assert round(psr.dm - 48.916, 7) == 0
# This needs to change when we start using the actual psrcat.db file
assert round(psr.orb.p - 8834.534998272, 7) == 0
print("success")

print("Testing spectralpower and spectralphase...", end=" ")
a = np.arange(5.0) + complex(0.0, 1.0)
assert np.allclose(presto.spectralpower(a), np.arange(5.0) ** 2.0 + 1)
assert np.allclose(
    presto.spectralphase(a),
    np.array([90.0, 45.0, 26.56505203, 18.43494797, 14.03624344]),
)
print("success")

print("Testing vector shifting / rotation...", end=" ")
a = np.arange(4, dtype=np.float32)
presto.frotate(a, 1)
assert np.allclose(a, np.array([1, 2, 3, 0]))
a = np.arange(4, dtype=np.float64)
presto.drotate(a, 1)
assert np.allclose(a, np.array([1, 2, 3, 0]))
print("success")

print("Testing orbit integration stuff...", end=" ")
orb = presto.orbitparams()
orb.p = 10000.0
orb.e = 0.1
orb.x = 1.0
orb.t = 1234.0
orb.w = 75.0
orb.pd = orb.wd = 0.0
E0 = presto.keplers_eqn(orb.t + 0.0, orb.p, orb.e, 1e-15)
E1 = presto.keplers_eqn(orb.t + 100.0, orb.p, orb.e, 1e-15)
E2 = presto.keplers_eqn(orb.t + 200.0, orb.p, orb.e, 1e-15)
E3 = presto.keplers_eqn(orb.t + 300.0, orb.p, orb.e, 1e-15)
Es = np.asarray([E0, E1, E2, E3])
Es_check = np.asarray([0.85050653, 0.9175909, 0.9842971, 1.05061346])
assert np.allclose(Es, Es_check)
Es_new = presto.dorbint(E0, 4, 100.0, orb)
assert np.allclose(Es_new, Es_check)
presto.E_to_v(Es, orb)
Vs_check = np.asarray([-112.15558594, -122.45299212, -131.9991447, -140.76659065])
assert np.allclose(Es, Vs_check)
minv, maxv = presto.binary_velocity(300.0, orb)
minv *= presto.SOL / 1000.0
maxv *= presto.SOL / 1000.0
assert round(minv - Vs_check.min(), 7) == 0
assert round(maxv - Vs_check.max(), 7) == 0
print("success")

print("Testing Fourier response generation...", end=" ")
numbetween = 16
z = 5.0
w = 40.0
bins_per_side = max(
    [
        presto.r_resp_halfwidth(presto.LOWACC),
        presto.z_resp_halfwidth(z, presto.LOWACC),
        presto.w_resp_halfwidth(z, w, presto.LOWACC),
    ]
)
nn = numbetween * bins_per_side * 2
rresp = presto.gen_r_response(0.0, numbetween, nn)
zresp = presto.gen_z_response(0.0, numbetween, nn, z)
wresp = presto.gen_w_response(0.0, numbetween, nn, z, w)
pr = presto.spectralpower(rresp)
pz = presto.spectralpower(zresp)
pw = presto.spectralpower(wresp)
rs = np.arange(float(nn)) / numbetween - bins_per_side
if False:
    import matplotlib.pyplot as plt

    plt.plot(rs, pr, "b-")
    plt.plot(rs, pz, "g-")
    plt.plot(rs, pw, "r-")
    plt.show()
assert rs[nn // 2] == 0.0
assert pr[nn // 2] == 1.0
assert round(pz[nn // 2] - 0.227675, 6) == 0
assert round(pw[nn // 2] - 0.019462, 6) == 0
print("success")

print("Testing Fourier 1-D response maximization...", end=" ")
N = 4096
f = 1023.564
data = 0.1 * np.random.standard_normal(4096) + np.cos(
    2 * np.pi * f * np.linspace(0.0, 1.0, N, endpoint=False)
)
ft = presto.rfft(data, -1)
derivs = presto.rderivs()
rmax, maxpow = presto.max_r_arr(ft, 1023.3, derivs)
predpow = N**2 / 4
assert 0.9 * predpow < maxpow < 1.1 * predpow
assert f - 0.03 < rmax < f + 0.03
print("success")

print("Testing Fourier f-fdot (r, z) maximization...", end=" ")
# These use large |z| on purpose.  The optimizers work internally in the
# scaled coordinate z/ZSCALE, and the interpolation kernel for the HIGHACC
# refinement pass must be sized from the *physical* Fourier f-dot.  Sizing it
# from the scaled coordinate instead leaves the kernel far too narrow, which
# truncates the acceleration response -- costing ~40% of a candidate's power
# at z=200 and ~72% at z=400, while being nearly harmless below |z| ~ 125
# (which is how PRESTO went 25 years without noticing; see issue #230).
N = 16384
predpow = N**2 / 4.0
us = np.arange(N, dtype=np.float64) / N


def ffdot_signal(r, z, numharm=1):
    """The FFT of numharm harmonics of a signal at Fourier (r, z)."""
    sig = np.zeros(N, dtype=np.float64)
    for h in range(1, numharm + 1):
        rh, zh = r * h, z * h
        sig += np.cos(2 * np.pi * (us * (us * (zh / 2.0) + (rh - 0.5 * zh))))
    return presto.rfft(sig, -1)


# Single harmonic: the truncated kernel shows up directly in the power, since
# max_rz_arr() returns the value of the function it was optimizing.
r = N / 4.0
for z in (0.0, 50.0, 200.0, 400.0):
    ft = ffdot_signal(r, z)
    for dr, dz in ((0.0, 0.0), (0.2, 1.3), (-0.35, -2.5)):
        maxpow, rmax, zmax, rd = presto.maximize_rz(ft, r + dr, z + dz, norm=1.0)
        assert maxpow > 0.99 * predpow
        assert abs(rmax - r) < 0.02
        assert abs(zmax - z) < 0.05

# Harmonic summing: maximize_rz_harmonics() recomputes each harmonic's power
# with get_derivs3d(), which sizes its own kernel, so the returned power does
# not see the truncation -- but the peak the optimizer walked to does.  A
# too-narrow kernel left |dz| of 0.02-0.04 here instead of < 0.001.
for numharm in (2, 4):
    r = N / (8.0 * numharm)
    ft = ffdot_signal(r, 100.0, numharm)
    for dr, dz in ((0.0, 0.0), (0.2, 1.3)):
        maxpow, rmax, zmax, rds = presto.maximize_rz_harmonics(
            ft, r + dr, 100.0 + dz, numharm, norm=1.0
        )
        assert maxpow > 0.99 * numharm * predpow
        assert abs(rmax - r) < 0.01
        assert abs(zmax - 100.0) < 0.01
print("success")

print("Testing angle functions...", end=" ")
dd1 = 15.25
dd2 = presto.dms2rad(*presto.deg2dms(dd1)) * presto.RADTODEG
assert round(dd1 - dd2, 12) == 0
dd1 = -0.5
dd2 = presto.dms2rad(*presto.deg2dms(dd1)) * presto.RADTODEG
assert round(dd1 - dd2, 12) == 0
hh1 = 12.125
hh2 = presto.hms2rad(*presto.hours2hms(hh1)) * presto.RADTODEG / 15.0
assert round(hh1 - hh2, 12) == 0
hh1 = -0.5
hh2 = presto.hms2rad(*presto.hours2hms(hh1)) * presto.RADTODEG / 15.0
assert round(hh1 - hh2, 12) == 0
ang = presto.sphere_ang_diff(
    10.0 * presto.DEGTORAD,
    10.0 * presto.DEGTORAD,
    50.0 * presto.DEGTORAD,
    -10.0 * presto.DEGTORAD,
)
assert round(160334.960 * presto.ARCSEC2RAD - ang, 7) == 0
print("success")

print("Testing get_baryv (barycenter)...", end=" ")
vavg1 = presto.get_baryv(
    "18:24:32.9520", "-24:52:12.0000", 56421.44222222222222, 214.5386496, obs="GB"
)
# The reference value is from the old TEMPO-based barycenter(); the
# ERFA-based version agrees with it to ~1e-9 in v/c (ephemeris difference)
vavg2 = -7.2069293455783169e-05
assert abs(vavg1 - vavg2) < 5e-9
print("success")

print("Testing simple folding code...", end=" ")
prof, phs = presto.fold(np.ones(10000), 0.001, 10, 1)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 200, 1)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 500, 1)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, tlo=-1.12313)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, startphs=0.65765)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, startphs=0.99)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
prof, phs = presto.fold(np.ones(10000), 0.001, 100, 1, startphs=0.99, standard=False)
assert np.allclose(prof, prof.mean())
assert np.all(prof > 0)
print("success")
