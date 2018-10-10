import numpy as np
from presto import presto
from presto import binary_psr
import matplotlib.pyplot as plt

N = 1000 # number of points in each orbit calc
ma = np.arange(float(N))*2.0/N

# This is going for Figure 1 in Hulse & Taylor 1975
psr1 = presto.psrepoch("B1913+16", 42320.0)
# unfortunatey, OMDOT is not in the binary part of the
# database correctly.  So we need to set that:
psr1.orb.w = 179.0
psr1.orb.t = 0.0
Eo = presto.keplers_eqn(psr1.orb.t, psr1.orb.p, psr1.orb.e, 1e-15)
Es = presto.dorbint(Eo, N, 2.0*psr1.orb.p/N, psr1.orb)
presto.E_to_v(Es, psr1.orb)
plt.plot(ma, Es, 'b-')
plt.xlabel("Orbital Phase")
plt.ylabel("Pulsar Velocity (km/s)")
plt.show()

# This is going for Figure 1 in Champion et al 2008
bpsr = binary_psr.binary_psr("1903+0327.par")
MJDs = bpsr.T0 + ma * bpsr.par.PB
xs, ys = bpsr.position(MJDs)
#cMJDs = bpsr.demodulate_TOAs(MJDs)
#cxs, cys = bpsr.position(cMJDs)

psr2 = presto.psrepoch("1903+0327.par", bpsr.T0)
psr2.orb.t = 0.0
Eo = presto.keplers_eqn(psr2.orb.t, psr2.orb.p, psr2.orb.e, 1e-15)
Es = presto.dorbint(Eo, N, 2.0*psr2.orb.p/N, psr2.orb)
# bt = Es.copy()
presto.E_to_phib(Es, psr2.orb)
# presto.E_to_phib_BT(bt, psr2.orb)

plt.plot(ma, Es, 'b-')
plt.plot(ma, -xs, 'r-')
#plt.plot(ma, Es - -xs, '-')
# plt.plot(ma, bt, 'g-')
# plt.plot(ma, -cxs, 'k-')
plt.xlabel("Orbital Phase")
plt.ylabel("Pulse Delay (s)")
plt.show()
