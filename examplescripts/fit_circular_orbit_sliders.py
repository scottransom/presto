#!/usr/bin/env python
import sys
import numpy as np
from presto import psr_utils as pu
from presto import psr_constants as pc
from presto import parfile
from presto import bestprof
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.optimize import leastsq

period = np.asarray([])
time = np.asarray([])

def parse_eph(filenm):
    global period, time
    suffix = filenm.split(".")[-1]
    if suffix=="bestprof":
        x = bestprof.bestprof(filenm)
        fs = pu.p_to_f(x.p0_bary, x.p1_bary, x.p2_bary)
        epoch = x.epochi_bary + x.epochf_bary
        T = x.T
    elif suffix=="par":
        x = parfile.psr_par(filenm)
        # Try to see how many freq derivs we have
        fs = [x.F0]
        for ii in range(1, 20):  # hopefully 20 is an upper limit!
            attrib = "F%d"%ii
            if hasattr(x, attrib):
                fs.append(getattr(x, attrib))
            else:
                break
        epoch = x.PEPOCH
        T = (x.FINISH - x.START) * 86400.0
    else:
        print("I don't recognize the file type for", filenm)
        sys.exit()
    newts = epoch + np.arange(int(T/10.0+0.5), dtype=float)/8640.0
    time = np.concatenate((time, newts))
    newps = 1.0 / pu.calc_freq(newts, epoch, *fs)
    period = np.concatenate((period, newps))
    print("%13.7f (%0.1f sec): " % (epoch, T), fs)


def orbeqn(Ppxt, times):
    # P = Ppsr, p = Porb, x = a*sin(i)/s, t = T_o
    PBsec = Ppxt[1] * 86400.0
    phi = pc.TWOPI*(times - Ppxt[3])*86400.0/PBsec
    return Ppxt[0]*(1.0+pc.TWOPI*Ppxt[2]/PBsec*np.cos(phi))


def funct(Ppxt, times, measured):
    return orbeqn(Ppxt, times) - measured


# The function to be called anytime a slider's value changes
def update(val):
    model.set_ydata(orbeqn([Ppsr_slider.val, Porb_slider.val,
                            Xorb_slider.val, T0_slider.val], model_time) - Ppsr_slider.val)
    data.set_ydata(period - Ppsr_slider.val)
    ax.set_ylabel("Pulsar Period - %.12f (s)"%(Ppsr_slider.val))
    fig.canvas.draw_idle()


if __name__ == '__main__':
    if len(sys.argv)==1:
        print("\nusage: fit_circular_orbit.py Porb Xorb T0 parfiles or bestprofs")
        exit(0)
    initPorb = float(sys.argv[1])
    initXorb = float(sys.argv[2])
    initT0 = float(sys.argv[3])
    # This fills the time and period arrays
    for infile in sys.argv[4:]:
        parse_eph(infile)
    period = np.asarray(period, dtype=float)
    initPpsr = period.mean()
    time = np.asarray(time, dtype=float)

    T = max(time)-min(time)
    model_time = np.arange(min(time)-0.1*T, max(time)+0.1*T, 0.01)

    fig, ax = plt.subplots()
    Ppsr, Porb, Xorb, T0 = initPpsr, initPorb, initXorb, initT0
    data, = ax.plot(time-model_time[0], period-initPpsr, '.')
    model, = ax.plot(model_time-model_time[0], orbeqn([Ppsr, Porb, Xorb, T0], model_time)-Ppsr, 'r')
    ax.set_xlabel("Days + %.7f"%model_time[0])
    ax.set_ylabel("Pulsar Period - %.12f (s)"%Ppsr)

    # adjust the main plot to make room for the sliders
    fig.subplots_adjust(left=0.25, bottom=0.25)

    # Make a horizontal slider to control the orbital period.
    axPorb = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    Porb_slider = Slider(
        ax=axPorb,
        label='Porb in days',
        valmin=initPorb / 3.0,
        valmax=initPorb * 3.0,
        valinit=initPorb,
    )
    # Make a horizontal slider to control the phase.
    axT0 = fig.add_axes([0.25, 0.05, 0.65, 0.03])
    T0_slider = Slider(
        ax=axT0,
        label='T0 in MJD',
        valmin=initT0 - initPorb,
        valmax=initT0 + initPorb,
        valinit=initT0,
    )

    # Make a vertically oriented slider to control the spin period
    axPpsr = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
    Ppsr_slider = Slider(
        ax=axPpsr,
        label="Spin Period in s",
        valmin=0.999 * initPpsr,
        valmax=1.001 * initPpsr,
        valinit=initPpsr,
        orientation="vertical"
    )
    # Make a vertically oriented slider to control the semi-major axis
    axXorb = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
    Xorb_slider = Slider(
        ax=axXorb,
        label="asini/c in lt-sec",
        valmin=initXorb / 4.0,
        valmax=initXorb * 4.0,
        valinit=initXorb,
        orientation="vertical"
    )
    # register the update function with each slider
    Ppsr_slider.on_changed(update)
    Porb_slider.on_changed(update)
    Xorb_slider.on_changed(update)
    T0_slider.on_changed(update)

    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')

    def reset(event):
        Ppsr_slider.reset()
        Porb_slider.reset()
        Xorb_slider.reset()
        T0_slider.reset()
    button.on_clicked(reset)

    plt.show()

    #ret = leastsq(funct, [Ppsr, Porb, Xorb, T0], args=(time, period))
    #To = ret[0][3]

    #if (ret[0][2] < 0.0):
    #    print("Modifying TO because of negative asini/c...")
    #    ret[0][3] += 0.5 * (ret[0][1]/86400.0)
    #    ret[0][2] = abs(ret[0][2])

    #print("P_orb = %.3f hrs" % (ret[0][1]/3600.0))
    #print("P0  %17.15g 1" % ret[0][0])
    #print("PB  %17.15g 1" % (ret[0][1]/86400.0))
    #print("A1  %17.15g 1" % ret[0][2])
    #print("T0  %17.15g 1" % ret[0][3])
    #print("E                 0.0")
    #print("OM                0.0")




