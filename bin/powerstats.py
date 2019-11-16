#!/usr/bin/env python
from __future__ import print_function
from builtins import input
from presto.events import *


def answer_yes(question):
    yes = ['', 'Y', 'y', 'Yes', 'yes', 'YES',
           'T', 't', 'True', 'true', 'TRUE']
    return input('\n'+question) in yes


def ask_float(question, default=None):
    while 1:
        ans = input('\n'+question)
        if not ans:
            ans = default
        try:
            return float(ans)
        except (ValueError, TypeError):
            print("\nThat was not a valid number.  Try again...\n")

def ask_int(question, default=None):
    while 1:
        ans = input('\n'+question)
        if not ans:
            ans = default
        try:
            return int(ans)
        except (ValueError, TypeError):
            print("\nThat was not a valid number.  Try again...\n")

if __name__ == '__main__':
    print("\nPower Statistics Calculation Routine")

    conf = ask_float(\
        "What confidence level would you like to use?  [0.99]  ", 0.99)
    Ntot = ask_int(\
        "How many data points were FFTd (N)?  ")
    dt = ask_float("What was the length in time (s) of each bin?  ")
    T = Ntot * dt
    P_max = ask_float(\
        "What was the maximum normalized power found?  ")
    rlo = 1
    rhi = Ntot / 2
    if answer_yes(\
        "Was this an RZW acceleration search (y/n)?  [y]  "):
        rlo = T * ask_float(\
                "What was the lowest freq searched (Hz)?   [1.0]  ", 1.0)
        rhi = T * ask_float(\
                "What was the highest freq searched (Hz)? [%.2f]  " %
                ((Ntot/2.0)/T), (Ntot/2.0)/T)
        zlo = ask_float(\
                "What was the lowest 'z' value searched?  [-100]  ", -100.0)
        zhi = ask_float(\
                "What was the highest 'z' value searched?  [100]  ", 100.0)
        Nsearch = (rhi - rlo) * (zhi - zlo + 1.0) / 6.95
    else:
        Nsearch = ask_int(\
                "How many independent bins were searched?  [N/2]  ", Ntot/2)
    if answer_yes(\
        "Was the data composed of binned counts (y/n)?  [y]  "):
        numphot = ask_int("How many counts (photons) were there?  ")
        lofreq, hifreq = rlo / T, rhi / T
        trial_freqs = (10.0**(Num.arange(7.0)-2.0)).tolist()
        trial_freqs = [x for x in trial_freqs if x > lofreq and x < hifreq]
        print("\nThe trial frequencies (Hz) are:", trial_freqs)
        if answer_yes(\
            "Would you like to add any more?  [y]  "):
            new_freq = ask_float(\
                "Enter a frequency (Hz) or '0' to stop.  ")
            while (new_freq):
                trial_freqs.append(new_freq)                
                new_freq = ask_float(\
                    "Enter a frequency (Hz) or '0' to stop.  ")
        trial_freqs.sort()
        print("\n\nCalculating...\n\n")
        print("")
        print("         Power Stats for Binned Data")
        print("     -----------------------------------")
        print("      Number of data points = %.0f" % Ntot)
        print("        Time per sample (s) = %g" % dt)
        print("    Total number of photons = %.0f" % numphot)
        print("           Confidence Level = %g%%" % (100 * conf))
        print(" Number of independent bins = %.2e" % Nsearch)
        print(" Threshold Power (P_detect) > %.2f" % \
              max_noise_power(Nsearch, 1, conf))
        ulim = required_signal_power(P_max, 1, conf)
        print("    Max Power Found (P_max) = %.2f" % P_max)
        print(" Max Signal Power (P_limit) < %.2f" % ulim)
        print("  Pulsed Fraction (P_limit) < %.3g" % \
              pulsed_fraction_limit(numphot, ulim))
        print("")
        sens = []
        ulim = []
        for f in trial_freqs:
            sens.append(binned_fft_sensitivity(Ntot, dt, f, Nsearch, 1, conf))
            ulim.append(required_signal_power(P_max, 1, conf))
        print("          Freq (Hz)  = ", end=' ')
        for f in trial_freqs:
            print(" f=%-7g" % (f), end=' ')
        print('\n                       '+'-'*len(trial_freqs)*11)
        print("  Power Sensitivity  > ", end=' ')
        for s in sens:
            print(" %-8.2f " % (s), end=' ')
        print('')
        pfract = []
        for s in sens:
            pfract.append(pulsed_fraction_limit(numphot, s))
        print("    Pulsed Fraction  < ", end=' ')
        for p in pfract:
            print(" %-8.3g " % (p), end=' ')
        print('\n')
    else:
        print("\n\nCalculating...\n\n")
        print("")
        print("         Power Stats for Normal Data")
        print("     -----------------------------------")
        print("      Number of data points = %.0f" % Ntot)
        print("           Confidence Level = %g%%" % (100 * conf))
        print(" Number of independent bins = %.2e" % Nsearch)
        print(" Threshold Power (P_detect) > %.2f" % \
              max_noise_power(Nsearch/2, 1, conf))
        sens = fft_sensitivity(Ntot, Nsearch, 1, conf)
        print("          Power Sensitivity > %.2f" % sens)
        ulim = required_signal_power(P_max, 1, conf)
        print("    Max Power Found (P_max) = %.2f" % P_max)
        print(" Max Signal Power (P_limit) < %.2f" % ulim)
        print("")
