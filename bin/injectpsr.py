#!/usr/bin/env python

"""Inject a fake pulsar into real data, creating
a filterbank file.

Patrick Lazarus, June 26, 2012
"""
import sys
import optparse
import warnings
import shutil

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import filterbank
import psr_utils

NUMSECS = 1.0 # Number of seconds of data to use to determine global scale
              # when repacking floating-point data into integers


class Profile(object):
    """A class to represent a generic pulse profile.
    """
    def __init__(self, prof_func, scale=1):
        """Construct a profile.

            Input:
                prof_func: A function of a single variable.
                    This function should:
                        1) Represent the pulse profile.
                        2) Expect input values of phase ranging between 0 and 1.
                        3) Work when provided with a numpy array.
                scale: An overall scaling factor to multiply
                    the profile by.

            Output:
                prof: The profile object.
        """
        self.scale = scale
        self.prof_func = lambda phs: prof_func(phs % 1)*self.scale

    def _get_profile(self):
        """Private method to get the pulse profile vs. phase
            function.
        """
        return self.prof_func

    def __call__(self, phs, end_phs=None):
        """Return the value of the profile at the given phase.

            Inputs:
                phs: The phase of the profile (between 0 and 1) where
                    the profile should be evaluated.
                end_phs: (option) If provided integrate the profile
                    between phs and end_phs and return the result.
                    phs and end_phs should have the same shape.

            Output:
                vals: The values of the profile at the requested phases.
        """
        prof = lambda ph: self.scale*self._get_profile()(ph)
        if end_phs is None:
            return prof(phs)
        else:
            vals = np.empty(len(phs))
            for ii, (phs0, phs1) in enumerate(zip(phs, end_phs)):
                vals[ii] = scipy.integrate.quad(prof, phs0, phs1)[0]
            return vals

    def plot(self, nbin=1024):
        fig = plt.figure()
        ax = plt.axes()
        x0 = np.linspace(0, 1.0, nbin, endpoint=False)
        x1 = np.linspace(0, 1.0, nbin+1, endpoint=True)[1:]
        plt.plot((x0+x1)/2.0, self(x0, end_phs=x1), 'k-', lw=3)
        plt.xlim(0,1)
        plt.xlabel("Phase")


class MultiComponentProfile(Profile):
    """A class to represent a pulse profile made up of 
        multiple components.
    """
    def __init__(self, components=None, scale=1):
        """Construct a multi-component profile.

            Input:
                components: A list of Profile object that serve
                    as the components of this MultiComponentProfile 
                    object. (Default: Create a multi-component profile
                    with no components.)
                scale: An overall scaling factor to multiply 
                    the profile by.

            Output:
                prof: The MultiComponentProfile object.
        """
        self.scale = scale
        if components is None:
            self.components = []
        else:
            self.components = components

    def _get_profile(self):
        """Private method to get the pulse profile vs. phase
            function.
        """
        if self.components:
            prof = lambda ph: np.sum([comp._get_profile()(ph) for comp \
                                        in self.components], axis=0)
        else:
            prof = lambda ph: 0
        return prof

    def add_component(self, comp):
        self.components.append(comp)

    def plot(self, nbin=1024):
        fig = plt.figure()
        ax = plt.axes()
        x0 = np.linspace(0, 1.0, nbin, endpoint=False)
        x1 = np.linspace(0, 1.0, nbin+1, endpoint=True)[1:]
        plt.plot((x0+x1)/2.0, self(x0, end_phs=x1), 'k-', lw=3)
        for comp in self.components:
            plt.plot((x0+x1)/2.0, comp(x0, end_phs=x1), lw=2)
        plt.xlim(0,1)
        plt.xlabel("Phase")


def vonmises_factory(amp,shape,loc):
    # Need to use a factory for the von Mises functions
    # to make sure the lambda uses amp,shape,loc from a local
    # scope. The values in a lambda function are stored by reference
    # and only looked up dynamically when the function is called.
    def vm(ph): 
        return amp*np.exp(shape*(np.cos(2*np.pi*(ph-loc))-1))
    return vm


def create_vonmises_components(vonmises_strs):
    if not vonmises_strs:
        warnings.warn("Using default von Mises profile (Amplitude=1.0 " \
                        "b=5, and phase=0.5)")
        vonmises = [(1.0, 5, 0.5)]
    else:
        vonmises_comps = []
        for vonmises_str in vonmises_strs:
            split = vonmises_str.split()
            if len(split) != 3:
                raise ValueError("Bad number of von Mises components " \
                        "should be 3, got %d" % len(split))
            amp = float(split[0])
            shape = float(split[1])
            loc = float(split[2])

            # Need to use a factory for the von Mises functions
            # to make sure the lambda uses amp,shape,loc from a local
            # scope. The values in a lambda function are stored by reference
            # and only looked up dynamically when the function is called.
            vonmises_comps.append(Profile(vonmises_factory(amp,shape,loc)))
    return vonmises_comps


def inject(fn, outfn, prof, period, dm):
    print "Injecting pulsar signal into: %s" % fn
    fil = filterbank.FilterbankFile(fn, read_only=True)
    nbin = int(np.round(period/fil.tsamp)) # Number of bins 
                                                   # across the profile
    # Because of manipulations performed on the arrays we need
    # to flip the delays to get a proper DM sweep.
    delays = psr_utils.delay_from_DM(dm, fil.frequencies)[::-1]
    delays_phs = delays/period
    delays_phs -= delays_phs[np.argmax(fil.frequencies)]
    delays_phs %= 1
    
    # for ii, (freq, delay) in enumerate(zip(fil.frequencies, delays_phs)):
    #     print "Channel %d: %g MHz, %g (phase)" % (ii, freq, delay)
   
    # Copy the input file to the output file name and modify the output
    # file in-place
    shutil.copy(fn, outfn) 

    # Read the first second of data to get the global scaling to use
    outfil = filterbank.FilterbankFile(outfn, read_only=False)
    onesec = outfil.get_timeslice(0, 1).copy()
    onesec_nspec = onesec.shape[0]
    phases = np.tile(np.arange(onesec_nspec)*fil.tsamp/period % 1, \
                        (fil.nchans,1)).T + delays_phs
    onesec += prof(phases)
    minimum = np.min(onesec)
    median = np.median(onesec)
    # Set median to 1/3 of dynamic range
    global_scale = (256.0/3.0) / median
    del onesec
    
    # Start an output file
    print "Creating out file: %s" % outfn
    nprofs = outfil.nspec/nbin
    remainder = outfil.nspec % nbin
    oldprogress = -1
    for iprof in np.arange(nprofs):
        spectra = outfil.get_spectra(iprof*nbin, (iprof+1)*nbin)
        phases = np.tile(np.arange(iprof*nbin, (iprof+1)*nbin)*fil.tsamp/period % 1, (fil.nchans,1)).T + delays_phs
        toinject = prof(phases)
        outfil.spectra[iprof*nbin:(iprof+1)*nbin] = \
                np.clip((spectra+toinject-minimum)*global_scale, 0, 256)
        outfil.spectra.flush()
        progress = int(100.0*((iprof+1)*nbin)/outfil.nspec)
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
    # Read all remaining spectra
    if remainder:
        spectra = outfil.get_spectra(-remainder, None)
        phases = np.tile(np.arange(nprofs*nbin, nprofs*nbin+remainder)*fil.tsamp/period % 1, (fil.nchans,1)).T + delays_phs
        toinject = prof(phases)
        outfil.spectra[-remainder:] = \
                np.clip((spectra+toinject-minimum)*global_scale, 0, 256)
        outfil.spectra.flush()
    sys.stdout.write("Done   \n")
    sys.stdout.flush()
    

def main():
    comps = create_vonmises_components(options.vonmises)
    prof = MultiComponentProfile(comps, scale=options.scale)
    if options.dryrun:
        print "Showing plot of profile to be injected..."
        prof.plot()
        plt.show()
        sys.exit()

    print "%d input files provided" % len(args)
    for fn in args:
        outfn = options.outname % fil.header 
        inject(fn, outfn, prof, options.period, options.dm)


def parse_model_file(modelfn):
    """Parse a pass model file (*.m) written by paas.
        Return a list of parameters describing each component.
        In particular (amplitude, shape, phase).

        Input:
            modelfn: The name of the model file.

        Outputs:
            params: List of parameters for each component.
                (i.e. "amplitude shape phase")
    """
    mfile = open(modelfn, 'r')
    return [" ".join(reversed(line.split())) \
                        for line in mfile.readlines()]


def parse_mfile_callback(option, opt_str, value, parser):
    vonmises = getattr(parser.values, 'vonmises')
    vonmises.extend(parse_model_file(value))


if __name__ == '__main__':
    parser = optparse.OptionParser(prog='injectpsr.py', \
                    version="v0.1 Patrick Lazarus (June 26, 2012)")
    parser.add_option("--dm", dest='dm', action='store', type='float', \
                    help="The DM of the (fake) injected pulsar signal. " \
                        "(This argument is required.", \
                    default=None)
    parser.add_option("-p", "--period", dest='period', action='store', \
                    default=None, type='float', \
                    help="The period (in seconds) of the (fake) injected " \
                        "pulsar signal. (This argument is required.)")
    parser.add_option("-s", "--scale", dest='scale', type='float', \
                    default=1, \
                    help="Overall scaling factor to multiply profile with. " \
                        "(Default: Don't scale.)")
    parser.add_option("-v", "--vonmises", dest='vonmises', action='append', \
                    help="A string of 3 parameters defining a vonmises " \
                        "component to be injected. Be sure to quote the " \
                        "3 parameters together. The params are: 'amplitude " \
                        "shape phase'. Amplitude is not related to SNR in " \
                        "any way. Also, phase should be between 0 and 1. " \
                        "(Default: if no compoments are provided " \
                        "a von Mises with amplitude=1.0, shape=5, and " \
                        "phase=0.5 will be used.)", \
                    default=[])
    parser.add_option("-m", "--model-file", dest="model_file", nargs=1, 
                    type='str', \
                    action="callback", callback=parse_mfile_callback, \
                    help="A model file (*.m) as written by 'paas'.")
    parser.add_option("-n", "--dryrun", dest="dryrun", action="store_true", \
                    help="Show the pulse profile to be injected and exit. " \
                        "(Default: do not show profile, inject it)")
    parser.add_option("-o", "--outname", dest='outname', action='store', \
                    default="injected.fil", \
                    help="The name of the output file.")
    (options, args) = parser.parse_args()
    if options.period is None or options.dm is None:
        raise ValueError("Both a period and a DM _must_ be provided!")
    main()
