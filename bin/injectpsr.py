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
import scipy.interpolate
import matplotlib.pyplot as plt

import filterbank
import psr_utils

NUMSECS = 1.0 # Number of seconds of data to use to determine global scale
              # when repacking floating-point data into integers
BLOCKSIZE = 1e5 # Number of spectra to manipulate at once

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
        self.spline = None
        self.scale = scale
        self.prof_func = lambda phs: prof_func(phs % 1)

    def _get_profile(self):
        """Private method to get the pulse profile vs. phase
            function.
        """
        prof = lambda ph: self.scale*self.prof_func(ph)
        return prof

    def _get_spline(self, npts=1024, remake=False, **spline_kwargs):
        """Private method to get an interpolating spline of the profile.
        """
        if remake or self.spline is None:
            phs = np.linspace(0,1, npts+1, endpoint=True)
            prof = self._get_profile()
            self.spline = scipy.interpolate.InterpolatedUnivariateSpline(phs, \
                                                prof(phs), **spline_kwargs)
        return self.spline

    def __call__(self, phs, end_phs=None, direct=False):
        """Return the value of the profile at the given phase.

            Inputs:
                phs: The phase of the profile (between 0 and 1) where
                    the profile should be evaluated.
                end_phs: (optional) If provided integrate the profile
                    between phs and end_phs and return the result.
                    phs and end_phs should have the same shape.
                direct: If True, evaluate the profile function directly.
                    If False, use an interpolating spline. 
                    (Default: use spline)

            Output:
                vals: The values of the profile at the requested phases.
        """
        if direct:
            prof = lambda ph: self._get_profile()(ph % 1)
        else:
            prof = lambda ph: self._get_spline()(ph % 1)

        if end_phs is None:
            vals = prof(phs.flatten())
        else:
            vals = np.empty(len(phs))
            for ii, (phs0, phs1) in enumerate(zip(phs, end_phs)):
                if direct:
                    vals[ii] = scipy.integrate.quad(prof, phs0, phs1)[0]
                else:
                    vals[ii] = prof.integral(phs0, phs1)
        vals.shape = phs.shape
        return vals

    def plot(self, nbin=1024):
        x0 = np.linspace(0, 1.0, nbin, endpoint=False)
        x1 = np.linspace(0, 1.0, nbin+1, endpoint=True)[1:]
        plt.plot(x0, self(x0))


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
        self.spline = None
        self.scale = scale
        self.components = []
        for component in components:
            component.scale = self.scale # Apply the same global
                                         # scaling to all components
            self.components.append(component)

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
        comp.scale = self.scale
        self.components.append(comp)
        self.spline = None

    def plot(self, nbin=1024):
        x0 = np.linspace(0, 1.0, nbin, endpoint=False)
        plt.plot(x0, self(x0), 'k-', lw=3)
        for comp in self.components:
            comp.plot(nbin=nbin)


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
        vonmises_comps = [Profile(vonmises_factory(1.0, 5, 0.5))]
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


def inject(infile, outfn, prof, period, dm, nbitsout=8, block_size=BLOCKSIZE):
    if isinstance(infile, filterbank.FilterbankFile):
        fil = infile
    else:
        fil = filterbank.FilterbankFile(infile, read_only=True)
    print "Injecting pulsar signal into: %s" % fil.filename
    delays = psr_utils.delay_from_DM(dm, fil.frequencies)
    delays -= delays[np.argmax(fil.frequencies)]
    print np.round(delays/fil.dt).astype('int')
   
    # Create the output filterbank file
    outfil = filterbank.create_filterbank_file(outfn, fil.header, nbits=nbitsout)

    if outfil.nbits == 8:
        # Read the first second of data to get the global scaling to use
        onesec = fil.get_timeslice(0, 1).copy()
        onesec_nspec = onesec.shape[0]
        times = np.atleast_2d(np.arange(onesec_nspec)*fil.tsamp).T+delays
        phases = times/period % 1
        onesec += prof(phases)
        minimum = np.min(onesec)
        median = np.median(onesec)
        # Set median to 1/3 of dynamic range
        global_scale = (256.0/3.0) / median
        toscale = True
        del onesec
    else:
        toscale = False

    # Start an output file
    print "Creating out file: %s" % outfn
    sys.stdout.write(" %3.0f %%\r" % 0)
    sys.stdout.flush()
    nblocks = int(fil.nspec/block_size)
    remainder = fil.nspec % block_size
    oldprogress = -1
    for iblock in np.arange(nblocks):
        lobin = iblock*block_size
        hibin = (iblock+1)*block_size
        spectra = fil.get_spectra(lobin, hibin)
        times = np.atleast_2d(np.arange(lobin, hibin)*fil.tsamp).T - delays
        phases = times/period % 1
        toinject = prof(phases)
        injected = spectra+toinject
        if toscale:
            scaled = (injected-minimum)*global_scale
        else:
            scaled = injected
        outfil.append_spectra(scaled)
        progress = int(100.0*hibin/fil.nspec)
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
    # Read all remaining spectra
    if remainder:
        spectra = fil.get_spectra(-remainder, None)
        times = np.atleast_2d(np.arange(nblocks*block_size, nblocks*block_size+remainder) * \
                            fil.tsamp).T - delays
        phases = times/period % 1
        toinject = prof(phases)
        injected = spectra+toinject
        if toscale:
            scaled = (injected-minimum)*global_scale
        else:
            scaled = injected
        outfil.append_spectra(scaled)
    sys.stdout.write("Done   \n")
    sys.stdout.flush()
    

def main():
    comps = create_vonmises_components(options.vonmises)
    print "Creating profile. Number of components: %d" % len(comps)
    prof = MultiComponentProfile(comps, scale=options.scale)
    if options.dryrun:
        print "Showing plot of profile to be injected..."
        prof.plot()
        plt.xlim(0,1)
        plt.xlabel("Phase")
        plt.show()
        sys.exit()

    print "%d input files provided" % len(args)
    for fn in args:
        fil = filterbank.FilterbankFile(fn, read_only=True)
        outfn = options.outname % fil.header 
        inject(fil, outfn, prof, options.period, options.dm, \
                nbitsout=options.output_nbits, block_size=options.block_size)


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
    parser.add_option("--block-size", dest='block_size', default=BLOCKSIZE, \
                    type='float', \
                    help="Number of spectra per block. This is the amount " \
                        "of data manipulated/written at a time. (Default: " \
                        " %d spectra)" % BLOCKSIZE)
    parser.add_option("--nbits", dest='output_nbits', default=8, type=int, \
                    help="Number of bits per same to use in output " \
                        "filterbank file. (Default: 8-bits)")
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
