#!/usr/bin/env python

"""Inject a fake pulsar into real data, creating
a filterbank file.

Patrick Lazarus, June 26, 2012
"""
import sys
import optparse
import warnings
import pickle
import copy

import numpy as np
import scipy.integrate
import scipy.interpolate
import matplotlib
matplotlib.use('agg') # Use a non-interactive backend
import matplotlib.pyplot as plt
import scipy.integrate

import filterbank
import psr_utils

NUMSECS = 1.0 # Number of seconds of data to use to determine global scale
              # when repacking floating-point data into integers
BLOCKSIZE = 1e4 # Number of spectra to manipulate at once


class Profile(object):
    """A class to represent a generic pulse profile.
    """
    def __init__(self, prof_func, scale=1):
        """Construct a profile.

            Inputs:
                prof_func: A function of a single variable.
                    This function should:
                        1) Represent the pulse profile.
                        2) Expect input values of phase ranging between 
                            0 and 1.
                        3) Work when provided with a numpy array.
                scale: An overall scaling factor to multiply
                    the profile by.

            Output:
                prof: The profile object.
        """
        self.prof_func = prof_func
        self.scale = scale

    def __call__(self, phs):
        """Return the value of the profile at the given phase.

            Inputs:
                phs: The phase of the profile (between 0 and 1) where
                    the profile should be evaluated.

            Output:
                vals: The values of the profile at the requested phases.
        """
        return self.scale*self.prof_func(phs)

    def plot(self, nbin=1024, scale=1):
        x0 = np.linspace(0, 1.0, nbin+1, endpoint=True)
        plt.plot(x0, self(x0)*scale)
        plt.xlabel("Phase")

    def get_area(self, npts=1024):
        """Return the area under the pulse in units of (intensity x phase).
            The area is calculated by evaluating the profile at many points
            and numerically integrated using the trapezoid rule.

            NOTE: the scale-factor will be applied.

            Input:
                npts: The number of points to use when evaluating the
                    profile.

            Ouput:
                area: The area under the pulse in units of (intensity x phase).
        """
#        phs = np.linspace(0, 1.0, npts+1, endpoint=True)
#        area = np.trapz(y=self(phs), x=phs)
        area, err = scipy.integrate.quadrature(self, 0, 1, maxiter=250)
        return area

    def get_max(self, npts=1024):
        """Return the maximum value of the profile.
            The profile is evaluated at many points. The quantity returned
            is the maximum value evaluated.

            NOTE: the scale-factor will be applied.

            Inputs:
                npts: The number of points to use when evaluating the
                    profile.

            Ouput:
                profmax: The profile maximum.
        """
        phs = np.linspace(0, 1.0, npts+1, endpoint=True)
        profmax = np.max(self(phs))
        return profmax
        
    def get_equivalent_width(self, npts=1024):
        """Determine and return the equivalent width of the profile, in phase.
            The equivalent width is the area under the pulse divided
            by the profile's maximum value.

            Input:
                npts: The number of points to use when evaluating the
                    profile.

            Ouput:
                weq: The equivalent width of the profile, in phase.
        """
        return self.get_area(npts=npts)/self.get_max(npts=npts)

    def set_scale(self, scale):
        """Set the profile's scaling factor.

            Input:
                scale: The scaling factor to use.

            Outputs:
                None
        """
        self.scale = scale

    def delay(self, phasedelay):
        """Delay the profile and return a new Profile object. 
 
            Input:
                phasedelay: The amount of phase to delay the profile by.

            Output:
                delayed: The delayed Profile.
        """
        delayed_prof = Profile(lambda ph: self((ph-phasedelay) % 1)/self.scale, \
                                scale=self.scale)
        return delayed_prof
    
    def convolve_with(self, other, npts=4096):
        """Convolve Profile with another. Return a SplineProfile
            with the requested number of points.

            Inputs:
                other: The Profile to convolve with.
                npts: The number of points to use when creating the
                    resulting convolution (i.e. a SplineProfile). 
                    (Default: 1024)

            Other:
                convolution: The convolution, a SplineProfile object.
        """
        phs = np.linspace(0, 1, npts, endpoint=False)
        conv_vals = np.fft.ifft(np.fft.fft(self(phs))*np.fft.fft(other(phs)))
#        window = np.linspace(0, 1, 4096, endpoint=False)
#        conv_vals = np.empty(npts)       
#        for ii, ph in enumerate(phs):
#            conv_vals[ii] = np.trapz(self(window)*other((ph-window)%1), x=window)
##            integrand = lambda xx: self(xx)*other((ph-xx)%1)
##            conv_vals[ii], err = scipy.integrate.quadrature(integrand, 0, 1, maxiter=200)
        convolution = SplineProfile(conv_vals/self.scale, scale=self.scale)
        return convolution

    def smear(self, smearphs, delayphs=0, npts=4096):
        """Smear the profile with a boxcar of width 'smearphs'. Return
            a SplineProfile object sampled with 'npts' points.

            Inputs:
                smearphs: The amount (in phase) to smear the profile by.
                delayphs: The amount (in phase) to delay the pulse by.
                    (Default: No delay)
                npts: The number of points to use when creating the
                    smeared SplineProfile. (Default: 4096)

            Output:
                smeared: The smeared Profile.
        """
        if smearphs < 0:
            raise ValueError("Amount of phase to smear by (%g) " \
                                "cannot be negative!" % scatterphs)
        elif smearphs == 0:
            # No scattering to do, return profile
            return copy.deepcopy(self) 
        bc = boxcar_factory(smearphs, delayphs)
        return self.convolve_with(bc, npts=npts)

    def scatter(self, scatterphs, npts=4096):
        """Scatter the profile with a one-sided exponential of width
            'scatterphs'. Return a SplineProfile object sampled
            with 'npts' points.

            Inputs:
                scatterphs: The time-scale (in phase) of the exponential
                    scattering function.
                npts: The number of points to use when creating the scattered
                    SplineProfile. (Default: 4096)

            Outputs:
                scattered: The scattered Profile.
        """
        if scatterphs < 0:
            raise ValueError("Amount of phase to scatter by (%g) " \
                                "cannot be negative!" % scatterphs)
        elif scatterphs == 0:
            # No scattering to do, return profile
            return copy.deepcopy(self)
        ex = exponential_factory(scatterphs)
        return self.convolve_with(ex, npts=npts)

class SplineProfile(Profile):
    def __init__(self, profvals, scale=1, **spline_kwargs):
        """Construct a profile that uses a spline to interpolate a function.

            Inputs:
                profvals: The values of the profile to be interpolated. 
                scale: An overall scaling factor to multiply
                    the profile by.
                **All additional keyword arguments are passed to the 
                    spline constructor.

            Output:
                prof: The profile object.
        """
        # TODO: Should we evaluate at the centre of the bins?
        phs = np.linspace(0,1, len(profvals)+1, endpoint=True)
        # Manually set value at phs=1.0 to the value at phs=0.0
        vals = np.concatenate((profvals, [profvals[0]]))
        # Create spline object and use it as the profile function
        spline = scipy.interpolate.InterpolatedUnivariateSpline(phs, \
                                                vals, **spline_kwargs)
        super(SplineProfile, self).__init__(spline, scale)

    def __call__(self, phs):
        """Return the value of the profile at the given phase.

            Inputs:
                phs: The phase of the profile (between 0 and 1) where
                    the profile should be evaluated.

            Output:
                vals: The values of the profile at the requested phases.
        """
        vals = super(SplineProfile, self).__call__(phs.flat)
        # Re-shape values because spline return flattened array.
        vals = np.atleast_1d(vals)
        vals.shape = phs.shape
        return vals


class MultiComponentProfile(Profile):
    """A class to represent a pulse profile made up of 
        multiple components.
    """
    def __init__(self, components=None, scale=1):
        """Construct a multi-component profile.

            Input:
                components: A list of Profile objects that serve
                    as the components of this MultiComponentProfile 
                    object. (Default: Create a multi-component profile
                    with no components.)
                scale: An overall scaling factor to multiply 
                    the profile by.

            Output:
                prof: The MultiComponentProfile object.
        """
        self.scale = scale
        self.components = []
        for component in components:
            self.add_component(component)
        super(MultiComponentProfile, self).__init__(self._get_profile(), scale)

    def _get_profile(self):
        """Private method to get the pulse profile vs. phase
            function.
        """
        if self.components:
            prof = lambda ph: np.sum([comp(ph) for comp \
                                        in self.components], axis=0)
        else:
            prof = lambda ph: 0
        return prof

    def add_component(self, comp):
        self.components.append(comp)

    def plot(self, nbin=1024):
        super(MultiComponentProfile, self).plot(nbin=nbin)
        for comp in self.components:
            comp.plot(nbin=nbin, scale=self.scale)



class VectorProfile(object):
    """A class to represent a vector of pulse profiles.
        This can be used to encode intrinsic profile variation, 
        or extrisinc smearing/scattering across the band
    """
    def __init__(self, profiles, scale=1):
        """Construct a vector of profiles.
            
            Inputs:
                profiles: A list of Profile objects.
                scale: An overall scaling factor to multiply
                    the profile by.

            Output:
                prof: The profile object.
        """
        self.profiles = profiles
        self.nprofs = len(profiles)
        self.scale = scale
    
    def __call__(self, phs):
        phs = np.atleast_1d(np.asarray(phs))
        nphs = phs.shape[0]
        vals = np.empty((nphs, self.nprofs))
        if phs.ndim == 1:
            # Evaluate all profiles at the same phases
            for ii, prof in enumerate(self.profiles):
                vals[:,ii] = prof(phs)*self.scale
        elif phs.ndim == 2:
            # Evaluate each profile at a different set of phases
            nphs_vecs = phs.shape[1]
            if nphs_vecs != self.nprofs:
                raise ValueError("Length of axis=1 of 'phs' (%d) must be " \
                                "equal to the number of profiles in the " \
                                "vector (%d)." % (nphs_vecs, self.nprofs))
            else:
                for ii, (prof, ph) in enumerate(zip(self.profiles, phs)):
                    vals[:,ii] = prof(ph)*self.scale
        else:
            raise ValueError("VectorProfile can only be evaluated with " \
                            "1D or 2D arrays")
        return vals

    def plot(self, nbin=1024, scale=1):
        phs = np.linspace(0, 1.0, nbin+1, endpoint=True)
        plt.imshow(self(phs).transpose(), interpolation='nearest', \
                    extent=(0, 1, 0, self.nprofs), aspect='auto')
        plt.xlabel("Phase")
        plt.ylabel("Channel number")


def apply_dm(inprof, period, dm, chan_width, freqs, tsamp, \
                do_smear=True, do_scatter=True):
    """Given a profile apply DM delays, smearing, and scattering 
        within each channel as is appropriate for the given params.

        Inputs:
            inprof: The profile to modify.
            period: The profiles period (in seconds)
            dm: The DM (in pc cm-3)
            chan_width: The width of each channel (in MHz)
            freqs: The list of frequencies (in MHz)
            tsamp: Sample time of the recipient filterbank file (in seconds).
            do_smear: Boolean, if True apply DM smearing to each channel.
                (Default: True)
            do_scatter: Boolean, if True apply scattering to each channel.
                (Default: True)

        Outputs:
            vecprof: The delayed and smeared VectorProfile.
    """
    weq = inprof.get_equivalent_width()
    nfreqs = len(freqs)
    print "Applying DM to profile (DM = %.2f; %d channels)..." % \
                (dm, nfreqs)
    # A list of profiles, one for each channel
    profiles = []

    # Prepare delays
    timedelays = psr_utils.delay_from_DM(dm, freqs)
    # Reference all delays to highest frequency channel, which remains
    # unchanged
    # TODO: Do we really want to refer to high freq?
    timedelays -= timedelays[np.argmax(freqs)]
    phasedelays = timedelays/period

    # Prepare for smear campaign
    smeartimes = psr_utils.dm_smear(dm, chan_width, freqs) # In seconds
    smearphases = smeartimes/period
    
    # Prepare to scatter
    scattertimes = psr_utils.pulse_broadening(dm, freqs)*1e-3 # In seconds
    scatterphases = scattertimes/period

    oldprogress = 0
    sys.stdout.write(" %3.0f %%\r" % oldprogress)
    sys.stdout.flush()
    for ii, (delayphs, smearphs, scattphs) in \
                enumerate(zip(phasedelays, smearphases, scatterphases)):
        #########
        # DEBUG: plot all profiles
        plt.clf()
        ax=plt.subplot(5,1,1)
        inprof.plot()
        if do_smear and not ((smearphs < 0.2*weq) or (smearphs < (tsamp/period))):
            # Only smear if requested and smearing-phase is large enough
            bc = boxcar_factory(smearphs, delayphs)
            plt.subplot(5,1,2,sharex=ax)
            bc.plot()
            tmpprof = inprof.smear(smearphs, delayphs)
        else:
            print "Not smearing"
            tmpprof = inprof.delay(delayphs)
        plt.subplot(5,1,3,sharex=ax)
        tmpprof.plot()
        if do_scatter and not ((scattphs < 0.2*weq) or (scattphs < (tsamp/period))):
            # Only scatter if requested and scattering-phase is large enough
            ex = exponential_factory(scattphs)
            plt.subplot(5,1,4,sharex=ax)
            ex.plot()
            tmpprof = tmpprof.scatter(scattphs)
        else:
            print "Not scattering"
        plt.subplot(5,1,5,sharex=ax)
        tmpprof.plot()
        profiles.append(tmpprof)
        plt.xlim(0,1)
        plt.xlabel("Phase")
        plt.suptitle("Prof %d (%f MHz)" % (ii, freqs[ii]))
        plt.savefig("prof%d.png" % ii)
        #########
        # Print progress to screen
        progress = int(100.0*ii/nfreqs)
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
    sys.stdout.write("Done   \n")
    sys.stdout.flush()
    vecprof = VectorProfile(profiles, scale=inprof.scale)
    return vecprof


def get_spline_profile(prof, npts=1024, **spline_kwargs):
    """Given a profile object evaluate it and return
        a SplineProfile object. If the input profile object
        is already an instance of SplineProfile, do nothing
        and return the input profile.

        Inputs:
            prof: The profile object to conver to a SplineProfile.
            npts: The number of points to use when evaluating the
                profile. (Default: 1024)
            **All additional keyword arguments are passed to the 
                spline constructor.

        Outputs:
            spline_prof: The resulting SplineProfile object.
    """
    if isinstance(prof, SplineProfile):
        # Input profile is already a SplineProfile. Do nothing. Return it.
        return prof
    else:
        phs = np.linspace(0,1, npts, endpoint=False)
        profvals = prof(phs)/prof.scale
        spline_prof = SplineProfile(profvals, scale=prof.scale, **spline_kwargs)
        return spline_prof


def vonmises_factory(amp,shape,loc):
    """Return a Profile with a single von Mises component.

        Inputs:
            amp: Amplitude of the von Mises function.
            shape: A parameter characterising the width 
                (NOTE: this value is not directly the width.)
            loc: Phase of the peak of the von Mises function.

        Output:
            vm_prof: A Profile object with a von Mises profile described
                by the input parameters.
    """
    # Need to use a factory for the von Mises functions
    # to make sure the lambda uses amp,shape,loc from a local
    # scope. The values in a lambda function are stored by reference
    # and only looked up dynamically when the function is called.
    def vm(ph): 
        return amp*np.exp(shape*(np.cos(2*np.pi*(ph-loc))-1))
    return Profile(vm)


def boxcar_factory(width, delay=0):
    """Return a boxcar Profile scaled to have unit area.

        Inputs:
            width: The width of the boxcar function in phase.
                NOTE: if width > 1, it will be folded
            delay: The delay, in phase, to apply to the boxcar.
                NOTE: positive delays cause the boxcar to be shifted
                (Default: No delay)
                to the right (i.e. larger pulse phase)

        Output:
            boxcar_prof: A boxcar Profile object with the given width.
    """
    width = float(width) # Make sure width is a floating-point number
    delay = delay % 1 # Make sure delay is in [0, 1)
    nwraps = int(width/1)
    rem = ((width%1)+delay)%1
    if delay < rem:
        # Boxcar does not straddle phase=1
        def bc(ph):    
            return (nwraps + ((delay<=ph) & (ph<rem)).astype('int'))/width
    else:
        # Boxcar is split
        def bc(ph):
            return (nwraps + ((ph<rem) | (delay<=ph)).astype('int'))/width
    return Profile(bc)


def exponential_factory(efold):
    """Return a one-sided exponential Profile.

        Inputs:
            efold: The e-folding phase of the exponential function.

        Output:
            exp_prof: A one-sided exponential Profile object.
    """
    denom = (1-np.exp(-1/efold))*efold
    def ex(ph):
        # Denominator comes from sum of geometric series
        return np.exp(-ph/efold)/denom
    return Profile(ex)


def create_vonmises_components(vonmises_strs):
    if not vonmises_strs:
        warnings.warn("Using default von Mises profile (Amplitude=1.0 " \
                        "b=5, and phase=0.5)")
        vonmises_comps = [vonmises_factory(1.0, 5, 0.5)]
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
            vonmises_comps.append(vonmises_factory(amp,shape,loc))
    return vonmises_comps


def determine_and_set_scale(fil, prof, smean, gain, tsys, rms):
    """Set the profile's scaling factor such that the simulated 
        injected pulsar signal will have the given Smean.

        Inputs:
            fil: A FilterbankFile object.
            prof: The Profile object representing the profile to inject.
            smean: The mean flux density to simulate, in mJy.
            gain: The telescope's gain, in K/Jy.
            tsys: The observing system's temperature, in K.
            rms: The RMS of the recipient file's DM=0 time series.

        Outputs:
            None
    """
    # Set scale to 1, just in case it's be altered already.
    prof.set_scale(1)
    
    # Characterise the pulse
    area = prof.get_area()
    profmax = prof.get_max()
    width = prof.get_equivalent_width() # in phase
    dutycycle = width # Because width is reported in phase, 
                      # it is actually the duty cycle

    # Characterise the recipient filterbank file
    tint = fil.nspec*fil.tsamp
    bw = np.abs(fil.foff*fil.nchans)

    # Target SNR
    warnings.warn("Assuming 2 (summed) polarizations.")
    snr = smean*gain*np.sqrt(2*tint*bw)/tsys*np.sqrt(1/dutycycle-1)
    scale = snr*rms/fil.nchans/np.sqrt(fil.nspec*profmax*area)

    print "Duty cycle: %f" % dutycycle
    print "Expected SNR of injected pulsar signal (after folding " \
            "and integrating over frequency): %g" % snr
    prof.set_scale(scale)
    print "Scale factor applied: %g" % scale


def inject(infile, outfn, prof, period, dm, nbitsout=None, block_size=BLOCKSIZE):
    if isinstance(infile, filterbank.FilterbankFile):
        fil = infile
    else:
        fil = filterbank.FilterbankFile(infile, read_only=True)
    print "Injecting pulsar signal into: %s" % fil.filename
    if False:
        delays = psr_utils.delay_from_DM(dm, fil.frequencies)
        delays -= delays[np.argmax(fil.frequencies)]
        get_phases = lambda times: (times-delays)/period % 1
    else:
        get_phases = lambda times: times/period % 1

    # Create the output filterbank file
    if nbitsout is None:
        nbitsout = fil.nbits
    outfil = filterbank.create_filterbank_file(outfn, fil.header, nbits=nbitsout)

    if outfil.nbits == 8:
        raise NotImplementedError("This code is out of date. 'delays' is not " \
                                    "done in this way anymore..")
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
        del onesec
    else:
        # No scaling to be performed
        # These values will cause scaling to keep data unchanged
        minimum = 0
        global_scale = 1

    # Start an output file
    print "Creating out file: %s" % outfn
    sys.stdout.write(" %3.0f %%\r" % 0)
    sys.stdout.flush()
    oldprogress = -1
    
    # Loop over data
    lobin = 0
    spectra = fil.get_spectra(0, block_size)
    numread = spectra.shape[0]
    while numread:
        hibin = lobin+numread
        times = (np.arange(lobin, hibin)+0.5)*fil.dt
        phases = get_phases(times)
        toinject = prof(phases)
        injected = spectra+toinject[:,np.newaxis]
        scaled = (injected-minimum)*global_scale
        outfil.append_spectra(scaled)
        
        # Print progress to screen
        progress = int(100.0*hibin/fil.nspec)
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
        
        # Prepare for next iteration
        lobin = hibin 
        spectra = fil.get_spectra(lobin, lobin+block_size)
        numread = spectra.shape[0]

    sys.stdout.write("Done   \n")
    sys.stdout.flush()
    

def main():
    plt.figure(figsize=(8,10))
    fn = args[0]
    fil = filterbank.FilterbankFile(fn, read_only=True)
    if options.inprof is not None:
        print "Loading profile from file (%s)" % options.inprof
        infile = open(options.inprof, 'rb')
        prof = pickle.load(infile)
        infile.close()
    else:
        comps = create_vonmises_components(options.vonmises)
        print "Creating profile. Number of components: %d" % len(comps)
        prof = MultiComponentProfile(comps)
        if options.use_spline:
            prof = get_spline_profile(prof)
        if options.scale is None:
            scale = determine_and_set_scale(fil, prof, smean=options.smean, \
                        gain=options.gain, tsys=options.tsys, rms=options.rms)
        else:
            prof.set_scale(options.scale)
   
        if options.apply_dm:
            prof = apply_dm(prof, options.period, options.dm, \
                            np.abs(fil.foff), fil.frequencies, fil.tsamp)
        if options.outprof is not None:
            print "Writing %s instance to file (%s)" % \
                    (type(prof).__name__, options.outprof)
            outfile = open(options.outprof, 'wb')
            pickle.dump(prof, outfile, protocol=pickle.HIGHEST_PROTOCOL)
            outfile.close()

    outfn = options.outname % fil.header 
    if options.dryrun:
        print "Showing plot of profile to be injected..."
        plt.clf()
        prof.plot()
        plt.xlim(0,1)
        plt.xlabel("Phase")
        plt.savefig(outfn+".ps")
        plt.ylim(0,15)
        plt.savefig(outfn+".zoom.ps")
        sys.exit()

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
    parser.add_option("-F", "--smean", dest='smean', type='float', \
                    default=None, \
                    help="Mean flux density (in mJy) of the injected pulsar " \
                        "signal. This will automatically set the scale factor " \
                        "and thus the '-F/--smean' option _cannot_ be " \
                        "combined with the '-s/--scale' option. " \
                        "(Default: Don't scale)")
    parser.add_option("-s", "--scale", dest='scale', type='float', \
                    default=None, \
                    help="Overall scaling factor to multiply profile with. " \
                        "(Default: Don't scale.)")
    parser.add_option("--rms", dest='rms', type='float', default=None, \
                    help="The RMS of the recipient file's DM=0 timeseries. " \
                        "This _must_ be provided if the '-F/--smean' option " \
                        "is used.")
    parser.add_option("--gain", dest='gain', type='float', default=None, \
                    help="The telescope's gain (in K/Jy). This _must_ be " \
                        "provided if the '-F/--smean' option is used.")
    parser.add_option("--tsys", dest='tsys', type='float', default=None, \
                    help="The observing system's temperature (in K). This " \
                        "_must_ be provided if the '-F/--smean' option is used.")
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
    parser.add_option("--nbits", dest='output_nbits', default=None, type=int, \
                    help="Number of bits per same to use in output " \
                        "filterbank file. (Default: same as input file)")
    parser.add_option("-n", "--dryrun", dest="dryrun", action="store_true", \
                    help="Show the pulse profile to be injected and exit. " \
                        "(Default: do not show profile, inject it)")
    parser.add_option("--use-spline", dest='use_spline', action='store_true', \
                    default=False, \
                    help="Evaluate the analytic pulse profile and interpolate " \
                        "with a spline. This is typically faster to execute, " \
                        "especially when the profile is made up of multiple " \
                        "components. (Default: Do not use spline.)")
    parser.add_option("--no-apply-dm", dest='apply_dm', action='store_false', \
                    default=True, \
                    help="Do not apply the DM (i.e. do not delay or smear " \
                        "the pulse; Default: Apply DM)")
    parser.add_option("--load-prof", dest="inprof", default=None, \
                    help="Load a profile object from file. (Default: " \
                        "create a fresh profile object.)")
    parser.add_option("--save-prof", dest='outprof', default=None, \
                    help="Save a profile object to file. (Default: " \
                        "do not save profile object.)")
    parser.add_option("-o", "--outname", dest='outname', action='store', \
                    default="injected.fil", \
                    help="The name of the output file.")
    (options, args) = parser.parse_args()
    if options.period is None or options.dm is None:
        raise ValueError("Both a period and a DM _must_ be provided!")
    if options.smean is not None and options.scale is not None:
        raise ValueError("Only one of '-F/--smean' and '-s/--scale' " \
                        "options may be provided!")
    if options.smean is not None and (options.rms is None or \
                                        options.gain is None or \
                                        options.tsys is None):
        raise ValueError("When automatically determining the scaling " \
                        "factor for a particular Smean, the '--rms', " \
                        "'--gain', and '--tsys' options must all be " \
                        "provided!")
    
    if (options.smean is not None or options.scale is not None) and \
            options.inprof is not None:
        raise ValueError("Loading a saved profile via the " \
                        "'--load-prof' options is incompatible " \
                        "with adjusting the profile's scale.")
    main()
