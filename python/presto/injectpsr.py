#!/usr/bin/env python

"""Inject a fake pulsar into real data, creating
a filterbank file.

Patrick Lazarus, June 26, 2012
"""
from __future__ import print_function
from builtins import zip
from builtins import object
import sys
import argparse
import warnings
import copy

import numpy as np
import scipy.integrate
import scipy.interpolate
import matplotlib
matplotlib.use('agg') # Use a non-interactive backend
import matplotlib.pyplot as plt
import scipy.integrate

from presto import filterbank
from presto import psr_utils

DEBUG = False # Print debugging messages

NUMSECS = 1.0 # Number of seconds of data to use to determine global scale
              # when repacking floating-point data into integers
BLOCKSIZE = 1e4 # Number of spectra to manipulate at once
NUMPOINTS = 512 # Number of points to use for spline profiles when applying DM
NINTEG_PER_BIN = 1 # Number of points to average integrate per time bin

def integrate_phase_bin(prof_func, phs, dphs, nseg=1):
    intervals = np.linspace(0, dphs, nseg+1, endpoint=True)
    edges = intervals+np.asarray(phs)[...,np.newaxis]
    return np.trapz(prof_func(edges), edges, axis=-1)


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
        profvals = self.prof_func(phs)
        return self.scale*profvals

    def plot(self, nbin=1024, scale=1):
        x0 = np.linspace(0, 1.0, nbin+1, endpoint=True)
        plt.plot(x0, self(x0)*scale)
        plt.xlabel("Phase")

    def get_area(self, npts=4096, **kwargs):
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
        phs = np.linspace(0, 1.0, npts+1, endpoint=True)
        area = np.trapz(y=self(phs), x=phs)
        #area, err = scipy.integrate.quadrature(self, 0, 1, maxiter=250)
        return area

    def get_max(self, npts=4096):
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
       
    def get_fwhm(self, npts=4096):
        """Determine and return the FWHM of the profile, in phase.
            This only works if two points in the profile are at half-maximum,
            and all points in between are larger than half-max.

            Input:
                npts: The number of points to use when evaluating the
                    profile.

            Ouput:
                fwhm: The full-width at half-maximum of the profile, in phase.
        """
        phs = np.linspace(0, 1.0, npts+1, endpoint=True)
        vals = self(phs)
        profmax = np.max(self(phs))
        halfmax = profmax/2.0
        ma = np.ma.masked_less(vals, halfmax)
        # Unmasked region has vals >= halfmax
        unmasked = np.ma.notmasked_contiguous(ma)
        if len(unmasked) == 1:
            sl = unmasked[0]
            return ((sl.stop-1) - (sl.start))/float(npts)
        elif len(unmasked) == 2:
            sl1 = unmasked[0]
            sl2 = unmasked[1]
            if sl1.start==0 and sl2.stop==len(phs):
                # profile straddles phase=0
                return (npts+(sl1.stop-1) - sl2.start)/float(npts)
        else:
            return None


    def get_equivalent_width(self, npts=4096):
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

    def set_scaling(self, scale):
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
    
    def convolve_with(self, other, npts=4096, conserve_area=True):
        """Convolve Profile with another. Return a SplineProfile
            with the requested number of points.

            Inputs:
                other: The Profile to convolve with.
                npts: The number of points to use when creating the
                    resulting convolution (i.e. a SplineProfile). 
                    (Default: 4096)
                conserve_area: If true, scale the covolution such that
                    it has the same area as the input profile (i.e. 'self')

            Other:
                convolution: The convolution, a SplineProfile object.
        """
        phs = np.linspace(0, 1, npts, endpoint=False)
        conv_vals = np.fft.irfft(np.fft.rfft(self(phs))*np.fft.rfft(other(phs)))
        if conserve_area:
            # scale convolution such that area is same as before
            # Isn't this already a property of convolutions? Should it always be true?
            scale = self.get_area()/np.trapz(conv_vals, phs)
        convolution = SplineProfile(conv_vals*scale)
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
        self.profvals = profvals
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
    def __init__(self, profiles, scales=1):
        """Construct a vector of profiles.
            
            Inputs:
                profiles: A list of Profile objects.
                scales: A numpy array of scaling factors to multiply
                    the profile vector by when evaluating.

            Output:
                prof: The profile object.
        """
        self.profiles = profiles
        self.nprofs = len(profiles)
        self.scales = scales
    
    def __call__(self, phs):
        phs = np.atleast_1d(np.asarray(phs))
        nphs = phs.shape[-1]
        vals = np.zeros((nphs, self.nprofs))
        if phs.ndim == 1:
            # Evaluate all profiles at the same phases
            for ii, (prof, scale) in enumerate(zip(self.profiles, self.scales)):
                # If scale is not a finite number don't inject anything for this
                # channel
                if np.isfinite(scale):
                    vals[:,ii] = prof(phs)*scale
        elif phs.ndim == 2:
            # Evaluate each profile at a different set of phases
            nphs_vecs = phs.shape[0]
            if nphs_vecs != self.nprofs:
                raise ValueError("Length of axis=1 of 'phs' (%d) must be " \
                                "equal to the number of profiles in the " \
                                "vector (%d)." % (nphs_vecs, self.nprofs))
            for ii, (prof, ph, scale) in enumerate(zip(self.profiles, phs, self.scales)):
                # If scale is not a finite number don't inject anything for this
                # channel
                if np.isfinite(scale):
                    vals[:,ii] = prof(ph)*scale
        else:
            raise ValueError("VectorProfile can only be evaluated with " \
                            "1D or 2D arrays")
        return vals
    
    def set_scaling(self, scales):
        """Set the profiles' scaling factors.

            Input:
                scales: The scaling factors to use.

            Outputs:
                None
        """
        arr = np.asarray(scales)
        if arr.ndim not in (0, 1):
            raise ValueError("Bad number of dimensions for 'scales'. " \
                        "Number of dimensions must be 0 or 1. " \
                        "Scales provided has %d dimensions." % arr.ndim) 
        self.scales = arr

    def get_area(self, *args, **kwargs):
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
        areas = np.zeros(len(self.profiles))
        for ii, prof in enumerate(self.profiles):
            areas[ii] = prof.get_area(*args, **kwargs)
        return areas

    def plot(self, nbin=1024, scale=1):
        phs = np.linspace(0, 1.0, nbin+1, endpoint=True)
        data = self(phs).transpose()
        imax = plt.axes((0.1, 0.1, 0.6, 0.8))
         
        plt.imshow(data, interpolation='nearest', \
                    extent=(0, 1, 0, self.nprofs), aspect='auto')
        plt.xlabel("Phase")
        plt.ylabel("Channel number")
        
        plt.axes((0.7, 0.1, 0.25, 0.8), sharey=imax)
        plt.plot(np.sum(data, axis=1)[::-1], np.arange(self.nprofs), 'k-')
        
        # Re-set axes for image
        imax.set_xlim(0, 1)
        imax.set_ylim(0, self.nprofs)


class DispersedProfile(VectorProfile):
    def __init__(self, profiles, dm, freqs, period, scales=1, intrinsic=None, delayed=True):
        super(DispersedProfile, self).__init__(profiles, scales)
        self.dm = dm
        self.freqs = freqs
        self.period = period
        self.intrinsic = intrinsic
        self.delayed = delayed

    def get_equivalent_width(self, *args, **kwargs):
        return np.ones_like(self.freqs)*self.intrinsic.get_equivalent_width(*args, **kwargs)
    
    def get_area(self, *args, **kwargs):
        return np.ones_like(self.freqs)*self.intrinsic.get_area(*args, **kwargs)

    def get_max(self, *args, **kwargs):
        return np.ones_like(self.freqs)*self.intrinsic.get_max(*args, **kwargs)

    def plot(self, nbin=1024, scale=1, dedisp=False):
        phs = np.linspace(0, 1.0, nbin+1, endpoint=True)
        data = self(phs).transpose()
        if self.delayed:
            delays = get_phasedelays(self.dm, self.freqs, self.period)
        else:
            delays = np.zeros_like(self.freqs)
        delayedphs = (phs+delays[:,np.newaxis]) % 1
        dedispdata = self(delayedphs).transpose() 

        imax = plt.axes((0.1, 0.1, 0.6, 0.6))
 
        if dedisp:
            plt.imshow(dedispdata, interpolation='nearest', \
                       extent=(0, 1, 0, self.nprofs), aspect='auto')
        else:
            plt.imshow(data, interpolation='nearest', \
                       extent=(0, 1, 0, self.nprofs), aspect='auto')
        plt.set_cmap('gist_yarg')
        plt.xlabel("Phase")
        plt.ylabel("Channel number")
 
        plt.axes((0.7, 0.1, 0.25, 0.6), sharey=imax)
        plt.plot(np.sum(data, axis=1)[::-1], np.arange(self.nprofs), 'k-')
 
        plt.axes((0.1, 0.7, 0.6, 0.25), sharex=imax)
        plt.plot(phs, np.sum(dedispdata, axis=0), ls='-', 
                 c='k', lw=1, label='Smeared and scattered', zorder=2)
        if self.intrinsic is not None:
            plt.plot(phs, self.intrinsic(phs)*np.ma.masked_invalid(self.scales).sum(), 
                     ls='-', c='#bbbbbb', lw=3, label='Input', zorder=1)
            plt.legend(loc='best', prop=dict(size='small'))
        plt.figtext(0.05, 0.05, "Period = %.3f ms" % (self.period*1000), size='xx-small')
        plt.figtext(0.05, 0.035, r"DM = %.3f cm$\mathrm{^{-3}}$pc" % self.dm, size='xx-small')

        # Re-set axes for image
        imax.set_xlim(0, 1)
        imax.set_ylim(0, self.nprofs)


def get_phasedelays(dm, freqs, period):
    """Return phase delays corresponding to a particular DM.

        Inputs:
            dm: DM (in pc cm-3)
            freqs: The list of frequencies (in MHz)
            period: The profiles period (in seconds)

        Outputs:
            phasedelays: The corresponding phase delays.
    """
    # Prepare delays
    timedelays = psr_utils.delay_from_DM(dm, freqs)
    # Reference all delays to highest frequency channel, which remains
    # unchanged
    # TODO: Do we really want to refer to high freq?
    timedelays -= timedelays[np.argmax(freqs)]
    phasedelays = timedelays/period
    return phasedelays


def apply_dm(inprof, period, dm, chan_width, freqs, tsamp, \
                do_delay=True, do_smear=True, do_scatter=True,
                verbose=True):
    """Given a profile apply DM delays, smearing, and scattering 
        within each channel as is appropriate for the given params.

        Inputs:
            inprof: The profile to modify.
            period: The profiles period (in seconds)
            dm: The DM (in pc cm-3)
            chan_width: The width of each channel (in MHz)
            freqs: The list of frequencies (in MHz)
            tsamp: Sample time of the recipient filterbank file (in seconds).
            do_delay: Boolean, if True apply DM delays to each channel.
                The highest freq channel is not shifted. (Default: True)
            do_smear: Boolean, if True apply DM smearing to each channel.
                (Default: True)
            do_scatter: Boolean, if True apply scattering to each channel.
                (Default: True)

        Outputs:
            vecprof: The delayed and smeared VectorProfile.
    """
    weq = inprof.get_equivalent_width()
    nfreqs = len(freqs)
    if verbose:
        print("Applying DM to profile (DM = %.2f; %d channels)..." % \
                (dm, nfreqs))
    # A list of profiles, one for each channel
    profiles = []

    if dm <= 0:
        warnings.warn("DM will not be applied because it is 0 (or smaller?!)")
        do_delay = False
        do_smear = False
        do_scatter = False

    if do_delay:
        phasedelays = get_phasedelays(dm, freqs, period)
    else:
        phasedelays = np.zeros(nfreqs)

    # Prepare for smear campaign
    smeartimes = psr_utils.dm_smear(dm, abs(chan_width), freqs) # In seconds
    smearphases = smeartimes/period
    
    # Prepare to scatter
    scattertimes = psr_utils.pulse_broadening(dm, freqs)*1e-3 # In seconds
    scatterphases = scattertimes/period

    if DEBUG:
        for ichan, (freq, smear, scatt, delay) in \
                enumerate(zip(freqs, smearphases, scatterphases, phasedelays)):
            print("    Chan #%d - Freq: %.3f MHz -- " \
                  "Smearing, scattering, delay (all in phase): " \
                  "%g, %g, %g" % (ichan, freq, smear, scatt, delay))
    oldprogress = 0
    sys.stdout.write(" %3.0f %%\r" % oldprogress)
    sys.stdout.flush()
#    ylim = None
#    ylim2 = None
#    ylim3 = None
#    ylim4 = None
#    ylim5 = None
    for ii, (delayphs, smearphs, scattphs) in \
                enumerate(zip(phasedelays, smearphases, scatterphases)):
        #########
        # DEBUG: plot all profiles
#        plt.clf()
#        ax = plt.subplot(5,1,1)
#        inprof.plot()
#        if ylim is not None:
#            ax.set_ylim(ylim)
#        else:
#            ylim = ax.get_ylim()

        if do_smear and not ((smearphs < 0.2*weq) or (smearphs < (tsamp/period))):
            # Only smear if requested and smearing-phase is large enough
#            bc = boxcar_factory(smearphs, delayphs)
#            ax2 = plt.subplot(5,1,2,sharex=ax)
#            bc.plot()
#            if ylim2 is not None:
#                ax2.set_ylim(ylim2)
#            else:
#                ylim2 = ax2.get_ylim()
            if DEBUG:
                print("Smearing")
            tmpprof = inprof.smear(smearphs, delayphs, npts=NUMPOINTS)
        else:
            tmpprof = inprof.delay(delayphs)
            phs = np.linspace(0, 1, NUMPOINTS+1)
            tmpprof = SplineProfile(tmpprof(phs))
#        ax3 = plt.subplot(5,1,3,sharex=ax)
#        if ylim3 is not None:
#            ax3.set_ylim(ylim3)
#        else:
#            ylim3 = ax3.get_ylim()
#        tmpprof.plot()
        if do_scatter and not ((scattphs < 0.2*weq) or (scattphs < (tsamp/period))):
            # Only scatter if requested and scattering-phase is large enough
#            ex = exponential_factory(scattphs)
#            ax4 = plt.subplot(5,1,4,sharex=ax)
#            ex.plot()
#            if ylim4 is not None:
#                ax4.set_ylim(ylim4)
#            else:
#                ylim4 = ax4.get_ylim()
            if DEBUG:
                print("Scattering")
            tmpprof = tmpprof.scatter(scattphs, npts=NUMPOINTS)
#        ax5 = plt.subplot(5,1,5,sharex=ax)
#        tmpprof.plot()
#        if ylim5 is not None:
#            ax5.set_ylim(ylim5)
#        else:
#            ylim5 = ax5.get_ylim()
        profiles.append(tmpprof)
#        plt.xlim(0,1)
#        plt.xlabel("Phase")
#        plt.suptitle("Prof %d (%f MHz)" % (ii, freqs[ii]))
#        plt.savefig("prof%d.png" % ii)
        #########
        # Print progress to screen
        progress = int(100.0*ii/nfreqs)
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
    sys.stdout.write("Done   \n")
    sys.stdout.flush()
    dispersedprof = DispersedProfile(profiles, dm=dm, freqs=freqs, period=period, 
                                     intrinsic=inprof, delayed=do_delay)
    return dispersedprof


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


def scale_from_snr(fil, prof, snr, rms):
    """Set the profile's scaling factor such that the simulated 
        injected pulsar signal will have the given Smean.

        Inputs:
            fil: A FilterbankFile object.
            prof: The Profile object representing the profile to inject.
            snr: The desired signal-to-noise ratio
            rms: The RMS of the recipient file's DM=0 time series.

        Outputs:
            scale: The scaling factor to apply to the profile.
    """
    # Set scale to 1, just in case it's be altered already.
    prof.set_scaling(1)
    
    # Characterise the pulse
    area = prof.get_area()
    profmax = prof.get_max()

    scale = snr*rms/fil.nchans/np.sqrt(fil.nspec*profmax*area)
    print("Average area %s, average profile maximum: %s" % \
            (np.mean(area), np.mean(profmax)))
    print("Average recommended scale factor: %s" % np.mean(scale))
    return scale


def snr_from_smean(fil, prof, smean, gain, tsys):
    """Set the profile's scaling factor such that the simulated 
        injected pulsar signal will have the given Smean.

        Inputs:
            fil: A FilterbankFile object.
            prof: The Profile object representing the profile to inject.
            smean: The mean flux density to simulate, in mJy.
            gain: The telescope's gain, in K/Jy.
            tsys: The observing system's temperature, in K.

        Outputs:
            snr: The target signal-to-noise ratio.
    """
    dutycycle = prof.get_equivalent_width() # in phase
                        # Because width is reported in phase, 
                        # it is actually the duty cycle

    # Characterise the recipient filterbank file
    tint = fil.nspec*fil.tsamp
    bw = np.abs(fil.foff*fil.nchans)

    # Target SNR
    warnings.warn("Assuming 2 (summed) polarizations.")
    snr = smean*gain*np.sqrt(2*tint*bw)/tsys*np.sqrt(1/dutycycle-1)
    print("Expected SNR of injected pulsar signal (after folding " \
            "and integrating over frequency): %s" % snr)
    return snr


def inject(infile, outfn, prof, period, dm, nbitsout=None, 
           block_size=BLOCKSIZE, pulsar_only=False, inplace=False):
    if isinstance(infile, filterbank.FilterbankFile):
        fil = infile
    elif inplace:
        fil = filterbank.FilterbankFile(infile, 'readwrite')
    else:
        fil = filterbank.FilterbankFile(infile, 'read')
    print("Injecting pulsar signal into: %s" % fil.filename)
    if False:
        delays = psr_utils.delay_from_DM(dm, fil.frequencies)
        delays -= delays[np.argmax(fil.frequencies)]
        get_phases = lambda times: (times-delays)/period % 1
    else:
        get_phases = lambda times: times/period % 1

    # Create the output filterbank file
    if nbitsout is None:
        nbitsout = fil.nbits
    if inplace:
        warnings.warn("Injecting pulsar signal *in-place*")
        outfil = fil
    else:
        # Start an output file
        print("Creating out file: %s" % outfn)
        outfil = filterbank.create_filterbank_file(outfn, fil.header, \
                                            nbits=nbitsout, mode='append')

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

    sys.stdout.write(" %3.0f %%\r" % 0)
    sys.stdout.flush()
    oldprogress = -1
    
    # Loop over data
    lobin = 0
    spectra = fil.get_spectra(0, block_size)
    numread = spectra.shape[0]
    while numread:
        if pulsar_only:
            # Do not write out data from input file
            # zero it out
            spectra *= 0
        hibin = lobin+numread
        # Sample at middle of time bin
        times = (np.arange(lobin, hibin, 1.0/NINTEG_PER_BIN)+0.5/NINTEG_PER_BIN)*fil.dt
        #times = (np.arange(lobin, hibin)+0.5)*fil.dt
        phases = get_phases(times)
        profvals = prof(phases)
        shape = list(profvals.shape)
        shape[1:1] = [NINTEG_PER_BIN]
        shape[0] /= NINTEG_PER_BIN
        profvals.shape = shape
        toinject = profvals.mean(axis=1)
        #toinject = profvals
        if np.ndim(toinject) > 1:
            injected = spectra+toinject
        else:
            injected = spectra+toinject[:,np.newaxis]
        scaled = (injected-minimum)*global_scale
        if inplace:
            outfil.write_spectra(scaled, lobin)
        else:
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


def load_profile(infn, verbose=True):
    if verbose:
        print("Loading profile from file (%s)" % infn)
    data = np.load(infn)
    profiles = []
    for key in sorted(data.keys()):
        profiles.append(SplineProfile(data[key]))
    prof = VectorProfile(profiles)
    data.close()
    return prof


def save_profile(prof, outfn, verbose=True):
    if verbose:
        print("Writing %s instance to file (%s)" % \
                (type(prof).__name__, outfn))
    outfile = open(outfn, 'wb')
    profvals = {}
    for ii, pp in enumerate(prof.profiles):
        profvals['chan%d' % ii] = pp.profvals
    np.savez(outfile, **profvals)
    outfile.close()


def parse_cfgstr(cfgstrs):
    cfgs = {}
    for cfgstr in cfgstrs:
        for cfg in cfgstr.split(','):
            key, val = cfg.split('=')
            cfgs[key] = val
    return cfgs


def get_scaling(fil, prof, cfgstrs):
    """Given a target filterbank file, a profile, and
        a configuration string return the corresponding 
        scaling factor.

        Inputs:
            fil: A filterbank.FilterbankFile object.
            prof: A Profile object.
            cfgstrs: A list of strings containing configurations.

        Output:
            scaling: The corresponding scaling.
    """
    cfgs = parse_cfgstr(cfgstrs)
    return float(cfgs['scale'])


def get_scaling_from_snr(fil, prof, cfgstrs):
    """Given a target filterbank file, a profile, and
        a configuration string compute the scaling factor
        given the SNR.

        Inputs:
            fil: A filterbank.FilterbankFile object.
            prof: A Profile object.
            cfgstrs: A list of strings containing configurations.

        Output:
            scaling: The corresponding scaling.
    """
    cfgs = parse_cfgstr(cfgstrs)
    snr = float(cfgs['snr'])
    rms = float(cfgs['rms'])
    scale = scale_from_snr(fil, prof, snr=snr, rms=rms)
    return scale


def get_scaling_from_smean(fil, prof, cfgstrs):
    """Given a target filterbank file, a profile, and
        a configuration string compute the scaling factor
        given the target mean flux density.

        Inputs:
            fil: A filterbank.FilterbankFile object.
            prof: A Profile object.
            cfgstrs: A list of strings containing configurations.

        Output:
            scaling: The corresponding scaling.
    """
    cfgs = parse_cfgstr(cfgstrs)
    smean = float(cfgs['smean'])
    rms = float(cfgs['rms'])
    gain = float(cfgs['gain'])
    tsys = float(cfgs['tsys'])
    snr = snr_from_smean(fil, prof, smean=smean, \
                            gain=gain, tsys=tsys)
    scale = scale_from_snr(fil, prof, snr, rms=rms)
    return scale
    

def get_scaling_from_file(fil, prof, cfgstrs):
    """Given a target filterbank file, a profile, and
        a configuration string read scaling factors from
        a text file. The file should have one floating point
        number per line. There should be as many lines as there
        are frequency channels in the filterbank file.

        Inputs:
            fil: A filterbank.FilterbankFile object.
            prof: A Profile object.
            cfgstrs: A list of strings containing configurations.

        Output:
            scaling: The corresponding scaling.
    """
    cfgs = parse_cfgstr(cfgstrs)
    fn = cfgs['file']
    smean = float(cfgs['smean'])
    col = int(cfgs.get('col', 0))

    area = prof.get_area()

    # 'conversion_factors' contains the per-channel conversion between mJy 
    # and receiver units (in units of mJy/rcvr)
    conversion_factors = np.loadtxt(fn, usecols=(col,))

    scales = smean/conversion_factors/area
    return scales


SCALE_METHODS = {'scale': get_scaling, \
                 'snr': get_scaling_from_snr, \
                 'radiometer': get_scaling_from_smean, \
                 'scalefile': get_scaling_from_file}

def scale_profile(prof, scale_name, scale_cfgstrs, fil, verbose=True):
    """Scale the profile.
        
        Inputs:
            prof: The profile object to scale.
            scale_name: The name of the type of scaler to use.
            scale_cfgstrs: A list of configuration strings for the scaler.
            fil: A FilterbankFile object.
            verbose: Print extra information. (Default: True)

        Outputs:
            None - the scalings are applied directly to the input profile.
    """
    # Determine scaling
    scale_getter = SCALE_METHODS[scale_name]
    scaling = scale_getter(fil, prof, scale_cfgstrs)
    if verbose:
        print("Band-averaged scale-factor: %g" % np.ma.masked_invalid(scaling).mean())
    prof.set_scaling(scaling)


def make_profile(vonmises, verbose=True):
    """Create the profile object to use.
        
        Inputs:
            vonmises: A list of strings defining von mises components.
            verbose: Print extra information. (Default: True)
    """
    comps = create_vonmises_components(vonmises)
    prof = MultiComponentProfile(comps)
    if verbose:
        print("Creating profile. Number of components: %d" % len(comps))
        print("Profile area (intensity x phase): %g" % prof.get_area())
        print("Equivalent width (phase): %g" % prof.get_equivalent_width())
        print("FWHM (phase): %g" % prof.get_fwhm())
        print("Profile maximum: %g" % prof.get_max())
    prof = get_spline_profile(prof)
    return prof


def main():
    fn = args.infile
    if args.inplace:
        fil = filterbank.FilterbankFile(fn, mode='readwrite')
    else:
        fil = filterbank.FilterbankFile(fn, mode='read')
    if args.inprof is not None:
        warnings.warn("Saved profiles already may be tuned to a particular " \
                        "DM, period and filterbank file (freq, nchans, " \
                        "tsamp, etc).")
        prof = load_profile(args.inprof)
    else:
        prof = make_profile(args.vonmises)
        prof = apply_dm(prof, args.period, args.dm, \
                        fil.foff, fil.frequencies, fil.tsamp)
        scale_profile(prof, args.scale_name, args.scale_cfgstrs, fil)
        if args.outprof is not None:
            save_profile(prof, args.outprof)

    outfn = args.outname % fil.header 
    print("Showing plot of profile to be injected...")
    plt.figure()
    plt.clf()
    prof.plot(dedisp=True)
    plt.xlim(0,1)
    plt.savefig(outfn+".ps")
    if args.dryrun:
        sys.exit()

    inject(fil, outfn, prof, args.period, args.dm, \
            nbitsout=args.output_nbits, block_size=args.block_size, \
            pulsar_only=args.pulsar_only, inplace=args.inplace)


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


class ParseMfileAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        vonmises = getattr(namespace, self.dest)
        vonmises.extend(parse_model_file(values))


class ScaleHelpAction(argparse.Action):
    def __call__(self, parser, namepsace, values, option_string=None):
        helptext = "Scaling methods:\n\n" + \
        "scale - Multiply injected signal by a scaling factor\n" + \
        "Configs: 'scale' - (float) the scaling factor\n\n" + \
        "snr - Scale the injected signal so the integrated \n"  + \
        "      profile has the given SNR\n" + \
        "Configs: 'snr' - (float) the target SNR\n" + \
        "         'rms' - (float) the RMS of the cleaned DM=0 timeseries\n\n" + \
        "radiometer - Scale the injected signal so the integrated \n" + \
        "             profile has the given mean flux density\n" + \
        "Configs: 'smean' - (float) the target mean flux density (in mJy)\n" + \
        "         'gain'  - (float) the telescope's gain (in K/Jy)\n" + \
        "         'tsys'  - (float) the observing system's temperature (in K)\n" + \
        "         'rms'   - (float) the RMS of the cleaned DM=0 timeseries\n\n" + \
        "scalefile - Scale the signal in each channel independently\n" + \
        "            according to the scaling factors in the file\n" + \
        "Configs: 'smean' - (float) the target mean flux density\n" + \
        "         'file'  - (string) a text file containing per-channel \n" + \
        "                   conversion factors from flux-density to 'receiver' \n" + \
        "                   units (in mJy/rcvr).\n" + \
        "                   Each row should have a single floating-point \n"  + \
        "                   number. The number and order of the rows should \n" + \
        "                   correspond to the input filterbank file.\n" + \
        "         'col'   - (int) The column number to read conversion factors \n" +\
        "                   from. Columns are numbered starting at 0.\n" + \
        "                   (Default: first column)\n\n"
        sys.stderr.write(helptext)
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='injectpsr.py', \
                    description="v0.7 Patrick Lazarus (Jan 16, 2014)")
    parser.add_argument("--dm", dest='dm', type=float, \
                    help="The DM of the (fake) injected pulsar signal. " \
                        "(This argument is required.", \
                    default=None)
    parser.add_argument("-p", "--period", dest='period', \
                    default=None, type=float, \
                    help="The *topocentric* period (in seconds) of the (fake) injected " \
                        "pulsar signal. (This argument is required.)")
    parser.add_argument("-c", "--scale-configs", dest='scale_cfgstrs', type=str, \
                    required=True, default=[], action='append', \
                    help="A string of comma-separated parameters to " \
                        "configure how the injected signal is scaled. " \
                        "Format is '<param-name>=<value>,...'")
    parser.add_argument("-s", "--scale-method", dest='scale_name', type=str, \
                    required=True, \
                    help="Name of a method for determining how the " \
                        "injected signal is scaled.")
    parser.add_argument("--scale-help", dest='show_scale_help', \
                    nargs=0, action=ScaleHelpAction, \
                    help="Show help text for scaling methods, parameters.")
    parser.add_argument("-v", "--vonmises", dest='vonmises', action='append', \
                    help="A string of 3 parameters defining a vonmises " \
                        "component to be injected. Be sure to quote the " \
                        "3 parameters together. The params are: 'amplitude " \
                        "shape phase'. Amplitude is not related to SNR in " \
                        "any way. Also, phase should be between 0 and 1. " \
                        "(Default: if no compoments are provided " \
                        "a von Mises with amplitude=1.0, shape=5, and " \
                        "phase=0.5 will be used.)", \
                    default=[])
    parser.add_argument("-m", "--model-file", dest="vonmises", \
                    type=str, action=ParseMfileAction, \
                    help="A model file (*.m) as written by 'paas'.")
    parser.add_argument("--block-size", dest='block_size', default=BLOCKSIZE, \
                    type=float, \
                    help="Number of spectra per block. This is the amount " \
                        "of data manipulated/written at a time. (Default: " \
                        " %d spectra)" % BLOCKSIZE)
    parser.add_argument("--nbits", dest='output_nbits', default=None, type=int, \
                    help="Number of bits per same to use in output " \
                        "filterbank file. (Default: same as input file)")
    parser.add_argument("-n", "--dryrun", dest="dryrun", action="store_true", \
                    help="Show the pulse profile to be injected and exit. " \
                        "(Default: do not show profile, inject it)")
    parser.add_argument("--no-apply-dm", dest='apply_dm', action='store_false', \
                    default=True, \
                    help="Do not apply the DM (i.e. do not delay or smear " \
                        "the pulse; Default: Apply DM)")
    parser.add_argument("--load-prof", dest="inprof", default=None, \
                    help="Load a profile object from file. (Default: " \
                        "create a fresh profile object.)")
    parser.add_argument("--save-prof", dest='outprof', default=None, \
                    help="Save a profile object to file. (Default: " \
                        "do not save profile object.)")
    parser.add_argument("-o", "--outname", dest='outname', \
                    default="injected.fil", \
                    help="The name of the output file.")
    parser.add_argument("--write-pulsar-only", dest='pulsar_only', \
                    action='store_true', \
                    help="Only write the pulsar signal to the output file. "
                         "That is, do not include the data from the input "
                         "file. This is useful for debugging. (Default: "
                         "write data from input file _and_ pulsar signal.)")
    parser.add_argument("--in-place", dest='inplace', \
                    action='store_true', \
                    help="Inject the pulsar signal in-place. " \
                         "THIS WILL OVERWRITE THE INPUT DATAFILE!" \
                         "(Default: Do _not_ inject in-place)")
    parser.add_argument("infile", \
                    help="File that will receive synthetic pulses.")
    args = parser.parse_args()
    if args.period is None or args.dm is None:
        raise ValueError("Both a period and a DM _must_ be provided!")
    if args.scale_name is not None and args.inprof is not None:
        raise ValueError("Loading a saved profile via the " \
                        "'--load-prof' args is incompatible " \
                        "with scaling the profile.")
    main()
