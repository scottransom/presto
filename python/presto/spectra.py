from builtins import str
from builtins import range
from builtins import object
import copy

import numpy as np
import scipy.signal
from presto import psr_utils

class Spectra(object):
    """A class to store spectra. This is mainly to provide
        reusable functionality.
    """
    def __init__(self, freqs, dt, data, starttime=0, dm=0):
        """Spectra constructor.
            
            Inputs:
                freqs: Observing frequencies for each channel.
                dt: Sample time (in seconds).
                data: A 2D numpy array containing pulsar data.
                        Axis 0 should contain channels. (e.g. data[0,:])
                        Axis 1 should contain spectra. (e.g. data[:,0])
                starttime: Start time (in seconds) of the spectra
                        with respect to the start of the observation.
                        (Default: 0).
                dm: Dispersion measure (in pc/cm^3). (Default: 0)

            Output:
                spectra_obj: Spectrum object.
        """
        self.numchans, self.numspectra = data.shape
        assert len(freqs)==self.numchans

        self.freqs = freqs
        self.data = data.astype('float')
        self.dt = dt
        self.starttime = starttime
        self.dm = 0

    def __str__(self):
        return str(self.data)

    def __getitem__(self, key):
        return self.data[key]
    
    def __setitem__(self, key, value):
        self.data[key] = value
    
    def get_chan(self, channum):
        return self.data[channum,:]

    def get_spectrum(self, specnum):
        return self.data[:,specnum]

    def shift_channels(self, bins, padval=0):
        """Shift each channel to the left by the corresponding 
            value in bins, an array.

            Inputs:
                bins: An array containing the number of bins
                    to shift each channel by.
                padval: Value to use when shifting near the edge
                    of a channel. This can be a numeric value,
                    'median', 'mean', or 'rotate'. 
                    
                    The values 'median' and 'mean' refer to the 
                    median and mean of the channel. The value 
                    'rotate' takes values from one end of the 
                    channel and shifts them to the other.

            Outputs:
                None

            *** Shifting happens in-place ***
        """
        assert self.numchans == len(bins)
        for ii in range(self.numchans):
            chan = self.get_chan(ii)
            # Use 'chan[:]' so update happens in-place
            # this way the change effects self.data
            chan[:] = psr_utils.rotate(chan, bins[ii])
            if padval!='rotate':
                # Get padding value
                if padval=='mean':
                    pad = np.mean(chan)
                elif padval=='median':
                    pad = np.median(chan)
                else:
                    pad = padval
                
                # Replace rotated values with padval
                if bins[ii]>0:
                    chan[-bins[ii]:] = pad
                elif bins[ii]<0:
                    chan[:-bins[ii]] = pad

    def subband(self, nsub, subdm=None, padval=0):
        """Reduce the number of channels to 'nsub' by subbanding.
            The channels within a subband are combined using the
            DM 'subdm'. 'padval' is passed to the call to
            'Spectra.shift_channels'.

            Inputs:
                nsub: Number of subbands. Must be a factor of 
                    the number of channels.
                subdm: The DM with which to combine channels within
                    each subband (Default: don't shift channels 
                    within each subband)
                padval: The padding value to use when shifting
                    channels during dedispersion. See documentation
                    of Spectra.shift_channels. (Default: 0)

            Outputs:
                None

            *** Subbanding happens in-place ***
        """
        assert (self.numchans % nsub) == 0
        assert (subdm is None) or (subdm >= 0)
        nchan_per_sub = self.numchans // nsub
        sub_hifreqs = self.freqs[np.arange(nsub)*nchan_per_sub]
        sub_lofreqs = self.freqs[(1+np.arange(nsub))*nchan_per_sub-1]
        sub_ctrfreqs = 0.5*(sub_hifreqs+sub_lofreqs)
        
        if subdm is not None:
            # Compute delays
            ref_delays = psr_utils.delay_from_DM(subdm-self.dm, sub_ctrfreqs)
            delays = psr_utils.delay_from_DM(subdm-self.dm, self.freqs)
            rel_delays = delays-ref_delays.repeat(nchan_per_sub) # Relative delay
            rel_bindelays = np.round(rel_delays/self.dt).astype('int')
            # Shift channels
            self.shift_channels(rel_bindelays, padval)

        # Subband
        self.data = np.array([np.sum(sub, axis=0) for sub in \
                                np.vsplit(self.data, nsub)])
        self.freqs = sub_ctrfreqs
        self.numchans = nsub

    def scaled(self, indep=False):
        """Return a scaled version of the Spectra object.
            When scaling subtract the median from each channel,
            and divide by global std deviation (if indep==False), or
            divide by std deviation of each row (if indep==True).

            Input:
                indep: Boolean. If True, scale each row
                    independantly (Default: False).

            Output:
                scaled_spectra: A scaled version of the
                    Spectra object.
        """
        other = copy.deepcopy(self)
        if not indep:
            std = other.data.std()
        for ii in range(other.numchans):
            chan = other.get_chan(ii)
            median = np.median(chan)
            if indep:
                std = chan.std()
            chan[:] = (chan-median)/std
        return other
    
    def scaled2(self, indep=False):
        """Return a scaled version of the Spectra object.
            When scaling subtract the min from each channel,
            and divide by global max (if indep==False), or
            divide by max of each row (if indep==True).

            Input:
                indep: Boolean. If True, scale each row
                    independantly (Default: False).

            Output:
                scaled_spectra: A scaled version of the
                    Spectra object.
        """
        other = copy.deepcopy(self)
        if not indep:
            max = other.data.max()
        for ii in range(other.numchans):
            chan = other.get_chan(ii)
            min = chan.min()
            if indep:
                max = chan.max()
            chan[:] = (chan-min)/max
        return other

    def masked(self, mask, maskval='median-mid80'):
        """Replace masked data with 'maskval'. Returns
            a masked copy of the Spectra object.
            
            Inputs:
                mask: An array of boolean values of the same size and shape
                    as self.data. True represents an entry to be masked.
                maskval: Value to use when masking. This can be a numeric
                    value, 'median', 'mean', or 'median-mid80'.

                    The values 'median' and 'mean' refer to the median and
                    mean of the channel, respectively. The value 'median-mid80'
                    refers to the median of the channel after the top and bottom
                    10% of the sorted channel is removed.
                    (Default: 'median-mid80')

            Output:
                maskedspec: A masked version of the Spectra object.
        """
        assert self.data.shape == mask.shape
        maskvals = np.ones(self.numchans)
        for ii in range(self.numchans):
            chan = self.get_chan(ii)
            # Use 'chan[:]' so update happens in-place
            if maskval=='mean':
                maskvals[ii]=np.mean(chan)
            elif maskval=='median':
                maskvals[ii]=np.median(chan)
            elif maskval=='median-mid80':
                n = int(np.round(0.1*self.numspectra))
                maskvals[ii]=np.median(sorted(chan)[n:-n])
            else:
                maskvals[ii]=maskval
            if np.all(mask[ii]):
                self.data[ii] = np.ones_like(self.data[ii])*(maskvals[:,np.newaxis][ii])
        return self

    def dedisperse(self, dm=0, padval=0):
        """Shift channels according to the delays predicted by
            the given DM.
            
            Inputs:
                dm: The DM (in pc/cm^3) to use.
                padval: The padding value to use when shifting
                    channels during dedispersion. See documentation
                    of Spectra.shift_channels. (Default: 0)

            Outputs:
                None

            *** Dedispersion happens in place ***
        """
        assert dm >= 0
        ref_delay = psr_utils.delay_from_DM(dm-self.dm, np.max(self.freqs))
        delays = psr_utils.delay_from_DM(dm-self.dm, self.freqs)
        rel_delays = delays-ref_delay # Relative delay
        rel_bindelays = np.round(rel_delays/self.dt).astype('int')
        # Shift channels
        self.shift_channels(rel_bindelays, padval)

        self.dm=dm

    def smooth(self, width=1, padval=0):
        """Smooth each channel by convolving with a top hat
            of given width. The height of the top had is
            chosen shuch that RMS=1 after smoothing. 
            Overlap values are determined by 'padval'.

            Inputs:
                width: Number of bins to smooth by (Default: no smoothing)
                padval: Padding value to use. Possible values are
                    float-value, 'mean', 'median', 'wrap'.
                    (Default: 0).

            Ouputs:
                None

            This bit of code is taken from Scott Ransom's
            PRESTO's single_pulse_search.py (line ~ 423).
            
            *** Smoothing is done in place. ***
        """
        if width > 1:
            kernel = np.ones(width, dtype='float32')/np.sqrt(width)
            for ii in range(self.numchans):
                chan = self.get_chan(ii)
                if padval=='wrap':
                    tosmooth = np.concatenate([chan[-width:], \
                                chan, chan[:width]])
                elif padval=='mean':
                    tosmooth = np.ones(self.numspectra+width*2) * \
                                np.mean(chan)
                    tosmooth[width:-width] = chan
                elif padval=='median':
                    tosmooth = np.ones(self.numspectra+width*2) * \
                                np.median(chan)
                    tosmooth[width:-width] = chan
                else: # padval is a float
                    tosmooth = np.ones(self.numspectra+width*2) * \
                                padval
                    tosmooth[width:-width] = chan
                    
                smoothed = scipy.signal.convolve(tosmooth, kernel, 'same')
                chan[:] = smoothed[width:-width]
                    
    def trim(self, bins=0):
        """Trim the end of the data by 'bins' spectra.
            
            Input:
                bins: Number of spectra to trim off the end of the observation.
                    If bins is negative trim spectra off the beginning of the
                    observation.

            Outputs:
                None

            *** Trimming is irreversible ***
        """
        assert bins < self.numspectra
        if bins == 0:
            return
        elif bins > 0:
            self.data = self.data[:,:-bins]
            self.numspectra = self.numspectra-bins
        elif bins < 0:
            self.data = self.data[:,bins:]
            self.numspectra = self.numspectra-bins
            self.starttime = self.starttime+bins*self.dt

    def downsample(self, factor=1, trim=True):
        """Downsample (in-place) the spectra by co-adding
            'factor' adjacent bins.

            Inputs:
                factor: Reduce the number of spectra by this
                    factor. Must be a factor of the number of
                    spectra if 'trim' is False.
                trim: Trim off excess bins.

            Ouputs:
                None

            *** Downsampling is done in place ***
        """
        assert trim or not (self.numspectra % factor)
        new_num_spectra = self.numspectra // factor
        num_to_trim = self.numspectra%factor
        self.trim(num_to_trim)
        self.data = np.array(np.column_stack([np.sum(subint, axis=1) for \
                        subint in np.hsplit(self.data,new_num_spectra)]))
        self.numspectra = new_num_spectra
        self.dt = self.dt*factor
