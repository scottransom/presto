from __future__ import absolute_import
from builtins import object
import numpy as np
from . import bary_and_topo


def topo_timeshift(bary_start_time, time_shift, topo):
    ind = np.where(topo == float(int(bary_start_time)/10*10))[0]
    return time_shift[ind]


def numsub(nchans, snr):
    if not nchans==960 and np.log2(nchans).is_integer(): #Puppi L-wide and GBNCC 
        if snr < 10:
            nsub = 32
        elif snr >= 10 and snr < 15:
            nsub = 64
        else:
            nsub = 128
    elif nchans == 960: #PALFA
        if snr < 10:
            nsub = 32
        elif snr >= 10 and snr < 15:
            nsub = 64
        else:
            nsub = 96
    else:
        nsub = nchans
    return nsub    


class params(object):
    """
    A class for input parameters for the spd plots.

    """
    def __init__(self):
        self.start_bin = 0
        self.dmfac = 0.0
        self.duration = 0.0
        self.nbins = 0
        self.nbinsextra = 0
        self.zerodm = None
        self.mask = False
        self.bandpass_corr = True
        self.nsub = 96
        self.subdm = 0.0
        self.dm = 0.0
        self.sweep_dm = 0.0
        self.downsamp = 1
        self.scaleindep = None
        self.width_bins = 1
        self.binratio = 50
        self.sigma = 0.0
        self.pulse_width = 0.0
        self.start = 0.0
        self.sweep_duration = 0.0
        self.bary_start_time = 0.0
        self.topo_start_time = 0.0
        self.sample_number = 0

    def read_from_file(self, params, tsamp, N, lofreq, hifreq, rawdatafile, loc_pulse = 0.5, dedisp = False, \
                       scaleindep = None, zerodm = None, mask = None, barytime = True, \
                       nsub = None, bandpass_corr = False): 
        """
           Set up parameters based on input from the groups.txt file.
           Input: params: list of parameters (DM, S/N, Time, Sample number, downfactor)
                  tsamp: sampling time (downsampled: 65.5 us for PALFA)
                  N    : The total number of time samples in the observation
                  lofreq: lowest observation frequency
                  hifreq: highest observation frequency
                  rawdatafile: supply a PSRFITS file instance
           optional arguments:
                  dedisp: Do you want to dedisperse?(Type: Boolean).
                  scaleindep:Do you want to scale each subband independently?(Type: Boolean)
                  zerodm:Do you want to use zero-DM filtering?(Type: Boolean)
                  mask: Do you want to use a rfifind mask? (Type: Boolean)
                  barytime: Is the given time(s) barycentric?
                  bandpass_corr:Would you like to remove the bandpass? (Type: Boolean)
        """
        self.subdm = params[0]
        self.sigma = params[1]
        topo, bary = bary_and_topo.bary_to_topo(rawdatafile.filename, rawdatafile=rawdatafile)
        time_shift = bary-topo
        if barytime:
            #### if the times in the groups.txt file are barycentric then get the corresponding
            #### topocentric times.
            self.bary_start_time = params[2]
            self.topo_start_time = self.bary_start_time - topo_timeshift(self.bary_start_time, \
                                                                         time_shift, topo)[0]
        else:
            #### Vice versa.
            self.topo_start_time = params[2]
            self.bary_start_time = self.topo_start_time + topo_timeshift(self.bary_start_time, \
                                                                         time_shift, topo)[0]
        self.sample_number = params[3]
        self.width_bins = params[4]
        sweep_posn = 0.0
        self.scaleindep = scaleindep
        self.downsamp = np.round((params[2]/self.sample_number/tsamp)).astype('int')
        self.duration = self.binratio * self.width_bins * tsamp * self.downsamp
        self.start = self.topo_start_time - (loc_pulse * self.duration)
        if (self.start<0.0):
            self.start = 0.0
        self.start_bin = np.round(self.start/tsamp).astype('int')
        self.pulse_width = self.width_bins*self.downsamp*tsamp
        if nsub is not None:
            self.nsub = nsub
        else:
            self.nsub = numsub(rawdatafile.nchan, self.sigma)
        self.zerodm = zerodm
        if dedisp:
            self.dm = self.subdm
            self.sweep_dm = None
            self.sweep_duration = 0.0
            self.dmfac = 4.15e3 * np.abs(1./lofreq**2 - 1./hifreq**2)
            self.nbins = np.round(self.duration/tsamp).astype('int')
            self.nbinsextra = np.round((self.duration + self.dmfac * self.dm)/tsamp).astype('int')
            if (self.start_bin+self.nbinsextra) > N-1:
                    self.nbinsextra = N-1-self.start_bin
            self.bandpass_corr = bandpass_corr
        else:
            self.dm = None
            self.sweep_dm = self.subdm
            self.sweep_duration = 4.15e3 * np.abs(1./lofreq**2-1./hifreq**2)*self.sweep_dm
            self.start = self.start + (loc_pulse*self.duration)
            self.start_bin = np.round(self.start/tsamp).astype('int')
            self.nbins = np.round(self.sweep_duration/tsamp).astype('int')
            self.nbinsextra = self.nbins
            if (self.start_bin+self.nbinsextra) > N-1:
                self.nbinsextra = N-1-self.start_bin
            self.bandpass_corr = False

    def manual_params(self, subdm, dm, sweep_dm, sigma, start_time, width_bins, downsamp, \
                      duration, nbins, nsub, tsamp, N, lofreq, hifreq, rawdatafile, loc_pulse=0.5, dedisp = False, \
                      scaleindep = None, zerodm = None, mask = False, barytime = True, \
                      bandpass_corr = False): 
        """
           Set up parameters based on input from the groups.txt file.
           Input:
                  subdm: DM to use when subbanding.
                  dm: DM to use when dedispersing data for plot. 
                  sweep_dm: Show the frequency sweep using this DM.
                  sigma: signal-to-noise of the pulse
                  start_time: start time of the data to be read in for waterfalling.
                  width_bins: Smooth each channel/subband with a boxcar width_bins wide.
                  downsamp: Factor to downsample in time by. Default: Don't downsample.
                  duration: duration of data to be waterfalled.
                  nbins: Number of time bins to plot. This option overrides
                          the duration argument. 
                  nsub: Number of subbands to use. Must be a factor of number of channels.
                  tsamp: sampling time (downsampled: 65.5 us for PALFA)
                  N: total number of samples in an observations
                  lofreq: lowest observation frequency
                  hifreq: highest observation frequency
                  rawdatafile: supply a psrfits file instance
           optional arguments:
                  dedisp: Do you want to dedisperse?(Type: Boolean).
                  scaleindep:Do you want to scale each subband independently?(Type: Boolean)
                  zerodm:Do you want to use zero-DM filtering?(Type: Boolean)
                  mask: Do you want to use a rfifind mask? (Type: Boolean)
                  barytime: Is the given time(s) barycentric?
                  bandpass_corr:Would you like to remove the bandpass? (Type: Boolean)
        """
        self.subdm = subdm
        self.mask = mask
        self.sigma = sigma
        topo, bary = bary_and_topo.bary_to_topo(rawdatafile.filename, rawdatafile=rawdatafile)
        time_shift = bary-topo
        if barytime:
            #### if the time is barycentric then get the corresponding topocentric time.
            self.bary_start_time = start_time
            self.topo_start_time = self.bary_start_time - topo_timeshift(self.bary_start_time, \
                                                                         time_shift, topo)[0]
        else:
            #### Vice versa.
            self.topo_start_time = start_time
            self.bary_start_time = self.topo_start_time + topo_timeshift(self.bary_start_time, \
                                                                         time_shift, topo)[0]

        self.sample_number = np.round(self.bary_start_time/tsamp).astype('int')
        self.width_bins = width_bins
        sweep_posn = 0.0
        if scaleindep:
            self.scaleindep = True
        self.downsamp = downsamp
        if nbins:
            self.nbins = nbins
            self.duration = self.nbins*tsamp
        if duration:
            self.duration = duration
            self.nbins = np.round(self.duration/tsamp).astype('int')
        self.start = self.topo_start_time - (loc_pulse * self.duration)
        if (self.start<0.0):
            self.start = 0.0
        self.start_bin = np.round(self.start/tsamp).astype('int')
        self.pulse_width = self.width_bins*self.downsamp*tsamp
        self.nsub = nsub
        self.zerodm = zerodm
        if dedisp:
            if dm:
                self.dm = dm
            else: 
                self.dm = self.subdm
            self.sweep_dm = None
            self.sweep_duration = 0.0
            self.dmfac = 4.15e3 * np.abs(1./lofreq**2 - 1./hifreq**2)
            self.nbinsextra = self.nbins + np.round((self.dmfac * self.dm)/tsamp).astype('int')
            if (self.start_bin+self.nbinsextra) > N-1:
                    self.nbinsextra = N-1-self.start_bin
            self.bandpass_corr = bandpass_corr
        else:
            self.dm = None
            if not sweep_dm:
                self.sweep_dm = self.subdm
            else:
                self.sweep_dm = sweep_dm
            self.sweep_duration = 4.15e3 * np.abs(1./lofreq**2-1./hifreq**2)*self.sweep_dm
            self.start = self.start + (loc_pulse*self.duration)
            self.start_bin = np.round(self.start/tsamp).astype('int')
            self.nbinsextra = self.nbins


            if (self.start_bin+self.nbinsextra) > N-1:
                self.nbinsextra = N-1-self.start_bin
            self.bandpass_corr = False
                
