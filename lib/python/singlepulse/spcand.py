#!/usr/bin/env python

import numpy as np
import bary_and_topo

def topo_timeshift(bary_start_time, time_shift, topo):
    ind = np.where(topo == float(int(bary_start_time)/10*10))[0]
    return time_shift[ind]


class params:
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

    def read_from_file(self, params, tsamp, N, lofreq, hifreq, inffile, dedisp = False, \
                       scaleindep = None, zerodm = None, mask = None, \
                       bandpass_corr = False): 
        """
           Set up parameters based on input from the groups.txt file.
           Input: params: list of parameters (DM, S/N, Time, Sample number, downfactor)
                  tsamp: sampling time (downsampled: 65.5 us for PALFA)
                  lofreq: lowest observation frequency
                  hifreq: highest observation frequency
                  inffile: supply a .inf file (I suggest .rfifind.inf) 
        """
        self.subdm = params[0]
        self.sigma = params[1]
        self.bary_start_time = params[2]
        topo, bary = bary_and_topo.bary_to_topo(inffile)
        time_shift = bary-topo
        self.topo_start_time = self.bary_start_time - topo_timeshift(self.bary_start_time, \
                                                                     time_shift, topo)[0]
        self.sample_number = params[3]
        self.width_bins = params[4]
        sweep_posn = 0.0
        self.scaleindep = scaleindep
        self.downsamp = np.round((params[2]/self.sample_number/tsamp)).astype('int')
        self.duration = self.binratio * self.width_bins * tsamp * self.downsamp
        self.start = self.topo_start_time - (0.25 * self.duration)
        if (self.start<0.0):
            self.start = 0.0
        self.start_bin = np.round(self.start/tsamp).astype('int')
        self.pulse_width = self.width_bins*self.downsamp*tsamp
        if self.sigma < 10:
            self.nsub = 32
        elif self.sigma >= 10 and self.sigma < 15:
            self.nsub = 64
        else:
            self.nsub = 96
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
            self.start = self.start + (0.25*self.duration)
            self.start_bin = np.round(self.start/tsamp).astype('int')
            self.nbins = np.round(self.sweep_duration/tsamp).astype('int')
            self.nbinsextra = self.nbins
            if (self.start_bin+self.nbinsextra) > N-1:
                self.nbinsextra = N-1-self.start_bin
            self.bandpass_corr = False

    def manual_params(self, subdm, dm, sweep_dm, sigma, start_time, width_bins, downsamp, \
                      duration, nbins, nsub, tsamp, N, lofreq, hifreq, inffile, dedisp = False, \
                      scaleindep = None, zerodm = None, mask = False, bandpass_corr = False): 
        """
           Set up parameters based on input from the groups.txt file.
           Input: params: list of parameters (DM, S/N, Time, Sample number, downfactor)
                  tsamp: sampling time (downsampled: 65.5 us for PALFA)
                  lofreq: lowest observation frequency
                  hifreq: highest observation frequency
                  inffile: supply a .inf file (I suggest .rfifind.inf) 
        """
        self.subdm = subdm
        self.mask = mask
        self.sigma = sigma
        self.bary_start_time = start_time
        topo, bary = bary_and_topo.bary_to_topo(inffile)
        time_shift = bary-topo
        self.topo_start_time = self.bary_start_time - topo_timeshift(self.bary_start_time, \
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
        self.start = self.topo_start_time - (0.25 * self.duration)
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
            self.start = self.start + (0.25*self.duration)
            self.start_bin = np.round(self.start/tsamp).astype('int')
            self.nbinsextra = self.nbins


            if (self.start_bin+self.nbinsextra) > N-1:
                self.nbinsextra = N-1-self.start_bin
            self.bandpass_corr = False
                
