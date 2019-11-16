from builtins import object
import numpy as _np
from scipy.special import erf as _erf
import sys


class spd(object):
    """
    A class for reading in single pulse files.

    A quick description of each item in the class:

     data_zerodm
        A 2D freq-vs-time array around the pulse, not dedispersed (zero-DM'd data)
     data_zerodm_dedisp
        A 2D freq-vs-time array around the pulse, dedispersed (zero-DM'd data)
     data_nozerodm
        A 2D freq-vs-time array around the pulse, not dedispersed (non-zero-DM'd data)
     data_nozerodm_dedisp
        A 2D freq-vs-time array around the pulse, dedispersed (non-zero-DM'd data)
     
     spfiles
        PRESTO single pulse files for the DM-vs-time scatterplot
     dmVt_this_dms
        DM values for this pulse in the DM-vs-time scatterplot
     dmVt_this_times
        Time values for this pulse in the DM-vs-time scatterplot
     dmVt_this_sigmas
        Sigma values (determining point size) for this pulse in the DM-vs-time scatterplot
     dmVt_this_widths
        box car values for this pulse.


     dmsweep_delays
        Delays corresponding to frequencies for drawn-in dispersion sweep
     dmsweep_freqs
        Frequencies corresponding to delays for drawn-in dispersion sweep

     filename
        Name of the observation file that was analyzed     
     telescope
        Which telescope was used
     ra
        Right ascension as hh:mm:ss.s string
     dec
        Declination as dd:mm:ss.s string
     ra_deg
        Right ascension in degrees
     dec_deg
        Declination in degrees
     mjd
        Observation MJD
     total_obs_time
        Total duration of the observation this pulse was found in, in seconds
     rank
        Single pulse sifting rank
     tsamp
        Sampling time of raw data in seconds
     best_dm
        Best determined dispersion measure for this event
     sigma
        Significance of this event
     pulse_peak_sample
        The sample number in the full dedispersed time series at which this event peaked
     pulse_peak_time
        The time in seconds in the full dedispersed time series at which this event peaked
     pulsewidth_bins
        The width of the boxcar filter used to optimally detect this event, in number of bins
     pulsewidth_seconds
        The width of the boxcar filter used to optimally detect this event, in seconds
     nsamp
        The number of original time series samples included in the (possibly downsampled) waterfall plot
     waterfall_duration
        The total duration of the dedispersed waterfall plot
     waterfall_start_time
        The time (in seconds) in the full dedispersed time series at which the waterfall plot begins
     waterfall_tsamp
        Sampling time of the waterfall plot in seconds
     waterfall_nbins
        The number of samples across the dedispersed waterfall plot
     waterfall_nsubs
        The number of frequency bins across the waterfall plot
     waterfall_prededisp_nbins
        The number of samples prior to dedispersing and cutting off the ends of the waterfall plot
     min_freq
        The lowest frequency plotted
     max_freq
        The highest frequency plotted
     sweep_duration
        The total duration of the dispersed pulse across the band
     sweep_start_time
        The time at which to start plotting the dispersed reference line
    """
    def __init__(self, spd_file):
        dd = dict(_np.load(spd_file))
        self.data_zerodm = dd['Data_zerodm']
        self.data_zerodm_dedisp = dd['Data_dedisp_zerodm']
        self.data_nozerodm = dd['Data_nozerodm']
        self.data_nozerodm_dedisp = dd['Data_dedisp_nozerodm']
         

        self.dmsweep_delays = dd['delays_nozerodm']
        self.dmsweep_freqs = dd['freqs_nozerodm']

        ll = dd['text_array']
       
        self.filename = ll[0]

        self.telescope = ll[1]
        self.ra = ll[2]
        self.dec = ll[3]
        self.ra_deg = _np.sum(_np.array(self.ra.split(":"), dtype=float) * _np.array([15., 15./60., 15./3600.]))
        dec_arr = _np.array(self.dec.split(":"), dtype=float)
        self.dec_deg = _np.sum(_np.abs(dec_arr) * _np.sign(dec_arr[0]) * _np.array([1., 1./60., 1./3600.]))
        self.mjd = float(ll[4])
        self.total_obs_time = float(ll[15])

        if ll[5]:
            self.rank = int(ll[5])
        else:
            self.rank=None
        self.tsamp = float(ll[14])
        self.best_dm = float(ll[8])
        if ll[9]:
            self.sigma = float(ll[9])
        else:
            self.sigma = None
        self.pulse_peak_sample = int(ll[10])
        self.pulse_peak_time = float(ll[16])
        self.pulsewidth_bins = int(ll[12])
        self.pulsewidth_seconds = float(ll[13])
        self.nsamp = int(ll[7])
        self.waterfall_duration = float(ll[11])
        self.waterfall_start_time = float(ll[17])
        self.waterfall_tsamp = float(ll[18])
        self.waterfall_nbins = self.data_zerodm_dedisp.shape[1]
        self.waterfall_nsubs = int(ll[6])
        self.waterfall_prededisp_nbins = int(ll[19])
        self.waterfall_downsamp = int(_np.round(self.waterfall_tsamp/self.tsamp))
        self.min_freq = float(ll[20])
        self.max_freq = float(ll[21])
        self.sweep_duration = float(ll[22])
        self.sweep_start_time = float(ll[23])
        self.bary_pulse_peak_time = float(ll[24])
        self.man_params = ll[25]

        #self.spfiles = dd['singlepulse_files']
        if not self.man_params:
            self.dmVt_this_dms = dd['dm_arr']
            self.dmVt_this_times = _np.array(dd['time_list'])
            self.dmVt_this_sigmas = dd['sigma_arr']
            self.dmVt_this_widths = dd['width_arr']
        
        # Get variance from the half of the waterfall plot that definitely should not contain the pulse
        # (which is 1/4 of the way into the plot)
        self.varprof = _np.var(self.data_zerodm_dedisp.sum(axis=0)[(self.waterfall_nbins/2):])

        #### Print the spd information ####
        sys.stdout.write("\n")
        sys.stdout.write("\t\tReading Header information:     \n")
        sys.stdout.write("\n")
        sys.stdout.write("                 filename  = %s\n"%self.filename)
        sys.stdout.write("                 telescope = %s\n"%self.telescope)
        sys.stdout.write("                        RA = %s\n"%self.ra)
        sys.stdout.write("                       DEC = %s\n"%self.dec)
        sys.stdout.write("                  min_freq = %0.2f\n"%self.min_freq)
        sys.stdout.write("                  max_freq = %0.2f\n"%self.max_freq)
        sys.stdout.write("                     tsamp = %f\n"%self.tsamp)
        sys.stdout.write("            Total_obs_time = %f\n"%self.total_obs_time)
        sys.stdout.write("                       MJD = %f\n"%self.mjd)
        sys.stdout.write("                        DM = %.2f\n"%self.best_dm)
        if self.sigma:
            sys.stdout.write("                     sigma = %.2f\n"%self.sigma)
        if self.rank:
            sys.stdout.write("                      rank = %i\n"%self.rank)
        sys.stdout.write("         pulse_peak_sample = %i\n"%self.pulse_peak_sample)
        sys.stdout.write("           pulse_peak_time = %0.2f\n"%self.pulse_peak_time)
        sys.stdout.write("          pulse_width_bins = %i\n"%self.pulsewidth_bins)
        sys.stdout.write("          pulse_width_secs = %0.2f\n"%self.pulsewidth_seconds)
        sys.stdout.write("                     nsamp = %i\n"%self.nsamp)
        sys.stdout.write("        waterfall_duration = %0.2f\n"%self.waterfall_duration)
        sys.stdout.write("      waterfall_start_time = %0.2f\n"%self.waterfall_start_time)
        sys.stdout.write("           waterfall_tsamp = %0.2f\n"%self.waterfall_tsamp)
        sys.stdout.write("           waterfall_nsubs = %i\n"%self.waterfall_nsubs)
        sys.stdout.write(" waterfall_prededisp_nbins = %i\n"%self.waterfall_prededisp_nbins)
        sys.stdout.write("            sweep_duration = %0.2f\n"%self.sweep_duration)
        sys.stdout.write("          sweep_start_time = %0.2f\n"%self.sweep_start_time)
        sys.stdout.write("      bary_pulse_peak_time = %0.2f\n"%self.bary_pulse_peak_time)

    def waterfall_time_axis(self, use_timeseries_time=False):
        """
        Generate a time axis for the waterfall plot in seconds, either beginning
        at zero or at the duration into the time series at which the plot actually
        begins.
        """
        self.waterfall_tsamp
        self.waterfall_start_time
        self.waterfall_nbins
        time_axis = _np.arange(0, self.waterfall_duration, self.waterfall_tsamp)[:self.waterfall_nbins]
        if use_timeseries_time: return time_axis + self.waterfall_start_time
        else: return time_axis

    def waterfall_freq_axis(self):
        """
        Generate a frequency axis for the waterfall plot.
        """
        return _np.linspace(self.min_freq, self.max_freq, self.waterfall_nsubs, endpoint=False)

def ddm_response(ddm, width_ms, band_MHz=(1214., 1537.)):
    """
    Returns the factor by which the S/N of a pulse of a given width observed
    in a particular radio band should decrease given an error in dispersion
    measure.

    ddm
       Difference from optimal dispersion measure in pc/cm^3
    width_ms
       Pulse width in milliseconds
    band
       The bottom and top of the observing band in MHz
        (default: the Arecibo Mock band)
    """
    if _np.isscalar(ddm):
        ddm = _np.array([ddm])
        scal = True
    else:
        ddm = _np.array(ddm)
        scal = False
    band_MHz = _np.array(band_MHz)
    zeta = 6.91e-3 * ddm * _np.diff(band_MHz)[0] / (width_ms * (_np.mean(band_MHz)/1000.)**3)
    result = _np.zeros_like(ddm)
    where_nonzero = _np.where(zeta != 0)
    result[where_nonzero] = 0.5*_np.sqrt(_np.pi)*_erf(zeta[where_nonzero])/zeta[where_nonzero]
    result[zeta == 0] = 1.
    if scal: return result[0]
    else: return result
