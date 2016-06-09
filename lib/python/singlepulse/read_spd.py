import numpy as _np
from scipy.special import erf as _erf

class spd:
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
        print "filename:", self.filename

        self.telescope = ll[1]
        print "telescope:", self.telescope
        self.ra = ll[2]
        print "RA:", self.ra
        self.dec = ll[3]
        print "DEC:", self.dec
        self.ra_deg = _np.sum(_np.array(self.ra.split(":"), dtype=float) * _np.array([15., 15./60., 15./3600.]))
        dec_arr = _np.array(self.dec.split(":"), dtype=float)
        self.dec_deg = _np.sum(_np.abs(dec_arr) * _np.sign(dec_arr[0]) * _np.array([1., 1./60., 1./3600.]))
        self.mjd = float(ll[4])
        print "MJD:", self.mjd
        self.total_obs_time = float(ll[15])
        print "Total_obs_time:", self.total_obs_time

        if ll[5]:
            self.rank = int(ll[5])
        else:
            self.rank=None
        print "rank:", self.rank
        self.tsamp = float(ll[14])
        print "tsamp:", self.tsamp
        self.best_dm = float(ll[8])
        print "DM:", self.best_dm
        if ll[9]:
            self.sigma = float(ll[9])
        else:
            self.sigma = None
        print "sigma:", self.sigma
        self.pulse_peak_sample = int(ll[10])
        print "pulse_peak_sample:", self.pulse_peak_sample
        self.pulse_peak_time = float(ll[16])
        print "p_p_time:", self.pulse_peak_time
        self.pulsewidth_bins = int(ll[12])
        print "pw_bins:", self.pulsewidth_bins
        self.pulsewidth_seconds = float(ll[13])
        print "pw_secs:", self.pulsewidth_seconds
        self.nsamp = int(ll[7])
        print "nsamp:", self.nsamp
        self.waterfall_duration = float(ll[11])
        print "waterfall_duration:", self.waterfall_duration
        self.waterfall_start_time = float(ll[17])
        print "waterfall_start_time:", self.waterfall_start_time
        self.waterfall_tsamp = float(ll[18])
        print "waterfall_tsamp:", self.waterfall_tsamp
        self.waterfall_nbins = self.data_zerodm_dedisp.shape[1]
        self.waterfall_nsubs = int(ll[6])
        print "waterfall_nsubs:", self.waterfall_nsubs
        self.waterfall_prededisp_nbins = int(ll[19])
        print "waterfall_prededisp_nbins:", self.waterfall_prededisp_nbins
        self.waterfall_downsamp = int(_np.round(self.waterfall_tsamp/self.tsamp))
        self.min_freq = float(ll[20])
        print "min_freq:", self.min_freq
        self.max_freq = float(ll[21])
        print "max_freq:", self.max_freq
        self.sweep_duration = float(ll[22])
        print "sweep_duration:", self.sweep_duration
        self.sweep_start_time = float(ll[23])
        print "sweep_start_time:", self.sweep_start_time
        self.bary_pulse_peak_time = float(ll[24])
        print "bary_pulse_peak_time:", self.bary_pulse_peak_time
        self.man_params = ll[25]
        print "man_params:", self.man_params

        #self.spfiles = dd['singlepulse_files']
        if not self.man_params:
            self.dmVt_this_dms = dd['dm_arr']
            self.dmVt_this_times = _np.array(dd['time_list'])
            self.dmVt_this_sigmas = dd['sigma_arr']
            self.dmVt_this_widths = dd['width_arr']
        
        # Get variance from the half of the waterfall plot that definitely should not contain the pulse
        # (which is 1/4 of the way into the plot)
        self.varprof = _np.var(self.data_zerodm_dedisp.sum(axis=0)[(self.waterfall_nbins/2):])

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
