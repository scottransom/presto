scripts = ['a2x.sh', 'binary_info.py', 'chooseN.py', 'combine_weights.py', 
  'compare_periods.py', 'dat2tim.py', 'DDplan.py', 'detrend_dat.py', 
  'downsample_filterbank.py', 'event_peak.py', 'fb_truncate.py', 
  'filter_zerolags.py', 'fit_circular_orbit.py', 'fitorb.py', 
  'GBNCC_search.py', 'GBT350_drift_prep.py', 'GBT350_drift_search.py', 
  'get_TOAs.py', 'gotocand.py', 'guppidrift2fil.py', 'GUPPI_drift_prep.py', 
  'injectpsr.py', 'make_spd.py', 'makezaplist.py', 'orbellipsefit.py', 
  'PALFA_presto_search.py', 'pfd2png.sh', 'pfd_for_timing.py', 'plot_spd.py', 
  'powerstats.py', 'psrfits2fil.py', 'psrfits_quick_bandpass.py', 
  'pulsestack.py', 'pygaussfit.py', 'pyplotres.py', 'quickffdots.py', 
  'quick_prune_cands.py', 'rfifind_stats.py', 'rrattrap.py', 
  'simple_zapbirds.py', 'single_pulse_search.py', 'sortwappfiles.py', 
  'subband_smearing.py', 'sum_profiles.py', 'tim2dat.py', 
  'waterfaller.py', 'weights_to_ignorechan.py']

# This puts all the scripts into the requested binary directory
foreach script: scripts
  install_data(script, install_dir: get_option('bindir'))
endforeach
