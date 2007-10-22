import sifting

# Basic parameters
# institution is one of: 'UBC', 'NRAOCV', "McGill", "Columbia", "Cornell"
institution           = "NRAOCV" 
scripts_loc           = "/homes/borgii/alfa/pipeline/PALFA_scripts/"
base_output_directory = "/home/sransom/results/ALFA"
db_pointing_file      = "/home/sransom/results/ALFA/PALFA_coords_table.txt"
pmw_module            = "/homes/borgii/alfa/processing/plotting" 
final_dir             = "/data/data6/PALFA/results/"  # this directory must exist
log_file              = "/data/data6/PALFA/results/PALFA_log.txt"
splash_image          = "/homes/borgii/alfa/pipeline/PALFA_scripts/mrpulsarplotter.gif"
known_radius          = 10 #arcmin

# Tunable parameters for searching and folding
# (you probably don't need to tune any of them)
rfifind_chunk_time    = 2**15 * 0.000064  # ~2.1 sec for dt = 64us 
singlepulse_threshold = 5.0  # threshold SNR for candidate determination
singlepulse_plot_SNR  = 6.0  # threshold SNR for singlepulse plot
singlepulse_maxwidth  = 0.1  # max pulse width in seconds
to_prepfold_sigma     = 6.0  # incoherent sum significance to fold candidates
max_cands_to_fold     = 50   # Never fold more than this many candidates
numhits_to_fold       = 2    # Number of DMs with a detection needed to fold
low_DM_cutoff         = 2.0  # Lowest DM to consider as a "real" pulsar
lo_accel_numharm      = 16   # max harmonics
lo_accel_sigma        = 2.0  # threshold gaussian significance
lo_accel_zmax         = 0    # bins
lo_accel_flo          = 2.0  # Hz
hi_accel_numharm      = 8    # max harmonics
hi_accel_sigma        = 3.0  # threshold gaussian significance
hi_accel_zmax         = 50   # bins
hi_accel_flo          = 1.0  # Hz
low_T_to_search       = 20.0 # sec

# Sifting specific parameters (don't touch without good reason!)
sifting.sigma_threshold = to_prepfold_sigma-1.0  # incoherent power threshold (sigma)
sifting.c_pow_threshold = 100.0                  # coherent power threshold
sifting.r_err           = 1.1    # Fourier bin tolerence for candidate equivalence
sifting.short_period    = 0.0005 # Shortest period candidates to consider (s)
sifting.long_period     = 15.0   # Longest period candidates to consider (s)
sifting.harm_pow_cutoff = 8.0    # Power required in at least one harmonic

