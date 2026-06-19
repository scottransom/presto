#include <ctype.h>
#include "presto.h"
#include "mask.h"
/* The pipeline prototypes below only need Cmdline as an opaque pointer, so   */
/* forward-declare the type rather than pulling in prepfold_cmd.h.  Including  */
/* the full header here would drag struct s_Cmdline into show_pfd.c, which     */
/* also includes show_pfd_cmd.h with its own (conflicting) struct s_Cmdline.   */
/* The .c files that touch Cmdline members include prepfold_cmd.h themselves.  */
typedef struct s_Cmdline Cmdline;

/* This causes the barycentric motion to be calculated once per second */

#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

/* Flags used to tweak the plots */

typedef struct plotflags {
  int events;
  int nosearch;
  int justprofs;
  int scaleparts;
  int allgrey;
  int fixchi;
  int samples;
  int showfold;
} plotflags;

/* structure used to pass information to plotting routine */

typedef struct PREPFOLDINFO {
  double *rawfolds;   /* Raw folds (nsub * npart * proflen points) */
  double *dms;        /* DMs used in the trials */
  double *periods;    /* Periods used in the trials */
  double *pdots;      /* P-dots used in the trials */
  foldstats *stats;   /* Statistics for the raw folds */
  int numdms;         /* Number of 'dms' */
  int numperiods;     /* Number of 'periods' */
  int numpdots;       /* Number of 'pdots' */
  int nsub;           /* Number of frequency subbands folded */
  int npart;          /* Number of folds in time over integration */
  int proflen;        /* Number of bins per profile */
  int numchan;        /* Number of channels for radio data */
  int pstep;          /* Minimum period stepsize in profile phase bins */
  int pdstep;         /* Minimum p-dot stepsize in profile phase bins */
  int dmstep;         /* Minimum DM stepsize in profile phase bins */
  int ndmfact;        /* 2*ndmfact*proflen+1 DMs to search */
  int npfact;         /* 2*npfact*proflen+1 periods and p-dots to search */
  char *filenm;       /* Filename of the folded data */
  char *candnm;       /* String describing the candidate */
  char *telescope;    /* Telescope where observation took place */
  char *pgdev;        /* PGPLOT device to use */
  char rastr[16];     /* J2000 RA  string in format hh:mm:ss.ssss */
  char decstr[16];    /* J2000 DEC string in format dd:mm:ss.ssss */
  double dt;          /* Sampling interval of the data */
  double startT;      /* Fraction of observation file to start folding */
  double endT;        /* Fraction of observation file to stop folding */
  double tepoch;      /* Topocentric eopch of data in MJD */
  double bepoch;      /* Barycentric eopch of data in MJD */
  double avgvoverc;   /* Average topocentric velocity */
  double lofreq;      /* Center of low frequency radio channel */
  double chan_wid;    /* Width of each radio channel in MHz */
  double bestdm;      /* Best DM */
  position topo;      /* Best topocentric p, pd, and pdd */
  position bary;      /* Best barycentric p, pd, and pdd */
  position fold;      /* f, fd, and fdd used to fold the initial data */
  orbitparams orb;    /* Barycentric orbital parameters used in folds */
} prepfoldinfo;

/* Per-candidate folding context.  Bundles the per-candidate working state    */
/* that the fold pipeline threads through its phase functions, so that         */
/* "fold one candidate" is expressed as a single foldcand* rather than ~30     */
/* loose locals plus in-place mutation of the parsed Cmdline.  Shared,         */
/* per-observation state (spectra_info, infodata, mask, observation timing,    */
/* barycentric time tables) is deliberately NOT stored here; it stays as       */
/* separate parameters owned by the driver.                                    */
/* NB: the resolved profile length lives in search.proflen (alongside          */
/* rawfolds/stats, which prepfoldinfo already owns); there is no separate      */
/* proflen member to avoid two sources of truth.                               */
typedef struct FOLDCAND {
  prepfoldinfo search; /* per-candidate output container (rawfolds, stats,    */
                       /* proflen, nsub, npart, candnm, pgdev, ...)           */
  /* spin parameters */
  double f, fd, fdd;
  double foldf, foldfd, foldfdd;
  double orig_foldf;
  /* dispersion (per-candidate DM is the crux of the multi-candidate goal) */
  double dm;
  int nsub;
  int npart;
  int *idispdts;
  double *obsf;
  /* binary / orbit */
  int binary;
  int numdelays;
  long long numbinpoints;
  double *Ep, *tp;
  /* folding control + scratch */
  int flags;
  double phs;          /* per-block fold phase offset (was cmd->phs)          */
  double *parttimes;
  /* polyco state */
  int polyco_index;
  double polyco_phase0;
  /* optimization results */
  foldstats beststats;
  double *bestprof;
  float *ppdot;
  /* output file names (the candidate name string lives in search.candnm) */
  char *outfilenm, *plotfilenm;
} foldcand;

/* Some function definitions */

int read_resid_rec(FILE * file, double *toa, double *obsf);

int read_floats(FILE *file, float *data, int numpts, int numchan);
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */

int read_shorts(FILE *file, float *data, int numpts, int numchan);
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains short integer data.  */
/* The equivalent floats are placed in *data.               */
/* It returns the number of points read.                    */

int read_PRESTO_subbands(FILE *infiles[], int numfiles, float *subbanddata, 
                         double timeperblk, int *maskchans, 
                         int *nummasked, mask *obsmask, float *padvals);
/* Read short int subband data written by prepsubband */

double *read_events(FILE *infile, int bin, int days, int *numevents,
		    double MJD0, double Ttot, double startfrac, double endfrac, 
		    double offset);
/* This routine reads a set of events from the open file 'infile'.     */
/* It returns a double precision vector of events in seconds from the  */
/* first event.  If 'bin' is true the routine treats the data as       */
/* binary double precision (otherwise text).  If 'days' is 1 then the  */
/* data is assumed to be in days since the 'inf' EPOCH (0 is sec from  */
/* EPOCH in 'inf').  If 'days' is 2, the data are assumed to be MJDs.  */
/* The number of events read is placed in 'numevents', and the raw     */
/* event is placed in 'firstevent'.  MJD0 is the time to use for the   */
/* reference time.  Ttot is the duration of the observation. 'start'   */
/* and 'end' are define the fraction of the observation22 that we are    */
/* interested in.  'offset' is a time offset to apply to the events.   */

void hunt(double *xx, int n, double x, int *jlo);

double switch_pfdot(double pf, double pfdot);

double switch_pfdotdot(double pf, double pfdot, double pfdotdot);

double fdot2phasedelay(double fdot, double time);

double phasedelay2fdot(double phasedelay, double time);

double fdotdot2phasedelay(double fdotdot, double time);

double phasedelay2fdotdot(double phasedelay, double time);

void double2float(double *in, float *out, int numpts);
/* Copy a double vector into a float vector */

void prepfold_plot(prepfoldinfo *in, plotflags *flags, int xwin, float *ppdot);
/* Make the beautiful 1 page prepfold output */

int bary2topo(double *topotimes, double *barytimes, int numtimes, 
	      double fb, double fbd, double fbdd, 
	      double *ft, double *ftd, double *ftdd);
/* Convert a set of barycentric pulsar spin parameters (fb, fbd, fbdd) */
/* into topocentric spin parameters (ft, ftd, ftdd) by performing      */
/* a linear least-squares fit (using LAPACK routine DGELS).  The       */
/* routine equates the pulse phase using topcentric parameters and     */
/* times to the pulse phase using barycentric parameters and times.    */

void init_prepfoldinfo(prepfoldinfo *in);
/* Set all values to 0 or NULL */

void delete_prepfoldinfo(prepfoldinfo *in);
/* Free all dynamic arrays in the prepfold array */

void print_prepfoldinfo(prepfoldinfo *in);
/* Print a prepfoldinfo data structure to STDOUT */

void write_prepfoldinfo(prepfoldinfo *in, char *filename);
/* Write a prepfoldinfo data structure to a binary file */

void read_prepfoldinfo(prepfoldinfo *in, char *filename);
/* Read a prepfoldinfo data structure from a binary file */

int cpgnice_output_2(char *out, double val, double err, int len);
/* Return a string that has the same formatting as       */
/* nice_output_2(), but for PGPLOT.  This way, exponents */
/* are actually in superscript!  Woo-hoo!                */

void correct_subbands_for_DM(double dm, prepfoldinfo *search,
			     double *ddprofs, foldstats *ddstats);
/* Calculate the DM delays and apply them to the subbands */
/* to create de-disopersed profiles.                      */

void normalize_stats(double *inprofs, foldstats * stats,
                    int numparts, int numsubbands, int proflen);
// Normalize the profiles and statistics for each fold, so that the
// profile average is ~0.0 and stdev ~1.0.  For weak pulses, this means
// that the off-pulse RMS is ~1, which should be OK for most pulsars
// as long as you are folding raw data with many parts and subbands

float estimate_offpulse_redchi2(double *inprofs, foldstats *stats,
                                int numparts, int numsubbands,
                                int proflen, int numtrials, double dofeff);
// Randomly offset each pulse profile in a .pfd data square or cube
// and combine them to estimate a "true" off-pulse level.  Do this
// numtrials times in order to improve the statistics.  Return the
// inverse of the average of the off-pulse reduced-chi^2 (i.e. the
// correction factor).  dofeff is the effective number of DOF as
// returned by DOF_corr().

/* prepfold_pipeline.c:  the folding pipeline split out of prepfold.c's main() */
/* into behavior-preserving phase functions.  Shared locals are passed by      */
/* pointer (Stage 1; a foldcand context will bundle them in Stage 2).          */

/* Forward declaration so the prototypes below can take a pointer to it;       */
/* the full definition lives in backend_common.h, included by the .c files.    */
struct spectra_info;

void init_foldcand(foldcand *fc);
/* Zero a foldcand and initialize its embedded prepfoldinfo. */

void free_foldcand(foldcand *fc, Cmdline *cmd, infodata *idata);
/* Free the per-candidate allocations owned by a foldcand (mirrors the     */
/* per-candidate portion of the original cleanup_fold).                    */

void identify_and_open_input(Cmdline *cmd, struct spectra_info *s, infodata *idata,
                             foldcand *fc, mask *obsmask, plotflags *pflags,
                             char **rootnm_out, int *ptsperrec_out,
                             long long *numrec_out, int *numchan_out,
                             int *useshorts_out, int *insubs_out,
                             double **events_out, int *numevents_out, double *T_out);

void resolve_output_names(Cmdline *cmd, foldcand *fc, char *rootnm);

/* Shared explicit-spin / output-name / epoch helpers, single-sourced so that  */
/* prepfold and prepfold_multi cannot drift.                                    */
char *make_autogen_candnm(int is_period, double x0);
void build_output_filenames(foldcand *fc, char *rootnm);
void resolve_default_proflen(Cmdline *cmd, prepfoldinfo *search);
double bary_epoch_to_infinite_freq(infodata *idata, double bepoch, double dm);

void compute_obs_timing(Cmdline *cmd, struct spectra_info *s, infodata *idata,
                        foldcand *fc, mask *obsmask, int insubs,
                        int useshorts, int ptsperrec, char *obs, char *ephem,
                        char *rastring, char *decstring, int *numchan_out,
                        int **maskchans_out, double *recdt_out, long long *lorec_out,
                        double *startTday_out, long long *numrec_out,
                        long *reads_per_part_out, double *T_out, double *N_out,
                        long *worklen_out);

void resolve_fold_params(Cmdline *cmd, infodata *idata, foldcand *fc,
                         int insubs, double T, double startTday, double recdt,
                         long long lorec, char *pname);

void setup_orbit_delays(Cmdline *cmd, infodata *idata, foldcand *fc,
                        int insubs, double T, double *N_out);

void print_fold_info(Cmdline *cmd, foldcand *fc, double N, double T, char *pname);

void allocate_fold_arrays(Cmdline *cmd, foldcand *fc);

void fold_events(Cmdline *cmd, infodata *idata, foldcand *fc, double T,
                 int numevents, double *events, double startTday);

void compute_bary_corrections(Cmdline *cmd, foldcand *fc, double T,
                              char *rastring, char *decstring, char *obs,
                              char *ephem, int *numbarypts_out,
                              double **barytimes_out, double **topotimes_out);

/* compute_bary_corrections split into a candidate-INDEPENDENT table builder    */
/* (one TEMPO barycenter() call) and a per-candidate fold-frequency correction. */
/* prepfold calls them via the compute_bary_corrections wrapper; prepfold_multi */
/* calls compute_bary_table once and apply_bary_fold_correction per candidate.  */
void compute_bary_table(double tepoch, double T, char *rastring, char *decstring,
                        char *obs, char *ephem, int *numbarypts_out,
                        double **barytimes_out, double **topotimes_out,
                        double *avgvoverc_out);
void apply_bary_fold_correction(foldcand *fc, int numbarypts,
                                double *barytimes, double *topotimes);

void compute_dispersion_delays(Cmdline *cmd, infodata *idata, foldcand *fc,
                               int numchan, int insubs);

void fold_timeseries(Cmdline *cmd, struct spectra_info *s, infodata *idata,
                     foldcand *fc, mask *obsmask, double T, double N,
                     long long lorec, int ptsperrec, double recdt, long worklen,
                     long reads_per_part, int numchan, double startTday,
                     int useshorts, int insubs, int *maskchans, char *obs,
                     char *ephem, char *rastring, char *decstring, float **data_out,
                     int *numbarypts_out, double **barytimes_out,
                     double **topotimes_out);

void fold_subband_block(Cmdline *cmd, foldcand *fc, float *data, long worklen,
                        long numread, long partnum, double fold_time0,
                        double foldf, double foldfd, double foldfdd,
                        double *buffers, double *phasesadded);
/* Fold one cleaned+dedispersed raw block (fc->nsub subbands, each `worklen`  */
/* points in `data`) into this candidate's rawfolds/stats for sub-integration */
/* `partnum`.  Shared by fold_timeseries (single-candidate) and prepfold_multi */
/* (read-once/dedisperse-many) so the fold math lives in exactly one place.   */

void optimize_candidate(Cmdline *cmd, infodata *idata, foldcand *fc,
                        plotflags *pflags, double T);

void write_results_and_plot(Cmdline *cmd, infodata *idata, foldcand *fc,
                            plotflags *pflags, double N, int insubs,
                            double *barytimes, double *topotimes, int numbarypts);

void cleanup_fold(Cmdline *cmd, mask *obsmask, int insubs, float *data,
                  char *rootnm, double *barytimes, double *topotimes);
