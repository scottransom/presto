#include "prepfold.h"

void prepfold_plot(prepfoldinfo *search);
/* Make the beautiful 1 page prepfold output */
{
  int ii, jj, pdelay, pddelay;
  double N=0.0, T, foldf=0.0, foldfd=0.0, dphase;
  double delay, parttime, *pdprofs;
  float *redchi, *parttimes;

  if (search->topo.pow){  /* Used topocentric period for folding */
    dphase = search->topo.p1 / search->proflen;  
    foldf = 1.0 / search->topo.p1;
    foldfd = -search->topo.p2 / (search->topo.p1 * search->topo.p1);
  } else {                /* Used barycentric period for folding */
    dphase = search->bary.p1 / search->proflen;  
    foldf = 1.0 / search->bary.p1;
    foldfd = -search->bary.p2 / (search->bary.p1 * search->bary.p1);
  }

  /* Find out how many total points were folded */

  for (ii = 0; ii < search->npart; ii++)
    N += search->stats[ii * search->nsub].numdata;

  /* Calculate the time per part and the total observation time */

  parttime = search->stats[0].numdata * search->dt;
  T = N * search->dt;

  /* Generate the array of times */

  parttimes = gen_freqs(search->npart, 0.0, parttime);

  /* Array of profiles corrected for DM and pdot */

  pdprofs = gen_dvect(search->npart * search->proflen);

  /* Correct profiles for best DM */

  if (search->nsub > 1){
    double *ddprofs, *subbanddelays, hif, hifdelay;
    foldstats *ddstats;

    /* Array of profiles corrected for DM and their stats */

    ddprofs = gen_dvect(search->npart * search->proflen);
    ddstats = (foldstats *)malloc(search->npart * sizeof(foldstats));
    dmdelays = gen_ivect(search->nsub);

    /* Doppler corrected hi freq and its delay at best DM */

    hif = doppler(search->lofreq + (search->numchan - 1.0) * 
		  search->chan_wid, search.avgvoverc);
    hifdelay = delay_from_dm(search->bestdm, hif);
    subbanddelays = subband_delays(search->numchan, search->nsub, 
				   search->bestdm, search->lofreq, 
				   search->chan_wid, search->avgvoverc);
    for (ii = 0; ii < cmd->nsub; ii++)
      dmdelays[ii] = ((int) ((subbanddelays[ii] - hifdelay) / 
			     dphase + 0.5)) % search->proflen;
    free(subbanddelays);
    combine_subbands(search->rawfolds, search->stats, search->npart, 
		     search->nsub, search->proflen, dmdelays, 
		     ddprofs, ddstats);

    /* Correct for best P-dot and Period */
    
    for (ii = 0; ii < cmd->npart; ii++){
      profindex = ii * search->proflen;
      pdelay = (int) (((currentfd - foldfd) * parttimes[ii] * 
		       parttimes[ii]) / (2.0 * dphase) + 0.5);
      pddelay = (int) (((currentfd - foldfd) * parttimes[ii] * 
			parttimes[ii]) / (2.0 * dphase) + 0.5);
      shift_prof(ddprofs + profindex, search.proflen, pddelay, 
		 pdprofs + profindex);
    }

	  /* Search over the periods */

	  for (kk = 0; kk < numtrials; kk++){
	    pdelay = kk - (numtrials - 1) / 2;
	    combine_profs(pdprofs, ddstats, cmd->npart, search.proflen, 
			  pdelay, currentprof, &currentstats);
	    if (currentstats.redchi > beststats.redchi){
	      search.bestdm = search.dms[ii];
	      if (idata.bary){
		search.bary.p1 = search.periods[kk];
		search.bary.p2 = search.pdots[jj];
	      } else {
		search.topo.p1 = search.periods[kk];
		search.topo.p2 = search.pdots[jj];
	      }
	      beststats = currentstats;

printf("%ld %ld:  dm = %f  p = %17.15f   pd = %12.6e  reduced chi = %f\n", 
       kk, jj, search.bestdm, search.topo.p1, search.topo.p2, 
       beststats.redchi);

	    }
	  }
	}


    free(ddprofs);
    free(ddstats);
    free(dmdelays);
  }


  /* Fold the best profile to get stats about the folds */
  /* as well as the Time vs. Phase, Time vs. RedChi and */
  /* the Best Profile plots.                            */

  free(parttimes);
  free(pdprofs);
  

  /*
   *  Now plot the results
   */

  /* Open and prep our device */
  cpgopen(in->pgdev);
  cpgpap(10.25, 8.5/11.0);
  cpgpage();
  cpgiden();
  cpgsch(0.8);

  /* Time versus phase */
  cpgsvp (0.06, 0.27, 0.09, 0.68);
  cpgswin(0.0, 1.9999, 0.0, npart * 120.0);
  cpgbox ('BCNST', 0.0, 0, 'BCNST', 0.0, 0);
  cpgmtxt('B', 2.6, 0.5, 0.5, "Phase");
  cpgmtxt('L', 2.1, 0.5, 0.5, "Time (s)");
  lo_col_ind, hi_col_ind = cpgqcol();
  lo_col_ind = lo_col_ind + 2;
  cpgscir(lo_col_ind, hi_col_ind);
  l = Numeric.array([0.0, 1.0]);
  r = Numeric.array([1.0, 0.0]);
  g = Numeric.array([1.0, 0.0]);
  b = Numeric.array([1.0, 0.0]);
  cpgctab(l,r,g,b);
  /* cpgimag(dmplots, numtrials, nsub, 0, numtrials-1, 0, nsub-1, 0.0, 1.0, asarray([0.0, 2.0 / nsub, 0.0, 0.0, 120.0, 0.0])); */
  /*  Time versus Reduced chisqr */
  cpgsvp (0.27, 0.36, 0.09, 0.68);
  cpgswin(0.0, 15.0001, 0.0, npart * 120.0);
  cpgbox ('BCNST', 0.0, 0, 'BCT', 0.0, 0);
  cpgmtxt('B', 2.6, 0.5, 0.5, "Reduced \gx\u2\d");
  cpgline(dmplot*14, time);
  /* Combined best profile */
  cpgsvp (0.06, 0.27, 0.68, 0.94);
  cpgswin(0.0, 2.0, 0.0, 1.1);
  cpgbox ('BC', 0.0, 0, 'BC', 0.0, 0);
  cpgmtxt('T', 1.0, 0.5, 0.5, "Best Profile");
  cpgline(phases, fullplot);
  /* DM vs reduced chisqr */
  cpgsvp (0.43, 0.66, 0.09, 0.22);
  cpgswin(min(dms), max(dms)+0.0001, min(dmchi), max(dmchi));
  cpgbox ('BCNST', 0.0, 0, 'BCNST', 0.0, 0);
  cpgmtxt('L', 2.0, 0.5, 0.5, "Reduced \gx\u2\d");
  cpgmtxt('B', 2.6, 0.5, 0.5, "DM");
  cpgline(dms, dmchi);
  /* Plots for each subband */
  cpgsvp (0.43, 0.66, 0.3, 0.68);
  cpgswin(0.0-0.01, 1.0+0.01, 0.0, nsub+1.0);
  cpgbox("BCNST", 0.25, 2, "BNST", 0.0, 0);
  cpgmtxt("L", 2.0, 0.5, 0.5, "Sub-band");
  cpgswin(0.0-0.01, 1.0+0.01, min(freqs)-2.0, max(freqs)+2.0);
  cpgbox("", 0.2, 2, "CMST", 0.0, 0);
  cpgmtxt("R", 2.3, 0.5, 0.5, "Frequency (MHz)");
  cpgmtxt("B", 2.5, 0.5, 0.5, "Phase");
  /* P P-dot image */
  cpgsvp (0.74, 0.94, 0.09, 0.29);
  cpgswin(0.0, 1.001, 0.0, 1.0);
  cpgbox("BNT", 0.0, 0, "BNT", 0.0, 0);
  cpgmtxt("L", 2.0, 0.5, 0.5, "P-dot (s/s)");
  cpgmtxt("B", 2.6, 0.5, 0.5, "Period (ms)");
  cpgswin(0.0, 1.001, 0.0, 1.0);
  cpgbox("CMT", 0.0, 0, "CMT", 0.0, 0);
  cpgmtxt("R", 2.5, 0.5, 0.5, "F-dot (Hz/s)");
  cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (Hz)");
  /* Period vs reduced chisqr */
  cpgsvp (0.74, 0.94, 0.41, 0.51);
  cpgswin(min(dms), max(dms)+0.0001, 0.0, max(dmchi));
  cpgbox ('BCNST', 0.0, 0, 'BCMST', 0.0, 0);
  cpgmtxt('B', 2.3, 0.5, 0.5, "Period (ms)");
  cpgmtxt('R', 2.4, 0.5, 0.5, "Reduced \gx\u2\d");
  cpgline(dms, dmchi);
  /* P-dot vs reduced chisqr */
  cpgsvp (0.74, 0.94, 0.58, 0.68);
  cpgswin(min(dms), max(dms)+0.0001, 0.0, max(dmchi));
  cpgbox ('BCNST', 0.0, 0, 'BCMST', 0.0, 0);
  cpgmtxt('B', 2.3, 0.5, 0.5, "P-dot (s/s)");
  cpgmtxt('R', 2.4, 0.5, 0.5, "Reduced \gx\u2\d");
  cpgline(dms, dmchi);
  /* Add the Data Info area */
  cpgsvp (0.27, 0.53, 0.68, 0.94);
  cpgswin(-0.1, 1.00, -0.1, 1.1);
  cpgmtxt('T', 1.0, 0.5, 0.5, "File:  Ter5_Full_.raw");
  cpgsch(0.7);
  cpgtext(0.0, 1.0, "Candidate:  PSR 1744-24A");
  cpgtext(0.0, 0.9, "Telescope:  Parkes Multibeam");
  cpgtext(0.0, 0.8, "Epoch\dtopo\u = 50000.12313123123");
  cpgtext(0.0, 0.7, "Epoch\dbary\u = 50000.12313123123");
  cpgtext(0.0, 0.6, "T\dsamp\u = 0.000125");
  cpgtext(0.0, 0.5, "N\dfolded\u = 1000000000");
  cpgtext(0.0, 0.4, "Data Avg = 100.0");
  cpgtext(0.0, 0.3, "\gs\ddata\u = 10.0");
  cpgtext(0.0, 0.2, "Bins/profile = 64");
  cpgtext(0.0, 0.1, "Prof Avg = 100.0");
  cpgtext(0.0, 0.0, "\gs\dprof\u = 10.0");
  /* Add the Fold Info area */
  cpgsvp (0.53, 0.94, 0.68, 0.94);
  cpgswin(-0.05, 1.05, -0.1, 1.1);
  cpgsch(0.8);
  cpgmtxt('T', 1.0, 0.5, 0.5, "Search Information");
  cpgsch(0.7);
  cpgtext(0.0, 1.0, "   Best Fit Parameters");
  cpgtext(0.0, 0.9, "Reduced \gx\u2\d = 17.876876");
  cpgtext(0.57, 0.9, "P(Noise) = 0.0000002132");
  cpgtext(0.0, 0.8, "DM = 242.21232");
  cpgtext(0.0, 0.7, "P\dtopo\u = 11.232312424562(88)");
  cpgtext(0.57, 0.7, "P\dbary\u = 11.232312424562(88)");
  cpgtext(0.0, 0.6, "P'\dtopo\u = 1.2345e-12");
  cpgtext(0.57, 0.6, "P'\dbary\u = 1.2345e-12");
  cpgtext(0.0, 0.5, "P''\dtopo\u = 1.2345e-12");
  cpgtext(0.57, 0.5, "P''\dbary\u = 1.2345e-12");
  cpgtext(0.0, 0.3, "   Binary Parameters");
  cpgtext(0.0, 0.2, "P\dorb\u (s) = 65636.123213");
  cpgtext(0.57, 0.2, "a\d1\usin(i)/c (s) = 1.132213");
  cpgtext(0.0, 0.1, "e = 0.0000");
  cpgtext(0.57, 0.1, "\gw (rad) = 1.232456");
  cpgtext(0.0, 0.0, "T\dperi\u = 50000.123234221323");
  cpgclos();

}
