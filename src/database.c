#include "presto.h"

static char num[41][5] =
{"0th", "1st", "2nd", "3rd", "4th", "5th", "6th", \
 "7th", "8th", "9th", "10th", "11th", "12th", \
 "13th", "14th", "15th", "16th", "17th", "18th", \
 "19th", "20th", "21st", "22nd", "23rd", "24th", \
 "25th", "26th", "27th", "28th", "29th", "30th", \
 "31st", "32nd", "33rd", "34th", "35th", "36th", \
 "37th", "38th", "39th", "40th"};


int num_psrs_in_database(void)
/* Returns the number of entries in the database */
{
  FILE *database;
  char databasenm[200], path[150] = DATABASE;
  int np, fsize;

  /* Open the binary data file */
  sprintf(databasenm, "%s/bincat.dat", path);
  database = chkfopen(databasenm, "rb");

  /* This is the length of the binary file in bytes */
  chkfread(&fsize, sizeof(int), 1, database);

  /* Go to EOF - 4 bytes and read the number of psrs */
  chkfileseek(database, (int) (fsize - sizeof(int)), 1, SEEK_CUR);
  chkfread(&np, sizeof(int), 1, database);

  /* Close the database */
  fclose(database);

  return np;
}

int read_database(psrdatabase * pdata)
/* Reads the full pulsar database into the structure *pdata */
{
  FILE *database;
  char databasenm[200], path[150] = DATABASE;
  int np, fsize;

  /* Open the binary data file */
  sprintf(databasenm, "%s/bincat.dat", path);
  database = chkfopen(databasenm, "rb");

  /* This is the length of the binary file in bytes */
  chkfread(&fsize, sizeof(int), 1, database);

  /* Read the pulsar data and the number of pulsars */
  chkfread(pdata, sizeof(*pdata), 1, database);
  chkfread(&np, sizeof(int), 1, database);

  /* Close the database */
  fclose(database);

  return np;
}


void collect_psrdata(psrdata **psrs, int *np)
/* Read the full pulsar database and place it in an array    */
/* of psrdata structures.  Return the number of psrs stored. */
{
  int i;
  psrdatabase pdata;

  /* Read the database. */
  *np = read_database(&pdata);

  /* Allocate the array or psrdata structures. */
  *psrs = (psrdata *)malloc(*np * sizeof(psrdata));

  /* Fill the structure with data */
  for(i=0; i<*np; i++) get_psr(i, *psrs+i, &pdata);
}


void get_psrdata_by_num(psrdata * psr, int pnum)
/* Read a full pulsar database entry for the pulsar pnum */
/* in the database.  Return the data in psrdata.         */
{
  static int np = 0;
  static psrdatabase pdata;

  /* Read the database */
  if (!np) {
    np = read_database(&pdata);
  }

  /* Fill the structure with data */
  if (pnum >= 0) get_psr(pnum, psr, &pdata);
  else {
    printf("Could not find the PSR in the database in get_psrdata().\n");
    exit(2);
  }
}


void get_psrdata(psrdata * psr, char * psrname)
/* Read a full pulsar database entry for pulsar psrname. */
/* Return the data in a psrdata structure.               */
{
  int np = 0, pnum = 0;
  psrdatabase pdata;

  /* Read the database. */
  np = read_database(&pdata);

  /* Find the pulsar of interest */
  if (np) pnum = psr_number_from_name(psrname, &pdata);
  else {
    printf("Problem reading the pulsar database in get_psrdata().\n");
    exit(1);
  }

  /* Fill the structure with data */
  if (pnum >= 0) get_psr(pnum, psr, &pdata);
  else {
    printf("Could not find the PSR in the database in get_psrdata().\n");
    exit(2);
  }
}


void get_psr(int psrnumber, psrdata *psr, psrdatabase * pdata)
/* Returns a full database entry for the pulsar #psrnumber in */
/* the database variable pdata.  Returns *psr completed.      */
{
  int i = psrnumber;
  
  if (i<0){
    printf("psrnumber < 0 in get_psr().  Exiting.\n\n");
    exit(1);
  }
  sprintf(psr->jname, "%.12s", pdata->jname + i * 12);
  sprintf(psr->bname, "%.8s", pdata->bname + i * 8);
  psr->ra2000 = pdata->ra2000[i];     
  psr->ra1950 = pdata->ra1950[i];     
  psr->rae = pdata->rae[i];        
  psr->dec2000 = pdata->dec2000[i];    
  psr->dec1950 = pdata->dec1950[i];    
  psr->dece = pdata->dece[i];       
  psr->nscode = pdata->nscode[i];     
  psr->dmin__ = pdata->dmin__[i];     
  psr->dmax__ = pdata->dmax__[i];     
  psr->dist = pdata->dist[i];       
  psr->ndflag = pdata->ndflag[i];     
  psr->lcode = pdata->lcode[i];      
  psr->ucode = pdata->ucode[i];      
  psr->ldeg = pdata->ldeg[i];       
  psr->bdeg = pdata->bdeg[i];       
  psr->pmra = pdata->pmra[i];       
  psr->pmrae = pdata->pmrae[i];      
  psr->pmdec = pdata->pmdec[i];      
  psr->pmdece = pdata->pmdece[i];     
  psr->posepoch = pdata->posepoch[i];   
  psr->p = pdata->p[i];          
  psr->pe = pdata->pe[i];         
  psr->pdot = pdata->pdot[i];       
  psr->pdote = pdata->pdote[i];      
  psr->f2 = pdata->f2[i];         
  psr->f2e = pdata->f2e[i];        
  psr->f3 = pdata->f3[i];         
  psr->f3e = pdata->f3e[i];        
  psr->epoch = pdata->epoch[i];      
  psr->dm = pdata->dm[i];         
  psr->dme = pdata->dme[i];        
  psr->rm = pdata->rm[i];         
  psr->rme = pdata->rme[i];        
  psr->we = pdata->we[i];         
  psr->w50 = pdata->w50[i];        
  psr->w10 = pdata->w10[i];        
  psr->s400 = pdata->s400[i];       
  psr->s600 = pdata->s600[i];       
  psr->s1400 = pdata->s1400[i];      
  psr->tau = pdata->tau[i];        
  psr->ntauflag = pdata->ntauflag[i];   
  psr->t408 = pdata->t408[i];       
  psr->ntype = pdata->ntype[i];      
  psr->modcode = pdata->modcode[i];    
  psr->limcode = pdata->limcode[i];    
  psr->distmod = pdata->distmod[i];    
  psr->lum = pdata->lum[i];        
  psr->bsurf = pdata->bsurf[i];      
  psr->age = pdata->age[i];        
  psr->edot = pdata->edot[i];       
  psr->ibin = pdata->ibin[i];       
  psr->pb = pdata->pb[i];         
  psr->pbe = pdata->pbe[i];        
  psr->a1 = pdata->a1[i];         
  psr->a1e = pdata->a1e[i];        
  psr->om = pdata->om[i];         
  psr->ome = pdata->ome[i];        
  psr->omdot = pdata->omdot[i];      
  psr->omdote = pdata->omdote[i];     
  psr->e = pdata->e[i];          
  psr->ee = pdata->ee[i];         
  psr->t0 = pdata->t0[i];         
  psr->t0e = pdata->t0e[i];        
  psr->gamma = pdata->gamma[i];      
  psr->gammae = pdata->gammae[i];     
  psr->pbdot = pdata->pbdot[i];      
  psr->pbdote = pdata->pbdote[i];     
  psr->si = pdata->si[i];         
  psr->sie = pdata->sie[i];        
  psr->r__ = pdata->r__[i];        
  psr->re = pdata->re[i];         
  psr->pb2 = pdata->pb2[i];        
  psr->pb2e = pdata->pb2e[i];       
  psr->a12 = pdata->a12[i];        
  psr->a12e = pdata->a12e[i];       
  psr->om2 = pdata->om2[i];        
  psr->om2e = pdata->om2e[i];       
  psr->omdot2 = pdata->omdot2[i];     
  psr->omdot2e = pdata->omdot2e[i];    
  psr->e2 = pdata->e2[i];         
  psr->e2e = pdata->e2e[i];        
  psr->t02 = pdata->t02[i];        
  psr->t02e = pdata->t02e[i];       
  psr->gamma2 = pdata->gamma2[i];     
  psr->gamma2e = pdata->gamma2e[i];    
  psr->pbdot2 = pdata->pbdot2[i];     
  psr->pbdot2e = pdata->pbdot2e[i];    
  psr->si2 = pdata->si2[i];        
  psr->si2e = pdata->si2e[i];       
  psr->r2 = pdata->r2[i];         
  psr->r2e = pdata->r2e[i];        
}


int psr_number_from_name(char * psrname, psrdatabase * pdata)
/* Returns the pulsar number of psrname from the database */
/* This number can be from zero to the total number       */
/* of pulsars minus 1.  This way you can use this number  */
/* as an index from the result of collect_psrdata().      */
/* Return -1 if no pulsar is found.                       */
{
  int i, psrnumber = -1;
  char psrjname[30], psrbname[30], tmp1[30], tmp2[30];

  sprintf(psrjname, "%-12s", psrname);				
  sprintf(psrbname, "%-8s", psrname);				
  
  /* Search for either the J-name or the B-name */		
  
  for (i = 0; i < NP; i++) {					
    sprintf(tmp1, "%.12s", pdata->jname + i * 12);
    sprintf(tmp2, "%.8s", pdata->bname + i * 8);
    if (!strcmp(tmp1, psrjname) || !strcmp(tmp2, psrbname)) {	
      psrnumber = i;
      break;
    }
  }
  
  /* Return the pulsar number */

  return psrnumber;
}


int return_psrparams_at_epoch(psrparams * psr, char * psrname, \
			      double epoch)
/* Reads info from the pulsar database and converts returned values */
/* to epoch 'epoch'.  Returned values go in psr (psrparams).        */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */
{
  static int np = 0;
  static psrdatabase pdata;

  /* If calling for the first time, read the database. */

  if (!np) {
    np = read_database(&pdata);
  }

  return get_psr_at_epoch(psrname, epoch, &pdata, psr);
}


int get_psr_at_epoch(char * psrname, double epoch, psrdatabase * pdata, \
		     psrparams * psr)
/* Converts info from the pulsar database to "current" epoch.       */
/* Returned values go in *psr.  The database data is in *pdata.     */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */
{
  int i;
  double difft, f, fd;

  /* Get the pulsar's number from the database */

  i = psr_number_from_name(psrname, pdata);

  /* Pulsar not in the database. */

  if (i < 0) return 0;

  /* Pulsar is in the database. */

  else {
    sprintf(psr->jname, "PSR J%.12s", pdata->jname + i * 12);
    sprintf(psr->bname, "PSR B%.8s", pdata->bname + i * 8);
    psr->ntype = pdata->ntype[i];
    psr->ra2000 = pdata->ra2000[i];
    psr->dec2000 = pdata->dec2000[i];
    psr->dm = pdata->dm[i];
    psr->dist = pdata->dist[i];
    psr->fwhm = pdata->w50[i];
    psr->timepoch = pdata->epoch[i];

    difft = SECPERDAY * (epoch - psr->timepoch);

    f = 1.0 / pdata->p[i];
    fd =  - pdata->pdot[i] * f * f;
    psr->fdd = pdata->f2[i];
    psr->f = f + fd * difft + 0.5 * psr->fdd * difft * difft;
    psr->fd = fd + psr->fdd * difft;
    psr->p = 1.0 / psr->f;
    psr->pd = - psr->fd * psr->p * psr->p;
    psr->pdd = (2.0 * (fd * fd) / f - psr->fdd) / (f * f);
    
    /* Selected pulsar is binary... */				
    
    if (psr->ntype & 8) {					
      psr->orb.t = pdata->t0[i];
      difft = SECPERDAY * (epoch - psr->orb.t);
      psr->orb.pd = pdata->pbdot[i] * 1.0E-12;
      psr->orb.p = pdata->pb[i] * SECPERDAY + psr->orb.pd * difft;
      psr->orb.x = pdata->a1[i];
      psr->orb.e = pdata->e[i];
      /* psr->orb.t is in seconds, _not_ MJD.  It represents the time    */
      /*      in sec _since_ the last periastron passage, _not_ when the */
      /*      next periastron will occur....                             */
      psr->orb.t = fmod(difft, psr->orb.p);
      if (psr->orb.t < 0.0)
	psr->orb.t += psr->orb.p;
      psr->orb.w = fmod((pdata->om[i] + difft * pdata->omdot[i] / \
			 SECPERJULYR), 360.0);
      psr->orb.wd = pdata->omdot[i];
    } else {
      psr->orb.p = 0.0;
      psr->orb.pd = 0.0;
      psr->orb.x = 0.0;
      psr->orb.e = 0.0;
      psr->orb.w = 0.0;
      psr->orb.t = 0.0;
    }
  }

  /* Return the number of the pulsar in the database. */

  return i;
}


int comp_psr_to_cand(fourierprops * cand, infodata * idata, \
		     char *output, int full)
{
  /* Compares a pulsar candidate defined by its properties found in   */
  /*   *cand, and *idata with all of the pulsars in the pulsar        */
  /*   database.  It returns a string (verbose if full==1) describing */
  /*   the results of the search in *output.                          */

  int i, j;
  static int np=0;
  static psrdatabase pdata;
  static infodata *old_idata;
  double theor, theoz, sidedr=20.0;
  double r_criteria, z_criteria, rdiff, zdiff, difft=0;
  static double T, beam2, ra, dec, epoch;
  char tmp1[80], tmp2[80], tmp3[80], shortout[30], psrname[20];
  rzwerrs rzws;

  /* If calling for the first time, read the database. */

  if (!np) {
    np = read_database(&pdata);
  }
  /* If calling for the first time for a certain data set, */
  /* initialize some values.                               */

  if (idata != old_idata) {
    /* Convert the beam width to radians */

    beam2 = 2.0 * ARCSEC2RAD * idata->fov;

    /* Convert RA and DEC to radians  (Use J2000) */

    ra = hms2rad(idata->ra_h, idata->ra_m, idata->ra_s);
    dec = dms2rad(idata->dec_d, idata->dec_m, idata->dec_s);
    T = idata->N * idata->dt;
    epoch = (double) idata->mjd_i + idata->mjd_f + T / (2.0 * SECPERDAY);

    /* Set up old_idata for next time */

    old_idata = idata;
  }
  /* Calculate the measured r, z, w's and their derivatives */

  calc_rzwerrs(cand, T, &rzws);

  /* Run through RAs in database looking for things close  */
  /* If find one, check the DEC as well (the angle between */
  /* the sources < 2*beam diam).  If find one, check its   */
  /* period.  If this matches within 2*perr, return the    */
  /* number of the pulsar.  If no matches, return 0.       */

  for (i = 0; i < np; i++) {

    /* See if we're close in RA */

    if (fabs(pdata.ra2000[i] - ra) < 5 * beam2) {

      /* See if we're close in RA and DEC */

      if (sphere_ang_diff(pdata.ra2000[i], pdata.dec2000[i], \
			  ra, dec) < 5 * beam2) {

	/* Predict the period of the pulsar at the observation MJD */

	difft = SECPERDAY * (epoch - pdata.epoch[i]);
	theor = T / (pdata.p[i] + pdata.pdot[i] * difft);
	theoz = -pdata.pdot[i] * theor * theor;

	/* Check the predicted period and its harmonics against the */
	/* measured period.                                         */

	for (j = 1; j < 41; j++) {

	  /* If the psr from the database is in a             */
	  /* binary orbit, loosen the match criteria.         */
	  /* This accounts for Doppler variations in period.  */

	  if (pdata.ntype[i] & 8) {
	    r_criteria = 0.001 * theor * j;  /* 0.1% fractional error   */
	    z_criteria = 9999999999.0;	     /* Always match for binary */
	    strcpy(tmp1, "?");
	    if (full) {
	      strcpy(tmp3, "Possibly (large error) ");
	    }
	  } else {
	    r_criteria = 5.0;                /* 5 bin error matching... */
	    z_criteria = 9999999999.0;	     /* Always match for binary */
	    /* z_criteria = 2.5 * cand->zerr; */
	    strcpy(tmp1, "");
	    if (full) {
	      strcpy(tmp3, "Looks like ");
	    }
	  }

	  if (theor*j > 1.5*cand->r)
	    break;

	  rdiff = fabs(theor * j - cand->r);
	  zdiff = fabs(theoz * j - cand->z);

	  if (rdiff < r_criteria && zdiff < z_criteria) {
	    sprintf(tmp2, "%.8s", pdata.bname + i * 8);
	    if (strcmp("        ", tmp2)==0){
	      sprintf(tmp2, "%.12s", pdata.jname + i * 12);
	      sprintf(psrname, "J%s", remove_whitespace(tmp2));
	    } else {
	      sprintf(psrname, "B%s", remove_whitespace(tmp2));
	    }
	    if (j == 1) {
	      if (full) {
		sprintf(tmp1, "the fundamental of ");
		sprintf(tmp2, "PSR %s. (predicted p = %11.7f s).\n", \
			psrname, T/theor);
		sprintf(output, "%s%s\n     %s", tmp3, tmp1, tmp2);
	      } else {
		sprintf(shortout, "PSR %s%s", psrname, tmp1);
		strncpy(output, shortout, 20);
	      }
	    } else {
	      if (full) {
		sprintf(tmp1, "the %s harmonic of ", num[j]);
		sprintf(tmp2, "PSR %s. (predicted p = %11.7f s).\n", \
			psrname, T/theor);
		sprintf(output, "%s%s\n     %s", tmp3, tmp1, tmp2);
	      } else {
		sprintf(shortout, "%s H %s%s", num[j], psrname, tmp1);
		strncpy(output, shortout, 20);
	      }
	    }
	    return i + 1;
	  } else if (rdiff < sidedr) {
	    sprintf(tmp2, "%.8s", pdata.bname + i * 8);
	    if (strcmp("        ", tmp2)==0){
	      sprintf(tmp2, "%.12s", pdata.jname + i * 12);
	      sprintf(psrname, "J%s", remove_whitespace(tmp2));
	    } else {
	      sprintf(psrname, "B%s", remove_whitespace(tmp2));
	    }
	    if (full) {
	      sprintf(tmp1, "a sidelobe of the %s harmonic of ", num[j]);
	      sprintf(tmp2, "PSR %s. (predicted p = %11.7f s).\n", \
		      psrname, T/theor);
	      sprintf(output, "%s%s\n     %s", tmp3, tmp1, tmp2);
	    } else {
	      sprintf(shortout, "SL H%d %s", j, psrname);
	      strncpy(output, shortout, 20);
	    }
	    return i + 1;
	  }
	}
      }
    }
  }

  /* Didn't find a match */

  if (full) {
    sprintf(output, \
	    "I don't recognize this candidate in the pulsar database.\n");
  } else {
    strncpy(output, "                       ", 20);
  }
  return 0;
}


int comp_bin_to_cand(binaryprops * cand, infodata * idata, \
		     char *output, int full)
  /* Compares a binary PSR candidate defined by its props found in    */
  /*   *cand, and *idata with all of the pulsars in the pulsar        */
  /*   database.  It returns a string (verbose if full==1) describing */
  /*   the results of the search in *output.                          */
{
  int i, j, k;
  static int np = 0;
  static psrdatabase pdata;
  double T, theop, ra, dec, beam2, difft = 0.0, epoch;
  double bmod, pmod;
  char tmp1[80], tmp2[80], tmp3[80], psrname[20], shortout[30];

  /* If calling for the first time, read the database. */

  if (!np) {
    np = read_database(&pdata);
  }
  /* Convert the beam width to radians */

  beam2 = 2.0 * ARCSEC2RAD * idata->fov;

  /* Convert RA and DEC to radians  (Use J2000) */

  ra = hms2rad(idata->ra_h, idata->ra_m, idata->ra_s);
  dec = dms2rad(idata->dec_d, idata->dec_m, idata->dec_s);

  /* Calculate the time related variables  */

  T = idata->N * idata->dt;
  epoch = (double) idata->mjd_i + idata->mjd_f;

  /* Run through RAs in database looking for things close  */
  /* If find one, check the DEC as well (the angle between */
  /* the sources < 2*beam diam).  If find one, check its   */
  /* period.  If this matches within 2*perr, return the    */
  /* number of the pulsar.  If no matches, return 0.       */

  for (i = 0; i < np; i++) {

    /* See if we're close in RA */

    if (fabs(pdata.ra2000[i] - ra) < beam2) {

      /* See if we're close in RA and DEC */

      if (sphere_ang_diff(pdata.ra2000[i], pdata.dec2000[i], \
			  ra, dec) < beam2) {

	/* Check that the psr in the database is in a binary   */

	if (pdata.ntype[i] & 8) {

	  /* Predict the period of the pulsar at the observation MJD */

	  difft = SECPERDAY * (epoch - pdata.epoch[i]);
	  theop = pdata.p[i] + pdata.pdot[i] * difft;

	  /* Check the predicted period and its harmonics against the */
	  /* measured period.  Use both pulsar and binary periods.    */

	  for (j = 1, pmod = 1.0; j < 41; j++, pmod = 1.0 / (double) j) {
	    if (fabs(theop * pmod - cand->ppsr) < \
		(4 * cand->ppsrerr)) {
	      for (k = 1, bmod = 1.0; k < 10; \
		   k++, bmod = 1.0 / (double) k) {
		if (fabs(pdata.pb[i] * bmod - cand->pbin / SECPERDAY) < \
		    (4 * cand->pbinerr / SECPERDAY)) {
		  sprintf(tmp2, "%.8s", pdata.bname + i * 8);
		  if (strcmp("        ", tmp2)==0){
		    sprintf(tmp2, "%.12s", pdata.jname + i * 12);
		    sprintf(psrname, "J%s", remove_whitespace(tmp2));
		  } else {
		    sprintf(psrname, "B%s", remove_whitespace(tmp2));
		  }
		  if (j > 1) {
		    sprintf(tmp1, "Possibly the %s phasemod harmonic ", num[k]);
		    if (full) {
		      sprintf(tmp2, "of the %s harmonic of PSR ", num[j]);
		      sprintf(tmp3, "%s (p = %11.7f s, pbin = %9.4f d).\n", \
			      psrname, theop, pdata.pb[i]);
		      sprintf(output, "%s%s%s", tmp1, tmp2, tmp3);
		    } else {
		      sprintf(shortout, "%s H %s", num[k], psrname);
		      strncpy(output, shortout, 20);
		    }
		  } else {
		    if (full) {
		      sprintf(tmp2, "of PSR %s (p = %11.7f s, pbin = %9.4f d).\n", \
			      psrname, theop, pdata.pb[i]);
		      sprintf(output, "%s%s", tmp1, tmp2);
		    } else {
		      sprintf(shortout, "PSR %s", psrname);
		      strncpy(output, shortout, 20);
		    }
		  }
		}
		return i + 1;
	      }
	    }
	  }
	}
      }
    }
  }

  /* Didn't find a match */

  if (full) {
    sprintf(output, "I don't recognize this candidate in the pulsar database.\n");
  } else {
    strncpy(output, "                       ", 20);
  }
  return 0;
}
