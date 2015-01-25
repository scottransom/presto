#include "presto.py"

/* #define USEMPI */

#ifdef USEMPI
#include "mpi.h"
#endif

void estimate_rz(psrparms psr, double T, double eo, double *r, double *z){
  double *e, *zs, *rs, dt = 1.0, avgr=0.0, avgz=0.0;
  int ii, numpts;

  numpts = (int)(T + 1);
  e = dorbint(eo, numpts, dt, psr->orb);
  rs = gen_dvect(numpts);
  zs = gen_dvect(numpts);
  for (ii=0 ; ii<numpts ; ii++){
    rs[ii] = e[ii];
    zs[ii] = e[ii];
  }
  E_to_z(zs, numpts, psr.p, T, psr->orb);
  E_to_p(rs, numpts, psr.p, psr->orb);
  for (ii=0 ; ii<numpts ; ii++){
    rs[ii] = T * ( 1.0 / psr.p - 1.0 / rs[ii]);
  }
  for (ii=0 ; ii<numpts ; ii++){
    avgr += rs[ii];
    avgz += zs[ii];
  }
  *r = avgr / numpts;
  *z = avgz / numpts;
  free(e);
  free(rs);
  free(zs);
}

double mass_funct2(double mp, double mc, double i){
  return pow(mc * sin(i), 3.0) / pow(mc + mp, 2.0);
}

double asini_c(double pb, double mf){
  return pow((mf * pb * pb / 8015123.37129), (1.0 / 3.0));
}

int main(int argc, char *argv[]){

  /* Admin parameters */
  char outdir[50] = "/home/ransom";
  char outnm[50] = "montebinresp";
  double pmass = 1.35;
  double cmass = 0.3;
  double ecc = 0.0;
  int orbsperpt = 10;

  /* Simulation parameters */
  int numTbyPb = 50;
  double minTbyPb = 0.01;
  double maxTbyPb = 10.0;
  int numppsr = 50;
  double minppsr = 0.0005;
  double maxppsr = 5.0;
  char ctype[3] = "WD";
  double Pb = 7200.0;
  double dt = 0.0001;
  char searchtype[10] = "ffdot";
  double maxTbyPb_ffdot = 1.0;
  double minTbyPb_sideband = 0.3;
  double fftlen_shortffts = 0.05;

  /* Normal variables */
  int ii, jj, kk, ll, itmp, myid=0, numprocs=1, psr_numbins, datalen, width;
  double dtmp, T, xb, eb, wb, tp, pwr, avgpwr, orbf, orbi, eo, tryr, tryz;
  char outfilenm[200];
  fcomplex *psr_resp;
  psrparams psr;
  FILE *outfile;

#ifdef USEMPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

  sprintf(outfilenm, "%s/%s_%d_%s_%s.out", \
	  outdir, outnm, myid, searchtype, ctype);
  printf("Output file is '%s'.\n", outfilenm);
  outfile = chkfopen(outfilenm, "w");

  TbyPb = gen_dvect(numTbyPb);
  logTbyPb = gen_dvect(numTbyPb);
  ppsr = gen_dvect(numppsr);
  logppsr = gen_dvect(numppsr);

  for (ii=0 ; ii<numTbyPb ; ii++){
    dtmp = (log(maxTbyPb) - log(minTbyPb)) / (numTbyPb - 1);
    logTbyPb[ii] = log(minTbyPb) + dtmp * ii;
    TbyPb[ii] = exp(logTbyPb[ii]);
  }
  for (ii=0 ; ii<numppsr ; ii++){
    dtmp = (log(maxppsr) - log(minppsr)) / (numppsr - 1);
    logppsr[ii] = log(minppsr) + dtmp * ii;
    ppsr[ii] = exp(logppsr[ii]);
  }

  /* The Simulation loops */

  /* Loop over T / Porb */
  for (ii=0 ; ii<numTbyPb ; ii++){
    T = Pb * TbyPb[ii];
    xb = asini_c(Pb, mass_funct2(pmass, cmass, pi / 3.0));

    /* Loop over ppsr */
    for (jj=0 ; jj<numppsr ; jj++){
      if !(jj % numprocs == myid) continue;
      else {
	pwr = 0.0;
	avgpwr = 0.0;
	if (((strcmp(searchtype, 'ffdot') == 0) &&
	     TbyPb[ii] < maxTbyPb_ffdot) ||
	    ((strcmp(searchtype, 'sideband') == 0) &&
	     TbyPb[ii] > minTbyPb_sideband) ||
	    (strcmp(searchtype, 'shortffts') == 0)){

	  /* Loop over the number of tries per point */
	  for (kk=0 ; kk<orbsperpt ; kk++){
	    if (ecc == 0.0){
	      wb = 0.0;
	      tp = kk * Pb / orbsperpt;
	    } else {
	      orbf = modf(kk / sqrt(orbsperpt), &orbi);
	      orbi /= sqrt(orbsperpt);
	      wb = orbf * 180.0;
	      tp = Pb * orbi;
	    }
	    psr.p = ppsr[jj];
	    psr.orb.p = Pb;
	    psr.orb.x = xb;
	    psr.orb.e = ecc;
	    psr.orb.t = wb;
	    psr.orb.w = tp;
	    psr_numbins = 2 * bin_resp_halfwidth(psr.p, T, psr.orb);
	    psr_resp = gen_bin_response(0.0, 1, psr.p, T, psr.orb,
					psr_numbins);

	    if (strcmp(searchtype,'ffdot')==0){
	      datalen = next2_to_n(psr_numbins * 2)
		if (datalen < 512) datalen = 512;
	      data = gen_cvect(datalen);
	      itmp = (len(data) - len(psr_resp)) / 2
	      for (ll=0 ; ll<psr_numbins ; ll++)
		data[ll+itmp] = psr_resp[ll];
	      eo = keplers_eqn(psr.orb.t, psr.orb.p, psr.orb.e, 1.0e-14)
	      estimate_rz(psr, T, eo, &tryr, &tryz);
	      tryr = tryr + datalen / 2.0;
	      width = 201;
	      ffd = ffdot_plane(data, tryr, 0.5, width, tryz, 2.0, width)
                        maxarg = argmax(spectralpower(ffd.flat))
                        peakr = ((maxarg % width) * 0.5 +
                                 int(tryr - (width * 0.5) / 2.0))
                        peakz = ((maxarg / width) * 2 +
                                 tryz - (width * 2) / 2)
                        if showplots:
                            show_ffdot_plane(data, tryr, tryz)
                            Pgplot.nextplotpage()
                        if debugout:
                            print ' peakr = '+`peakr`+' peakz = '+`peakz`
                        [pows[kk], rmax, zmax, rd] = \
                                   maximize_rz(data, peakr, peakz, \
                                               norm=1.0)
                        if debugout:
                            print ' rmax = '+`rmax`+' zmax = '+`zmax`
                    elif searchtype == 'sideband':
                        psr_pows = spectralpower(psr_resp)
                        data = zeros(next2_to_n(psr_numbins * 2), 'd')
                        data[0:psr_numbins] = psr_pows
                    if debugout:
                        print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                        print `ppsr[y]`+'  '+`pows[kk]`
                    tim = clock() - stim
                    if debugout:
                        print 'Time for this point was ',tim, ' s.'
            if debugout:
                print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                print `ppsr[y]`+'  '+`average(pows)`+'  ',
                print `max(pows)`+'  '+`min(pows)`
            file.write(`x`+'  '+`y`+'  '+`TbyPb[x]`+'  ')
            file.write(`ppsr[y]`+'  '+`average(pows)`+'  ')
            file.write(`max(pows)`+'  '+`min(pows)`+'\n')
            file.flush()
  free(TbyPb);
  free(logTbyPb);
  free(ppsr);
  free(logppsr);
#ifdef USEMPI
  MPI_Finalize();
#endif
      }
