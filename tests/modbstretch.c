#include <stdio.h>
#include <math.h>
#include "../include/fftapps.h"

main(int argc, char *argv[]) {

/* bstretch: moves topocentric time series to barycentric */
/* call 'bstretch raw_datafile out_datafile timefile dt nskip' */
/* dt must be given in compressed units (e.g., 20 kHz, but c=3, dt=0.0004) */

  int i,j,k,ct2,ct3,*idat,*ivector(),ict,comp,fct,nsk,kct,ifl;
  long rct;
  float *dat,*dout,*vector();
  double dt,dtb,fact,t1,t2,tb0,tb1,tb2;
  double a1,a2,tel,telb1,telb2,flim;
  FILE *fin,*fout,*bary;
  char nam[50];

  dout=gen_fvect(16384);
  dat=gen_fvect(8192);
  idat=gen_ivect(2048);
  ifl=0;

  flim=pow(2.0,31);

  fin=fopen(argv[1],"r");
  sprintf(nam,"%s/junk.000",argv[2]);
  printf("opening %s\n",nam);
  fout=fopen(nam,"w");

  bary=fopen(argv[3],"r");

  dt=atof(argv[4])/(3600.0*24.0);
  nsk=atoi(argv[5]);
  kct=0;
  rct=0;
  for (i=0;i<16384;++i)
    dout[i]=0.0;

  /* tb0 is barycentric time - double precision */
  fscanf(bary,"%lf",&tb0);
  tb2=tb0;
  t1=0.0;
  t2=dt*8192.0;
  printf("initial time = %lf\n",tb0);


/*** MAIN LOOP ****/

  ict=0;
  fct=0;
  while(fread(dat,sizeof(float),8192,fin)==8192) {

    /* get new dt-barycentered and time offset */
    tb1=tb2;
    fscanf(bary,"%lf",&tb2);
    telb2=tb2;
    telb1=tb1;

    dtb=(telb2-telb1)/8192.0;

    if ((ict%100)==0)
      printf("time elapsed= %d %17.10f %17.10f\n",ict,telb2,dtb);
    ++ict;

    /* loop over new data segment */
    for (i=0;i<8192;++i) {

      tel=i*dtb+telb1;

      if (tel>=t2) {
	if (kct>nsk)
	  fwrite(dout,sizeof(float),8192,fout);
	++rct;

	if (rct>=65472 && ifl==0) {
	  fclose(fout);
	  sprintf(nam,"%s/junk.001",argv[2]);
	  printf("opening %s\n",nam);
	  fout=fopen(nam,"w");
	  ifl=1;
	}
	for (j=0;j<8192;++j)
	  dout[j]=dout[j+8192];
	for (j=8192;j<16384;++j)
	  dout[j]=0.0;
	t1=t2;
	t2+=dt*8192;
      }

      ++kct;
      j=(int)((tel-t1)/dt);
      k=(int)((tel-t1+dtb)/dt);
      /*      printf("%lf %lf %d %d\n",(tel-t1)/dt,(tel+dtb-t1)/dt,j,k);*/

      if (j==k)
	dout[j]+=dat[i];

      if ((k-j)==1) {

	fact=k*dt+t1-tel;
	dout[j]+=dat[i]*fact/dtb;
	dout[k]+=dat[i]*(1.0-fact/dtb);
      }

      if ((k-j)==2) {
/*	printf("hi there %d\n",ict);*/

	fact=(j+1)*dt+t1-tel;
	dout[j]+=dat[i]*fact/dtb;
	dout[j+1]+=dat[i]*dt/dtb;
	dout[k]+=dat[i]*(1.0-(fact+dt)/dtb);
      }
    }
/*    if (ict>2600)
      printf("loop done\n");*/
  }

  fclose(fin);
  fclose(bary);

  fwrite(dout,sizeof(float),8192,fout);

  fclose(fout);

  free(dout);
  free(dat);
  free(idat);


}


