#ifndef SLALIBHDEF
#define SLALIBHDEF
/*
**  Author:
**    Patrick Wallace  (ptw@tpsoft.demon.co.uk)
**
**  License:
**    This program is free software; you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation; either version 2 of the License, or
**    (at your option) any later version.
**
**    This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program; if not, write to the Free Software 
**    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  
**    USA.
**
**  Last revision:   10 December 2002
**
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

void slaAddet ( double rm, double dm, double eq, double *rc, double *dc );

void slaAfin ( char *string, int *iptr, float *a, int *j );

double slaAirmas ( double zd );

void slaAltaz ( double ha, double dec, double phi,
                double *az, double *azd, double *azdd,
                double *el, double *eld, double *eldd,
                double *pa, double *pad, double *padd );

void slaAmp ( double ra, double da, double date, double eq,
              double *rm, double *dm );

void slaAmpqk ( double ra, double da, double amprms[21],
                double *rm, double *dm );

void slaAop ( double rap, double dap, double date, double dut,
              double elongm, double phim, double hm, double xp,
              double yp, double tdk, double pmb, double rh,
              double wl, double tlr,
              double *aob, double *zob, double *hob,
              double *dob, double *rob );

void slaAoppa ( double date, double dut, double elongm, double phim,
                double hm, double xp, double yp, double tdk, double pmb,
                double rh, double wl, double tlr, double aoprms[14] );

void slaAoppat ( double date, double aoprms[14] );

void slaAopqk ( double rap, double dap, double aoprms[14],
                double *aob, double *zob, double *hob,
                double *dob, double *rob );

void slaAtmdsp ( double tdk, double pmb, double rh, double wl1,
                 double a1, double b1, double wl2, double *a2, double *b2 );

void slaAv2m ( float axvec[3], float rmat[3][3] );

float slaBear ( float a1, float b1, float a2, float b2 );

void slaCaf2r ( int ideg, int iamin, float asec, float *rad, int *j );

void slaCaldj ( int iy, int im, int id, double *djm, int *j );

void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j );

void slaCc2s ( float v[3], float *a, float *b );

void slaCc62s ( float v[6], float *a, float *b, float *r,
                float *ad, float *bd, float *rd );

void slaCd2tf ( int ndp, float days, char *sign, int ihmsf[4] );

void slaCldj ( int iy, int im, int id, double *djm, int *j );

void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat );

void slaCombn ( int nsel, int ncand, int list[], int *j );

void slaCr2af ( int ndp, float angle, char *sign, int idmsf[4] );

void slaCr2tf ( int ndp, float angle, char *sign, int ihmsf[4] );

void slaCs2c ( float a, float b, float v[3] );

void slaCs2c6 ( float a, float b, float r, float ad,
                float bd, float rd, float v[6] );

void slaCtf2d ( int ihour, int imin, float sec, float *days, int *j );

void slaCtf2r ( int ihour, int imin, float sec, float *rad, int *j );

void slaDaf2r ( int ideg, int iamin, double asec, double *rad, int *j );

void slaDafin ( char *string, int *iptr, double *a, int *j );

double slaDat ( double dju );

void slaDav2m ( double axvec[3], double rmat[3][3] );

double slaDbear ( double a1, double b1, double a2, double b2 );

void slaDbjin ( char *string, int *nstrt,
                double *dreslt, int *jf1, int *jf2 );

void slaDc62s ( double v[6], double *a, double *b, double *r,
                double *ad, double *bd, double *rd );

void slaDcc2s ( double v[3], double *a, double *b );

void slaDcmpf ( double coeffs[6], double *xz, double *yz, double *xs,
                double *ys, double *perp, double *orient );

void slaDcs2c ( double a, double b, double v[3] );

void slaDd2tf ( int ndp, double days, char *sign, int ihmsf[4] );

void slaDe2h ( double ha, double dec, double phi,
               double *az, double *el );

void slaDeuler ( char *order, double phi, double theta, double psi,
                 double rmat[3][3] );

void slaDfltin ( char *string, int *nstrt, double *dreslt, int *jflag );

void slaDh2e ( double az, double el, double phi, double *ha, double *dec);

void slaDimxv ( double dm[3][3], double va[3], double vb[3] );

void slaDjcal ( int ndp, double djm, int iymdf[4], int *j );

void slaDjcl ( double djm, int *iy, int *im, int *id, double *fd, int *j );

void slaDm2av ( double rmat[3][3], double axvec[3] );

void slaDmat ( int n, double *a, double *y, double *d, int *jf, int *iw );

void slaDmoon ( double date, double pv[6] );

void slaDmxm ( double a[3][3], double b[3][3], double c[3][3] );

void slaDmxv ( double dm[3][3], double va[3], double vb[3] );

double slaDpav ( double v1[3], double v2[3] );

void slaDr2af ( int ndp, double angle, char *sign, int idmsf[4] );

void slaDr2tf ( int ndp, double angle, char *sign, int ihmsf[4] );

double slaDrange ( double angle );

double slaDranrm ( double angle );

void slaDs2c6 ( double a, double b, double r, double ad, double bd,
                double rd, double v[6] );

void slaDs2tp ( double ra, double dec, double raz, double decz,
                double *xi, double *eta, int *j );

double slaDsep ( double a1, double b1, double a2, double b2 );

double slaDsepv ( double v1[3], double v2[3] );

double slaDt ( double epoch );

void slaDtf2d ( int ihour, int imin, double sec, double *days, int *j );

void slaDtf2r ( int ihour, int imin, double sec, double *rad, int *j );

void slaDtp2s ( double xi, double eta, double raz, double decz,
                double *ra, double *dec );

void slaDtp2v ( double xi, double eta, double v0[3], double v[3] );

void slaDtps2c ( double xi, double eta, double ra, double dec,
                 double *raz1, double *decz1,
                 double *raz2, double *decz2, int *n );

void slaDtpv2c ( double xi, double eta, double v[3],
                 double v01[3], double v02[3], int *n );

double slaDtt ( double dju );

void slaDv2tp ( double v[3], double v0[3], double *xi, double *eta, int *j );

double slaDvdv ( double va[3], double vb[3] );

void slaDvn ( double v[3], double uv[3], double *vm );

void slaDvxv ( double va[3], double vb[3], double vc[3] );

void slaE2h ( float ha, float dec, float phi, float *az, float *el );

void slaEarth ( int iy, int id, float fd, float posvel[6] );

void slaEcleq ( double dl, double db, double date, double *dr, double *dd );

void slaEcmat ( double date, double rmat[3][3] );

void slaEcor ( float rm, float dm, int iy, int id, float fd,
               float *rv, float *tl );

void slaEg50 ( double dr, double dd, double *dl, double *db );

void slaEl2ue ( double date, int jform, double epoch, double orbinc,
                double anode, double perih, double aorq, double e,
                double aorl, double dm, double u[], int *jstat );

double slaEpb ( double date );

double slaEpb2d ( double epb );

double slaEpco ( char k0, char k, double e );

double slaEpj ( double date );

double slaEpj2d ( double epj );

void slaEqecl ( double dr, double dd, double date, double *dl, double *db );

double slaEqeqx ( double date );

void slaEqgal ( double dr, double dd, double *dl, double *db );

void slaEtrms ( double ep, double ev[3] );

void slaEuler ( char *order, float phi, float theta, float psi,
                float rmat[3][3] );

void slaEvp ( double date, double deqx,
              double dvb[3], double dpb[3],
              double dvh[3], double dph[3] );

void slaFitxy ( int itype, int np, double xye[][2], double xym[][2],
                double coeffs[6], int *j );

void slaFk425 ( double r1950, double d1950, double dr1950,
                double dd1950, double p1950, double v1950,
                double *r2000, double *d2000, double *dr2000,
                double *dd2000, double *p2000, double *v2000 );

void slaFk45z ( double r1950, double d1950, double bepoch,
                double *r2000, double *d2000 );

void slaFk524 ( double r2000, double d2000, double dr2000,
                double dd2000, double p2000, double v2000,
                double *r1950, double *d1950, double *dr1950,
                double *dd1950, double *p1950, double *v1950 );

void slaFk52h ( double r5, double d5, double dr5, double dd5,
                double *dr, double *dh, double *drh, double *ddh );

void slaFk54z ( double r2000, double d2000, double bepoch,
                double *r1950, double *d1950,
                double *dr1950, double *dd1950 );

void slaFk5hz ( double r5, double d5, double epoch,
                double *rh, double *dh );

void slaFlotin ( char *string, int *nstrt, float *reslt, int *jflag );

void slaGaleq ( double dl, double db, double *dr, double *dd );

void slaGalsup ( double dl, double db, double *dsl, double *dsb );

void slaGe50 ( double dl, double db, double *dr, double *dd );

void slaGeoc ( double p, double h, double *r, double *z );

double slaGmst ( double ut1 );

double slaGmsta ( double date, double ut1 );

void slaH2e ( float az, float el, float phi, float *ha, float *dec );

void slaH2fk5 ( double dr, double dh, double drh, double ddh,
                double *r5, double *d5, double *dr5, double *dd5 );

void slaHfk5z ( double rh, double dh, double epoch,
                double *r5, double *d5, double *dr5, double *dd5 );

void slaImxv ( float rm[3][3], float va[3], float vb[3] );

void slaInt2in ( char *string, int *nstrt, int *ireslt, int *jflag );

void slaIntin ( char *string, int *nstrt, long *ireslt, int *jflag );

void slaInvf ( double fwds[6], double bkwds[6], int *j );

void slaKbj ( int jb, double e, char *k, int *j );

void slaM2av ( float rmat[3][3], float axvec[3] );

void slaMap ( double rm, double dm, double pr, double pd,
              double px, double rv, double eq, double date,
              double *ra, double *da );

void slaMappa ( double eq, double date, double amprms[21] );

void slaMapqk ( double rm, double dm, double pr, double pd,
                double px, double rv, double amprms[21],
                double *ra, double *da );

void slaMapqkz ( double rm, double dm, double amprms[21],
                 double *ra, double *da );

void slaMoon ( int iy, int id, float fd, float posvel[6] );

void slaMxm ( float a[3][3], float b[3][3], float c[3][3] );

void slaMxv ( float rm[3][3], float va[3], float vb[3] );

void slaNut ( double date, double rmatn[3][3] );

void slaNutc ( double date, double *dpsi, double *deps, double *eps0 );

void slaNutc80 ( double date, double *dpsi, double *deps, double *eps0 );

void slaOap ( char *type, double ob1, double ob2, double date,
              double dut, double elongm, double phim, double hm,
              double xp, double yp, double tdk, double pmb,
              double rh, double wl, double tlr,
              double *rap, double *dap );

void slaOapqk ( char *type, double ob1, double ob2, double aoprms[14],
                double *rap, double *dap );

void slaObs ( int n, char *c, char *name, double *w, double *p, double *h );

double slaPa ( double ha, double dec, double phi );

double slaPav ( float v1[3], float v2[3] );

void slaPcd ( double disco, double *x, double *y );

void slaPda2h ( double p, double d, double a,
                double *h1, int *j1, double *h2, int *j2 );

void slaPdq2h ( double p, double d, double q,
                double *h1, int *j1, double *h2, int *j2 );

void slaPermut ( int n, int istate[], int iorder[], int *j );

void slaPertel (int jform, double date0, double date1,
                double epoch0, double orbi0, double anode0,
                double perih0, double aorq0, double e0, double am0,
                double *epoch1, double *orbi1, double *anode1,
                double *perih1, double *aorq1, double *e1, double *am1,
                int *jstat );

void slaPertue ( double date, double u[], int *jstat );

void slaPlanel ( double date, int jform, double epoch, double orbinc,
                 double anode, double perih, double aorq,  double e,
                 double aorl, double dm, double pv[6], int *jstat );

void slaPlanet ( double date, int np, double pv[6], int *j );

void slaPlante ( double date, double elong, double phi, int jform,
                 double epoch, double orbinc, double anode, double perih,
                 double aorq, double e, double aorl, double dm,
                 double *ra, double *dec, double *r, int *jstat );

void slaPlantu ( double date, double elong, double phi, double u[],
                 double *ra, double *dec, double *r, int *jstat );

void slaPm ( double r0, double d0, double pr, double pd,
             double px, double rv, double ep0, double ep1,
             double *r1, double *d1 );

void slaPolmo ( double elongm, double phim, double xp, double yp,
                double *elong, double *phi, double *daz );

void slaPrebn ( double bep0, double bep1, double rmatp[3][3] );

void slaPrec ( double ep0, double ep1, double rmatp[3][3] );

void slaPrecl ( double ep0, double ep1, double rmatp[3][3] );

void slaPreces ( char sys[3], double ep0, double ep1,
                 double *ra, double *dc );

void slaPrenut ( double epoch, double date, double rmatpn[3][3] );

void slaPv2el ( double pv[], double date, double pmass, int jformr,
                int *jform, double *epoch, double *orbinc,
                double *anode, double *perih, double *aorq, double *e,
                double *aorl, double *dm, int *jstat );

void slaPv2ue ( double pv[], double date, double pmass,
                double u[], int *jstat );

void slaPvobs ( double p, double h, double stl, double pv[6] );

void slaPxy ( int np, double xye[][2], double xym[][2],
              double coeffs[6],
              double xyp[][2], double *xrms, double *yrms, double *rrms );

float slaRange ( float angle );

float slaRanorm ( float angle );

double slaRcc ( double tdb, double ut1, double wl, double u, double v );

void slaRdplan ( double date, int np, double elong, double phi,
                 double *ra, double *dec, double *diam );

void slaRefco ( double hm, double tdk, double pmb, double rh,
                double wl, double phi, double tlr, double eps,
                double *refa, double *refb );

void slaRefcoq ( double tdk, double pmb, double rh, double wl,
                double *refa, double *refb );

void slaRefro ( double zobs, double hm, double tdk, double pmb,
                double rh, double wl, double phi, double tlr, double eps,
                double *ref );

void slaRefv ( double vu[3], double refa, double refb, double vr[3] );

void slaRefz ( double zu, double refa, double refb, double *zr );

float slaRverot ( float phi, float ra, float da, float st );

float slaRvgalc ( float r2000, float d2000 );

float slaRvlg ( float r2000, float d2000 );

float slaRvlsrd ( float r2000, float d2000 );

float slaRvlsrk ( float r2000, float d2000 );

void slaS2tp ( float ra, float dec, float raz, float decz,
               float *xi, float *eta, int *j );

float slaSep ( float a1, float b1, float a2, float b2 );

float slaSepv ( float v1[3], float v2[3] );

void slaSmat ( int n, float *a, float *y, float *d, int *jf, int *iw );

void slaSubet ( double rc, double dc, double eq,
                double *rm, double *dm );

void slaSupgal ( double dsl, double dsb, double *dl, double *db );

void slaSvd ( int m, int n, int mp, int np,
              double *a, double *w, double *v, double *work,
              int *jstat );

void slaSvdcov ( int n, int np, int nc,
                 double *w, double *v, double *work, double *cvm );

void slaSvdsol ( int m, int n, int mp, int np,
                 double *b, double *u, double *w, double *v,
                 double *work, double *x );

void slaTp2s ( float xi, float eta, float raz, float decz,
               float *ra, float *dec );

void slaTp2v ( float xi, float eta, float v0[3], float v[3] );

void slaTps2c ( float xi, float eta, float ra, float dec,
                float *raz1, float *decz1,
                float *raz2, float *decz2, int *n );

void slaTpv2c ( float xi, float eta, float v[3],
                float v01[3], float v02[3], int *n );

void slaUe2el ( double u[], int jformr,
                int *jform, double *epoch, double *orbinc,
                double *anode, double *perih, double *aorq, double *e,
                double *aorl, double *dm, int *jstat );

void slaUe2pv ( double date, double u[], double pv[], int *jstat );

void slaUnpcd ( double disco, double *x, double *y );

void slaV2tp ( float v[3], float v0[3], float *xi, float *eta, int *j );

float slaVdv ( float va[3], float vb[3] );

void slaVn ( float v[3], float uv[3], float *vm );

void slaVxv ( float va[3], float vb[3], float vc[3] );

void slaXy2xy ( double x1, double y1, double coeffs[6],
                double *x2, double *y2 );

double slaZd ( double ha, double dec, double phi );

#ifdef __cplusplus
}
#endif

#endif
