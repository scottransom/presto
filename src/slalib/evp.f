      SUBROUTINE sla_EVP (DATE, DEQX, DVB, DPB, DVH, DPH)
*+
*     - - - -
*      E V P
*     - - - -
*
*  Barycentric and heliocentric velocity and position of the Earth
*
*  All arguments are double precision
*
*  Given:
*
*     DATE          TDB (loosely ET) as a Modified Julian Date
*                                         (JD-2400000.5)
*
*     DEQX          Julian Epoch (e.g. 2000.0D0) of mean equator and
*                   equinox of the vectors returned.  If DEQX .LE. 0D0,
*                   all vectors are referred to the mean equator and
*                   equinox (FK5) of epoch DATE.
*
*  Returned (all 3D Cartesian vectors):
*
*     DVB,DPB       barycentric velocity, position (AU/s, AU)
*     DVH,DPH       heliocentric velocity, position (AU/s, AU)
*
*  Called:  sla_EPJ, sla_PREC
*
*  Notes:
*
*  1  This routine is accurate enough for many purposes but faster and
*     more compact than the sla_EPV routine.  The maximum deviations
*     from the JPL DE96 ephemeris are as follows:
*
*       barycentric velocity         0.42  m/s
*       barycentric position         6900  km
*
*       heliocentric velocity        0.42  m/s
*       heliocentric position        1600  km
*
*  2  The routine is adapted from the BARVEL and BARCOR subroutines of
*     Stumpff (1980).  Most of the changes are merely cosmetic and do
*     not affect the results at all.  However, some adjustments have
*     been made so as to give results that refer to the IAU 1976 'FK5'
*     equinox and precession, although the differences these changes
*     make relative to the results from Stumpff's original 'FK4' version
*     are smaller than the inherent accuracy of the algorithm.  One
*     minor shortcoming in the original routines that has NOT been
*     corrected is that better numerical accuracy could be achieved if
*     the various polynomial evaluations were nested.
*
*  Reference:
*
*    Stumpff, P., Astron.Astrophys.Suppl.Ser. 41, 1-8 (1980).
*
*  Last revision:   7 April 2005
*
*  Copyright P.T.Wallace.  All rights reserved.
*
*  License:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA
*
*-

      IMPLICIT NONE

      DOUBLE PRECISION DATE,DEQX,DVB(3),DPB(3),DVH(3),DPH(3)

      INTEGER IDEQ,I,J,K

      REAL CC2PI,CCSEC3,CCSGD,CCKM,CCMLD,CCFDI,CCIM,T,TSQ,A,PERTL,
     :     PERTLD,PERTR,PERTRD,COSA,SINA,ESQ,E,PARAM,TWOE,TWOG,G,
     :     PHI,F,SINF,COSF,PHID,PSID,PERTP,PERTPD,TL,SINLM,COSLM,
     :     SIGMA,B,PLON,POMG,PECC,FLATM,FLAT

      DOUBLE PRECISION DC2PI,DS2R,DCSLD,DC1MME,DT,DTSQ,DLOCAL,DML,
     :                 DEPS,DPARAM,DPSI,D1PDRO,DRD,DRLD,DTL,DSINLS,
     :                 DCOSLS,DXHD,DYHD,DZHD,DXBD,DYBD,DZBD,DCOSEP,
     :                 DSINEP,DYAHD,DZAHD,DYABD,DZABD,DR,
     :                 DXH,DYH,DZH,DXB,DYB,DZB,DYAH,DZAH,DYAB,
     :                 DZAB,DEPJ,DEQCOR,B1950

      REAL SN(4),CCSEL(3,17),CCAMPS(5,15),CCSEC(3,4),CCAMPM(4,3),
     :     CCPAMV(4),CCPAM(4),FORBEL(7),SORBEL(17),SINLP(4),COSLP(4)
      EQUIVALENCE (SORBEL(1),E),(FORBEL(1),G)

      DOUBLE PRECISION DCFEL(3,8),DCEPS(3),DCARGS(2,15),DCARGM(2,3),
     :                 DPREMA(3,3),W,VW(3)

      DOUBLE PRECISION sla_EPJ

      PARAMETER (DC2PI=6.2831853071796D0,CC2PI=6.283185)
      PARAMETER (DS2R=0.7272205216643D-4)
      PARAMETER (B1950=1949.9997904423D0)

*
*   Constants DCFEL(I,K) of fast changing elements
*                     I=1                I=2              I=3
      DATA DCFEL/ 1.7400353D+00, 6.2833195099091D+02, 5.2796D-06,
     :            6.2565836D+00, 6.2830194572674D+02,-2.6180D-06,
     :            4.7199666D+00, 8.3997091449254D+03,-1.9780D-05,
     :            1.9636505D-01, 8.4334662911720D+03,-5.6044D-05,
     :            4.1547339D+00, 5.2993466764997D+01, 5.8845D-06,
     :            4.6524223D+00, 2.1354275911213D+01, 5.6797D-06,
     :            4.2620486D+00, 7.5025342197656D+00, 5.5317D-06,
     :            1.4740694D+00, 3.8377331909193D+00, 5.6093D-06/

*
*   Constants DCEPS and CCSEL(I,K) of slowly changing elements
*                      I=1           I=2           I=3
      DATA DCEPS/  4.093198D-01,-2.271110D-04,-2.860401D-08 /
      DATA CCSEL/  1.675104E-02,-4.179579E-05,-1.260516E-07,
     :             2.220221E-01, 2.809917E-02, 1.852532E-05,
     :             1.589963E+00, 3.418075E-02, 1.430200E-05,
     :             2.994089E+00, 2.590824E-02, 4.155840E-06,
     :             8.155457E-01, 2.486352E-02, 6.836840E-06,
     :             1.735614E+00, 1.763719E-02, 6.370440E-06,
     :             1.968564E+00, 1.524020E-02,-2.517152E-06,
     :             1.282417E+00, 8.703393E-03, 2.289292E-05,
     :             2.280820E+00, 1.918010E-02, 4.484520E-06,
     :             4.833473E-02, 1.641773E-04,-4.654200E-07,
     :             5.589232E-02,-3.455092E-04,-7.388560E-07,
     :             4.634443E-02,-2.658234E-05, 7.757000E-08,
     :             8.997041E-03, 6.329728E-06,-1.939256E-09,
     :             2.284178E-02,-9.941590E-05, 6.787400E-08,
     :             4.350267E-02,-6.839749E-05,-2.714956E-07,
     :             1.348204E-02, 1.091504E-05, 6.903760E-07,
     :             3.106570E-02,-1.665665E-04,-1.590188E-07/

*
*   Constants of the arguments of the short-period perturbations
*   by the planets:   DCARGS(I,K)
*                       I=1               I=2
      DATA DCARGS/ 5.0974222D+00,-7.8604195454652D+02,
     :             3.9584962D+00,-5.7533848094674D+02,
     :             1.6338070D+00,-1.1506769618935D+03,
     :             2.5487111D+00,-3.9302097727326D+02,
     :             4.9255514D+00,-5.8849265665348D+02,
     :             1.3363463D+00,-5.5076098609303D+02,
     :             1.6072053D+00,-5.2237501616674D+02,
     :             1.3629480D+00,-1.1790629318198D+03,
     :             5.5657014D+00,-1.0977134971135D+03,
     :             5.0708205D+00,-1.5774000881978D+02,
     :             3.9318944D+00, 5.2963464780000D+01,
     :             4.8989497D+00, 3.9809289073258D+01,
     :             1.3097446D+00, 7.7540959633708D+01,
     :             3.5147141D+00, 7.9618578146517D+01,
     :             3.5413158D+00,-5.4868336758022D+02/

*
*   Amplitudes CCAMPS(N,K) of the short-period perturbations
*           N=1          N=2          N=3          N=4          N=5
      DATA CCAMPS/
     : -2.279594E-5, 1.407414E-5, 8.273188E-6, 1.340565E-5,-2.490817E-7,
     : -3.494537E-5, 2.860401E-7, 1.289448E-7, 1.627237E-5,-1.823138E-7,
     :  6.593466E-7, 1.322572E-5, 9.258695E-6,-4.674248E-7,-3.646275E-7,
     :  1.140767E-5,-2.049792E-5,-4.747930E-6,-2.638763E-6,-1.245408E-7,
     :  9.516893E-6,-2.748894E-6,-1.319381E-6,-4.549908E-6,-1.864821E-7,
     :  7.310990E-6,-1.924710E-6,-8.772849E-7,-3.334143E-6,-1.745256E-7,
     : -2.603449E-6, 7.359472E-6, 3.168357E-6, 1.119056E-6,-1.655307E-7,
     : -3.228859E-6, 1.308997E-7, 1.013137E-7, 2.403899E-6,-3.736225E-7,
     :  3.442177E-7, 2.671323E-6, 1.832858E-6,-2.394688E-7,-3.478444E-7,
     :  8.702406E-6,-8.421214E-6,-1.372341E-6,-1.455234E-6,-4.998479E-8,
     : -1.488378E-6,-1.251789E-5, 5.226868E-7,-2.049301E-7, 0.0E0,
     : -8.043059E-6,-2.991300E-6, 1.473654E-7,-3.154542E-7, 0.0E0,
     :  3.699128E-6,-3.316126E-6, 2.901257E-7, 3.407826E-7, 0.0E0,
     :  2.550120E-6,-1.241123E-6, 9.901116E-8, 2.210482E-7, 0.0E0,
     : -6.351059E-7, 2.341650E-6, 1.061492E-6, 2.878231E-7, 0.0E0/

*
*   Constants of the secular perturbations in longitude
*   CCSEC3 and CCSEC(N,K)
*                      N=1           N=2           N=3
      DATA CCSEC3/-7.757020E-08/,
     :     CCSEC/  1.289600E-06, 5.550147E-01, 2.076942E+00,
     :             3.102810E-05, 4.035027E+00, 3.525565E-01,
     :             9.124190E-06, 9.990265E-01, 2.622706E+00,
     :             9.793240E-07, 5.508259E+00, 1.559103E+01/

*   Sidereal rate DCSLD in longitude, rate CCSGD in mean anomaly
      DATA DCSLD/1.990987D-07/,
     :     CCSGD/1.990969E-07/

*   Some constants used in the calculation of the lunar contribution
      DATA CCKM/3.122140E-05/,
     :     CCMLD/2.661699E-06/,
     :     CCFDI/2.399485E-07/

*
*   Constants DCARGM(I,K) of the arguments of the perturbations
*   of the motion of the Moon
*                       I=1               I=2
      DATA DCARGM/  5.1679830D+00, 8.3286911095275D+03,
     :              5.4913150D+00,-7.2140632838100D+03,
     :              5.9598530D+00, 1.5542754389685D+04/

*
*   Amplitudes CCAMPM(N,K) of the perturbations of the Moon
*            N=1          N=2           N=3           N=4
      DATA CCAMPM/
     :  1.097594E-01, 2.896773E-07, 5.450474E-02, 1.438491E-07,
     : -2.223581E-02, 5.083103E-08, 1.002548E-02,-2.291823E-08,
     :  1.148966E-02, 5.658888E-08, 8.249439E-03, 4.063015E-08/

*
*   CCPAMV(K)=A*M*DL/DT (planets), DC1MME=1-MASS(Earth+Moon)
      DATA CCPAMV/8.326827E-11,1.843484E-11,1.988712E-12,1.881276E-12/
      DATA DC1MME/0.99999696D0/

*   CCPAM(K)=A*M(planets), CCIM=INCLINATION(Moon)
      DATA CCPAM/4.960906E-3,2.727436E-3,8.392311E-4,1.556861E-3/
      DATA CCIM/8.978749E-2/




*
*   EXECUTION
*   ---------

*   Control parameter IDEQ, and time arguments
      IDEQ = 0
      IF (DEQX.GT.0D0) IDEQ=1
      DT = (DATE-15019.5D0)/36525D0
      T = REAL(DT)
      DTSQ = DT*DT
      TSQ = REAL(DTSQ)

*   Values of all elements for the instant DATE
      DO K=1,8
         DLOCAL = MOD(DCFEL(1,K)+DT*DCFEL(2,K)+DTSQ*DCFEL(3,K), DC2PI)
         IF (K.EQ.1) THEN
            DML = DLOCAL
         ELSE
            FORBEL(K-1) = REAL(DLOCAL)
         END IF
      END DO
      DEPS = MOD(DCEPS(1)+DT*DCEPS(2)+DTSQ*DCEPS(3), DC2PI)
      DO K=1,17
         SORBEL(K) = MOD(CCSEL(1,K)+T*CCSEL(2,K)+TSQ*CCSEL(3,K),
     :                   CC2PI)
      END DO

*   Secular perturbations in longitude
      DO K=1,4
         A = MOD(CCSEC(2,K)+T*CCSEC(3,K), CC2PI)
         SN(K) = SIN(A)
      END DO

*   Periodic perturbations of the EMB (Earth-Moon barycentre)
      PERTL =  CCSEC(1,1)          *SN(1) +CCSEC(1,2)*SN(2)+
     :        (CCSEC(1,3)+T*CCSEC3)*SN(3) +CCSEC(1,4)*SN(4)
      PERTLD = 0.0
      PERTR = 0.0
      PERTRD = 0.0
      DO K=1,15
         A = SNGL(MOD(DCARGS(1,K)+DT*DCARGS(2,K), DC2PI))
         COSA = COS(A)
         SINA = SIN(A)
         PERTL = PERTL + CCAMPS(1,K)*COSA+CCAMPS(2,K)*SINA
         PERTR = PERTR + CCAMPS(3,K)*COSA+CCAMPS(4,K)*SINA
         IF (K.LT.11) THEN
            PERTLD = PERTLD+
     :               (CCAMPS(2,K)*COSA-CCAMPS(1,K)*SINA)*CCAMPS(5,K)
            PERTRD = PERTRD+
     :               (CCAMPS(4,K)*COSA-CCAMPS(3,K)*SINA)*CCAMPS(5,K)
         END IF
      END DO

*   Elliptic part of the motion of the EMB
      ESQ = E*E
      DPARAM = 1D0-DBLE(ESQ)
      PARAM = REAL(DPARAM)
      TWOE = E+E
      TWOG = G+G
      PHI = TWOE*((1.0-ESQ*0.125)*SIN(G)+E*0.625*SIN(TWOG)
     :          +ESQ*0.54166667*SIN(G+TWOG) )
      F = G+PHI
      SINF = SIN(F)
      COSF = COS(F)
      DPSI = DPARAM/(1D0+DBLE(E*COSF))
      PHID = TWOE*CCSGD*((1.0+ESQ*1.5)*COSF+E*(1.25-SINF*SINF*0.5))
      PSID = CCSGD*E*SINF/SQRT(PARAM)

*   Perturbed heliocentric motion of the EMB
      D1PDRO = 1D0+DBLE(PERTR)
      DRD = D1PDRO*(DBLE(PSID)+DPSI*DBLE(PERTRD))
      DRLD = D1PDRO*DPSI*(DCSLD+DBLE(PHID)+DBLE(PERTLD))
      DTL = MOD(DML+DBLE(PHI)+DBLE(PERTL), DC2PI)
      DSINLS = SIN(DTL)
      DCOSLS = COS(DTL)
      DXHD = DRD*DCOSLS-DRLD*DSINLS
      DYHD = DRD*DSINLS+DRLD*DCOSLS

*   Influence of eccentricity, evection and variation on the
*   geocentric motion of the Moon
      PERTL = 0.0
      PERTLD = 0.0
      PERTP = 0.0
      PERTPD = 0.0
      DO K=1,3
         A = SNGL(MOD(DCARGM(1,K)+DT*DCARGM(2,K), DC2PI))
         SINA = SIN(A)
         COSA = COS(A)
         PERTL = PERTL +CCAMPM(1,K)*SINA
         PERTLD = PERTLD+CCAMPM(2,K)*COSA
         PERTP = PERTP +CCAMPM(3,K)*COSA
         PERTPD = PERTPD-CCAMPM(4,K)*SINA
      END DO

*   Heliocentric motion of the Earth
      TL = FORBEL(2)+PERTL
      SINLM = SIN(TL)
      COSLM = COS(TL)
      SIGMA = CCKM/(1.0+PERTP)
      A = SIGMA*(CCMLD+PERTLD)
      B = SIGMA*PERTPD
      DXHD = DXHD+DBLE(A*SINLM)+DBLE(B*COSLM)
      DYHD = DYHD-DBLE(A*COSLM)+DBLE(B*SINLM)
      DZHD =     -DBLE(SIGMA*CCFDI*COS(FORBEL(3)))

*   Barycentric motion of the Earth
      DXBD = DXHD*DC1MME
      DYBD = DYHD*DC1MME
      DZBD = DZHD*DC1MME
      DO K=1,4
         PLON = FORBEL(K+3)
         POMG = SORBEL(K+1)
         PECC = SORBEL(K+9)
         TL = MOD(PLON+2.0*PECC*SIN(PLON-POMG), CC2PI)
         SINLP(K) = SIN(TL)
         COSLP(K) = COS(TL)
         DXBD = DXBD+DBLE(CCPAMV(K)*(SINLP(K)+PECC*SIN(POMG)))
         DYBD = DYBD-DBLE(CCPAMV(K)*(COSLP(K)+PECC*COS(POMG)))
         DZBD = DZBD-DBLE(CCPAMV(K)*SORBEL(K+13)*COS(PLON-SORBEL(K+5)))
      END DO

*   Transition to mean equator of date
      DCOSEP = COS(DEPS)
      DSINEP = SIN(DEPS)
      DYAHD = DCOSEP*DYHD-DSINEP*DZHD
      DZAHD = DSINEP*DYHD+DCOSEP*DZHD
      DYABD = DCOSEP*DYBD-DSINEP*DZBD
      DZABD = DSINEP*DYBD+DCOSEP*DZBD

*   Heliocentric coordinates of the Earth
      DR = DPSI*D1PDRO
      FLATM = CCIM*SIN(FORBEL(3))
      A = SIGMA*COS(FLATM)
      DXH = DR*DCOSLS-DBLE(A*COSLM)
      DYH = DR*DSINLS-DBLE(A*SINLM)
      DZH =          -DBLE(SIGMA*SIN(FLATM))

*   Barycentric coordinates of the Earth
      DXB = DXH*DC1MME
      DYB = DYH*DC1MME
      DZB = DZH*DC1MME
      DO K=1,4
         FLAT = SORBEL(K+13)*SIN(FORBEL(K+3)-SORBEL(K+5))
         A = CCPAM(K)*(1.0-SORBEL(K+9)*COS(FORBEL(K+3)-SORBEL(K+1)))
         B = A*COS(FLAT)
         DXB = DXB-DBLE(B*COSLP(K))
         DYB = DYB-DBLE(B*SINLP(K))
         DZB = DZB-DBLE(A*SIN(FLAT))
      END DO

*   Transition to mean equator of date
      DYAH = DCOSEP*DYH-DSINEP*DZH
      DZAH = DSINEP*DYH+DCOSEP*DZH
      DYAB = DCOSEP*DYB-DSINEP*DZB
      DZAB = DSINEP*DYB+DCOSEP*DZB

*   Copy result components into vectors, correcting for FK4 equinox
      DEPJ=sla_EPJ(DATE)
      DEQCOR = DS2R*(0.035D0+0.00085D0*(DEPJ-B1950))
      DVH(1) = DXHD-DEQCOR*DYAHD
      DVH(2) = DYAHD+DEQCOR*DXHD
      DVH(3) = DZAHD
      DVB(1) = DXBD-DEQCOR*DYABD
      DVB(2) = DYABD+DEQCOR*DXBD
      DVB(3) = DZABD
      DPH(1) = DXH-DEQCOR*DYAH
      DPH(2) = DYAH+DEQCOR*DXH
      DPH(3) = DZAH
      DPB(1) = DXB-DEQCOR*DYAB
      DPB(2) = DYAB+DEQCOR*DXB
      DPB(3) = DZAB

*   Was precession to another equinox requested?
      IF (IDEQ.NE.0) THEN

*     Yes: compute precession matrix from MJD DATE to Julian epoch DEQX
         CALL sla_PREC(DEPJ,DEQX,DPREMA)

*     Rotate DVH
         DO J=1,3
            W=0D0
            DO I=1,3
               W=W+DPREMA(J,I)*DVH(I)
            END DO
            VW(J)=W
         END DO
         DO J=1,3
            DVH(J)=VW(J)
         END DO

*     Rotate DVB
         DO J=1,3
            W=0D0
            DO I=1,3
               W=W+DPREMA(J,I)*DVB(I)
            END DO
            VW(J)=W
         END DO
         DO J=1,3
            DVB(J)=VW(J)
         END DO

*     Rotate DPH
         DO J=1,3
            W=0D0
            DO I=1,3
               W=W+DPREMA(J,I)*DPH(I)
            END DO
            VW(J)=W
         END DO
         DO J=1,3
            DPH(J)=VW(J)
         END DO

*     Rotate DPB
         DO J=1,3
            W=0D0
            DO I=1,3
               W=W+DPREMA(J,I)*DPB(I)
            END DO
            VW(J)=W
         END DO
         DO J=1,3
            DPB(J)=VW(J)
         END DO
      END IF

      END
