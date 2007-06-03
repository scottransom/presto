      SUBROUTINE sla_MOON (IY, ID, FD, PV)
*+
*     - - - - -
*      M O O N
*     - - - - -
*
*  Approximate geocentric position and velocity of the Moon
*  (single precision).
*
*  Given:
*     IY       i       year
*     ID       i       day in year (1 = Jan 1st)
*     FD       r       fraction of day
*
*  Returned:
*     PV       r(6)    Moon position & velocity vector
*
*  Notes:
*
*  1  The date and time is TDB (loosely ET) in a Julian calendar
*     which has been aligned to the ordinary Gregorian
*     calendar for the interval 1900 March 1 to 2100 February 28.
*     The year and day can be obtained by calling sla_CALYD or
*     sla_CLYD.
*
*  2  The Moon 6-vector is Moon centre relative to Earth centre,
*     mean equator and equinox of date.  Position part, PV(1-3),
*     is in AU;  velocity part, PV(4-6), is in AU/sec.
*
*  3  The position is accurate to better than 0.5 arcminute
*     in direction and 1000 km in distance.  The velocity
*     is accurate to better than 0.5"/hour in direction and
*     4 m/s in distance.  (RMS figures with respect to JPL DE200
*     for the interval 1960-2025 are 14 arcsec and 0.2 arcsec/hour in
*     longitude, 9 arcsec and 0.2 arcsec/hour in latitude, 350 km and
*     2 m/s in distance.)  Note that the distance accuracy is
*     comparatively poor because this routine is principally intended
*     for computing topocentric direction.
*
*  4  This routine is only a partial implementation of the original
*     Meeus algorithm (reference below), which offers 4 times the
*     accuracy in direction and 30 times the accuracy in distance
*     when fully implemented (as it is in sla_DMOON).
*
*  Reference:
*     Meeus, l'Astronomie, June 1984, p348.
*
*  Called:  sla_CS2C6
*
*  P.T.Wallace   Starlink   8 December 1994
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
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

      INTEGER IY,ID
      REAL FD,PV(6)

      INTEGER ITP(4,4),ITL(4,39),ITB(4,29),I,IY4,N
      REAL D2R,RATCON,ERADAU
      REAL ELP0,ELP1,ELP1I,ELP1F
      REAL EM0,EM1,EM1F
      REAL EMP0,EMP1,EMP1I,EMP1F
      REAL D0,D1,D1I,D1F
      REAL F0,F1,F1I,F1F
      REAL TL(39)
      REAL TB(29)
      REAL TP(4)
      REAL YI,YF,T,ELP,EM,EMP,D,F,EL,ELD,COEFF,CEM,CEMP
      REAL CD,CF,THETA,THETAD,B,BD,P,PD,SP,R,RD
      REAL V(6),EPS,SINEPS,COSEPS

*  Degrees to radians
      PARAMETER (D2R=1.745329252E-2)

*  Rate conversion factor:  D2R**2/(86400*365.25)
      PARAMETER (RATCON=9.652743551E-12)

*  Earth radius in AU:  6378.137/149597870
      PARAMETER (ERADAU=4.2635212653763E-5)

*
*  Coefficients for fundamental arguments
*
*  Fixed term (deg), term in T (deg & whole revs + fraction per year)
*
*  Moon's mean longitude
      DATA ELP0,ELP1,ELP1I,ELP1F /
     :            270.434164, 4812.678831, 4680., 132.678831 /
*
*  Sun's mean anomaly
      DATA EM0,EM1,EM1F /
     :            358.475833,  359.990498,        359.990498 /
*
*  Moon's mean anomaly
      DATA EMP0,EMP1,EMP1I,EMP1F /
     :            296.104608, 4771.988491, 4680.,  91.988491 /
*
*  Moon's mean elongation
      DATA D0,D1,D1I,D1F /
     :            350.737486,  4452.671142, 4320., 132.671142 /
*
*  Mean distance of the Moon from its ascending node
      DATA F0,F1,F1I,F1F /
     :             11.250889, 4832.020251, 4680., 152.020251 /

*
*  Coefficients for Moon position
*
*   T(N)       = coefficient of term (deg)
*   IT(N,1-4)  = coefficients of M, M', D, F in argument
*
*  Longitude
*                                         M   M'  D   F
      DATA TL( 1)/            +6.288750                     /,
     :     (ITL(I, 1),I=1,4)/             0, +1,  0,  0     /
      DATA TL( 2)/            +1.274018                     /,
     :     (ITL(I, 2),I=1,4)/             0, -1, +2,  0     /
      DATA TL( 3)/            +0.658309                     /,
     :     (ITL(I, 3),I=1,4)/             0,  0, +2,  0     /
      DATA TL( 4)/            +0.213616                     /,
     :     (ITL(I, 4),I=1,4)/             0, +2,  0,  0     /
      DATA TL( 5)/            -0.185596                     /,
     :     (ITL(I, 5),I=1,4)/            +1,  0,  0,  0     /
      DATA TL( 6)/            -0.114336                     /,
     :     (ITL(I, 6),I=1,4)/             0,  0,  0, +2     /
      DATA TL( 7)/            +0.058793                     /,
     :     (ITL(I, 7),I=1,4)/             0, -2, +2,  0     /
      DATA TL( 8)/            +0.057212                     /,
     :     (ITL(I, 8),I=1,4)/            -1, -1, +2,  0     /
      DATA TL( 9)/            +0.053320                     /,
     :     (ITL(I, 9),I=1,4)/             0, +1, +2,  0     /
      DATA TL(10)/            +0.045874                     /,
     :     (ITL(I,10),I=1,4)/            -1,  0, +2,  0     /
      DATA TL(11)/            +0.041024                     /,
     :     (ITL(I,11),I=1,4)/            -1, +1,  0,  0     /
      DATA TL(12)/            -0.034718                     /,
     :     (ITL(I,12),I=1,4)/             0,  0, +1,  0     /
      DATA TL(13)/            -0.030465                     /,
     :     (ITL(I,13),I=1,4)/            +1, +1,  0,  0     /
      DATA TL(14)/            +0.015326                     /,
     :     (ITL(I,14),I=1,4)/             0,  0, +2, -2     /
      DATA TL(15)/            -0.012528                     /,
     :     (ITL(I,15),I=1,4)/             0, +1,  0, +2     /
      DATA TL(16)/            -0.010980                     /,
     :     (ITL(I,16),I=1,4)/             0, -1,  0, +2     /
      DATA TL(17)/            +0.010674                     /,
     :     (ITL(I,17),I=1,4)/             0, -1, +4,  0     /
      DATA TL(18)/            +0.010034                     /,
     :     (ITL(I,18),I=1,4)/             0, +3,  0,  0     /
      DATA TL(19)/            +0.008548                     /,
     :     (ITL(I,19),I=1,4)/             0, -2, +4,  0     /
      DATA TL(20)/            -0.007910                     /,
     :     (ITL(I,20),I=1,4)/            +1, -1, +2,  0     /
      DATA TL(21)/            -0.006783                     /,
     :     (ITL(I,21),I=1,4)/            +1,  0, +2,  0     /
      DATA TL(22)/            +0.005162                     /,
     :     (ITL(I,22),I=1,4)/             0, +1, -1,  0     /
      DATA TL(23)/            +0.005000                     /,
     :     (ITL(I,23),I=1,4)/            +1,  0, +1,  0     /
      DATA TL(24)/            +0.004049                     /,
     :     (ITL(I,24),I=1,4)/            -1, +1, +2,  0     /
      DATA TL(25)/            +0.003996                     /,
     :     (ITL(I,25),I=1,4)/             0, +2, +2,  0     /
      DATA TL(26)/            +0.003862                     /,
     :     (ITL(I,26),I=1,4)/             0,  0, +4,  0     /
      DATA TL(27)/            +0.003665                     /,
     :     (ITL(I,27),I=1,4)/             0, -3, +2,  0     /
      DATA TL(28)/            +0.002695                     /,
     :     (ITL(I,28),I=1,4)/            -1, +2,  0,  0     /
      DATA TL(29)/            +0.002602                     /,
     :     (ITL(I,29),I=1,4)/             0, +1, -2, -2     /
      DATA TL(30)/            +0.002396                     /,
     :     (ITL(I,30),I=1,4)/            -1, -2, +2,  0     /
      DATA TL(31)/            -0.002349                     /,
     :     (ITL(I,31),I=1,4)/             0, +1, +1,  0     /
      DATA TL(32)/            +0.002249                     /,
     :     (ITL(I,32),I=1,4)/            -2,  0, +2,  0     /
      DATA TL(33)/            -0.002125                     /,
     :     (ITL(I,33),I=1,4)/            +1, +2,  0,  0     /
      DATA TL(34)/            -0.002079                     /,
     :     (ITL(I,34),I=1,4)/            +2,  0,  0,  0     /
      DATA TL(35)/            +0.002059                     /,
     :     (ITL(I,35),I=1,4)/            -2, -1, +2,  0     /
      DATA TL(36)/            -0.001773                     /,
     :     (ITL(I,36),I=1,4)/             0, +1, +2, -2     /
      DATA TL(37)/            -0.001595                     /,
     :     (ITL(I,37),I=1,4)/             0,  0, +2, +2     /
      DATA TL(38)/            +0.001220                     /,
     :     (ITL(I,38),I=1,4)/            -1, -1, +4,  0     /
      DATA TL(39)/            -0.001110                     /,
     :     (ITL(I,39),I=1,4)/             0, +2,  0, +2     /
*
*  Latitude
*                                         M   M'  D   F
      DATA TB( 1)/            +5.128189                     /,
     :     (ITB(I, 1),I=1,4)/             0,  0,  0, +1     /
      DATA TB( 2)/            +0.280606                     /,
     :     (ITB(I, 2),I=1,4)/             0, +1,  0, +1     /
      DATA TB( 3)/            +0.277693                     /,
     :     (ITB(I, 3),I=1,4)/             0, +1,  0, -1     /
      DATA TB( 4)/            +0.173238                     /,
     :     (ITB(I, 4),I=1,4)/             0,  0, +2, -1     /
      DATA TB( 5)/            +0.055413                     /,
     :     (ITB(I, 5),I=1,4)/             0, -1, +2, +1     /
      DATA TB( 6)/            +0.046272                     /,
     :     (ITB(I, 6),I=1,4)/             0, -1, +2, -1     /
      DATA TB( 7)/            +0.032573                     /,
     :     (ITB(I, 7),I=1,4)/             0,  0, +2, +1     /
      DATA TB( 8)/            +0.017198                     /,
     :     (ITB(I, 8),I=1,4)/             0, +2,  0, +1     /
      DATA TB( 9)/            +0.009267                     /,
     :     (ITB(I, 9),I=1,4)/             0, +1, +2, -1     /
      DATA TB(10)/            +0.008823                     /,
     :     (ITB(I,10),I=1,4)/             0, +2,  0, -1     /
      DATA TB(11)/            +0.008247                     /,
     :     (ITB(I,11),I=1,4)/            -1,  0, +2, -1     /
      DATA TB(12)/            +0.004323                     /,
     :     (ITB(I,12),I=1,4)/             0, -2, +2, -1     /
      DATA TB(13)/            +0.004200                     /,
     :     (ITB(I,13),I=1,4)/             0, +1, +2, +1     /
      DATA TB(14)/            +0.003372                     /,
     :     (ITB(I,14),I=1,4)/            -1,  0, -2, +1     /
      DATA TB(15)/            +0.002472                     /,
     :     (ITB(I,15),I=1,4)/            -1, -1, +2, +1     /
      DATA TB(16)/            +0.002222                     /,
     :     (ITB(I,16),I=1,4)/            -1,  0, +2, +1     /
      DATA TB(17)/            +0.002072                     /,
     :     (ITB(I,17),I=1,4)/            -1, -1, +2, -1     /
      DATA TB(18)/            +0.001877                     /,
     :     (ITB(I,18),I=1,4)/            -1, +1,  0, +1     /
      DATA TB(19)/            +0.001828                     /,
     :     (ITB(I,19),I=1,4)/             0, -1, +4, -1     /
      DATA TB(20)/            -0.001803                     /,
     :     (ITB(I,20),I=1,4)/            +1,  0,  0, +1     /
      DATA TB(21)/            -0.001750                     /,
     :     (ITB(I,21),I=1,4)/             0,  0,  0, +3     /
      DATA TB(22)/            +0.001570                     /,
     :     (ITB(I,22),I=1,4)/            -1, +1,  0, -1     /
      DATA TB(23)/            -0.001487                     /,
     :     (ITB(I,23),I=1,4)/             0,  0, +1, +1     /
      DATA TB(24)/            -0.001481                     /,
     :     (ITB(I,24),I=1,4)/            +1, +1,  0, +1     /
      DATA TB(25)/            +0.001417                     /,
     :     (ITB(I,25),I=1,4)/            -1, -1,  0, +1     /
      DATA TB(26)/            +0.001350                     /,
     :     (ITB(I,26),I=1,4)/            -1,  0,  0, +1     /
      DATA TB(27)/            +0.001330                     /,
     :     (ITB(I,27),I=1,4)/             0,  0, -1, +1     /
      DATA TB(28)/            +0.001106                     /,
     :     (ITB(I,28),I=1,4)/             0, +3,  0, +1     /
      DATA TB(29)/            +0.001020                     /,
     :     (ITB(I,29),I=1,4)/             0,  0, +4, -1     /
*
*  Parallax
*                                         M   M'  D   F
      DATA TP( 1)/            +0.051818                     /,
     :     (ITP(I, 1),I=1,4)/             0, +1,  0,  0     /
      DATA TP( 2)/            +0.009531                     /,
     :     (ITP(I, 2),I=1,4)/             0, -1, +2,  0     /
      DATA TP( 3)/            +0.007843                     /,
     :     (ITP(I, 3),I=1,4)/             0,  0, +2,  0     /
      DATA TP( 4)/            +0.002824                     /,
     :     (ITP(I, 4),I=1,4)/             0, +2,  0,  0     /



*  Whole years & fraction of year, and years since J1900.0
      YI=FLOAT(IY-1900)
      IY4=MOD(MOD(IY,4)+4,4)
      YF=(FLOAT(4*(ID-1/(IY4+1))-IY4-2)+4.0*FD)/1461.0
      T=YI+YF

*  Moon's mean longitude
      ELP=D2R*MOD(ELP0+ELP1I*YF+ELP1F*T,360.0)

*  Sun's mean anomaly
      EM=D2R*MOD(EM0+EM1F*T,360.0)

*  Moon's mean anomaly
      EMP=D2R*MOD(EMP0+EMP1I*YF+EMP1F*T,360.0)

*  Moon's mean elongation
      D=D2R*MOD(D0+D1I*YF+D1F*T,360.0)

*  Mean distance of the moon from its ascending node
      F=D2R*MOD(F0+F1I*YF+F1F*T,360.0)

*  Longitude
      EL=0.0
      ELD=0.0
      DO N=39,1,-1
         COEFF=TL(N)
         CEM=FLOAT(ITL(1,N))
         CEMP=FLOAT(ITL(2,N))
         CD=FLOAT(ITL(3,N))
         CF=FLOAT(ITL(4,N))
         THETA=CEM*EM+CEMP*EMP+CD*D+CF*F
         THETAD=CEM*EM1+CEMP*EMP1+CD*D1+CF*F1
         EL=EL+COEFF*SIN(THETA)
         ELD=ELD+COEFF*COS(THETA)*THETAD
      END DO
      EL=EL*D2R+ELP
      ELD=RATCON*(ELD+ELP1/D2R)

*  Latitude
      B=0.0
      BD=0.0
      DO N=29,1,-1
         COEFF=TB(N)
         CEM=FLOAT(ITB(1,N))
         CEMP=FLOAT(ITB(2,N))
         CD=FLOAT(ITB(3,N))
         CF=FLOAT(ITB(4,N))
         THETA=CEM*EM+CEMP*EMP+CD*D+CF*F
         THETAD=CEM*EM1+CEMP*EMP1+CD*D1+CF*F1
         B=B+COEFF*SIN(THETA)
         BD=BD+COEFF*COS(THETA)*THETAD
      END DO
      B=B*D2R
      BD=RATCON*BD

*  Parallax
      P=0.0
      PD=0.0
      DO N=4,1,-1
         COEFF=TP(N)
         CEM=FLOAT(ITP(1,N))
         CEMP=FLOAT(ITP(2,N))
         CD=FLOAT(ITP(3,N))
         CF=FLOAT(ITP(4,N))
         THETA=CEM*EM+CEMP*EMP+CD*D+CF*F
         THETAD=CEM*EM1+CEMP*EMP1+CD*D1+CF*F1
         P=P+COEFF*COS(THETA)
         PD=PD-COEFF*SIN(THETA)*THETAD
      END DO
      P=(P+0.950724)*D2R
      PD=RATCON*PD

*  Transform parallax to distance (AU, AU/sec)
      SP=SIN(P)
      R=ERADAU/SP
      RD=-R*PD/SP

*  Longitude, latitude to x,y,z (AU)
      CALL sla_CS2C6(EL,B,R,ELD,BD,RD,V)

*  Mean obliquity
      EPS=D2R*(23.45229-0.00013*T)
      SINEPS=SIN(EPS)
      COSEPS=COS(EPS)

*  Rotate Moon position and velocity into equatorial system
      PV(1)=V(1)
      PV(2)=V(2)*COSEPS-V(3)*SINEPS
      PV(3)=V(2)*SINEPS+V(3)*COSEPS
      PV(4)=V(4)
      PV(5)=V(5)*COSEPS-V(6)*SINEPS
      PV(6)=V(5)*SINEPS+V(6)*COSEPS

      END
