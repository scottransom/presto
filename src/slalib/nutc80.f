      SUBROUTINE sla_NUTC80 (DATE, DPSI, DEPS, EPS0)
*+
*     - - - - - - -
*      N U T C 8 0
*     - - - - - - -
*
*  Nutation:  longitude & obliquity components and mean obliquity,
*  using the IAU 1980 theory (double precision)
*
*  Given:
*     DATE        d     TDB (loosely ET) as Modified Julian Date
*                                            (JD-2400000.5)
*  Returned:
*     DPSI,DEPS   d     nutation in longitude,obliquity
*     EPS0        d     mean obliquity
*
*  Called:  sla_DRANGE
*
*  Notes:
*
*  1  Earth attitude predictions made by combining the present nutation
*     model with IAU 1976 precession are accurate to 0.35 arcsec over
*     the interval 1900-2100.  (The accuracy is much better near the
*     middle of the interval.)
*
*  2  The sla_NUTC routine is the equivalent of the present routine
*     but using the Shirai & Fukushima 2001 forced nutation theory
*     (SF2001).  The newer theory is more accurate than IAU 1980,
*     within 1 mas (with respect to the ICRF) for a few decades around
*     2000.  The improvement is mainly because of the corrections to the
*     IAU 1976 precession that the SF2001 theory includes.
*
*  References:
*     Final report of the IAU Working Group on Nutation,
*      chairman P.K.Seidelmann, 1980.
*     Kaplan,G.H., 1981, USNO circular no. 163, pA3-6.
*
*  P.T.Wallace   Starlink   7 October 2001
*
*  Copyright (C) 2001 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION DATE,DPSI,DEPS,EPS0

      DOUBLE PRECISION T2AS,AS2R,U2R
      DOUBLE PRECISION T,EL,EL2,EL3
      DOUBLE PRECISION ELP,ELP2
      DOUBLE PRECISION F,F2,F4
      DOUBLE PRECISION D,D2,D4
      DOUBLE PRECISION OM,OM2
      DOUBLE PRECISION DP,DE
      DOUBLE PRECISION A

      DOUBLE PRECISION sla_DRANGE


*  Turns to arc seconds
      PARAMETER (T2AS=1296000D0)
*  Arc seconds to radians
      PARAMETER (AS2R=0.484813681109535994D-5)
*  Units of 0.0001 arcsec to radians
      PARAMETER (U2R=AS2R/1D4)




*  Interval between basic epoch J2000.0 and current epoch (JC)
      T=(DATE-51544.5D0)/36525D0

*
*  FUNDAMENTAL ARGUMENTS in the FK5 reference system
*

*  Mean longitude of the Moon minus mean longitude of the Moon's perigee
      EL=sla_DRANGE(AS2R*(485866.733D0+(1325D0*T2AS+715922.633D0
     :                   +(31.310D0+0.064D0*T)*T)*T))

*  Mean longitude of the Sun minus mean longitude of the Sun's perigee
      ELP=sla_DRANGE(AS2R*(1287099.804D0+(99D0*T2AS+1292581.224D0
     :                   +(-0.577D0-0.012D0*T)*T)*T))

*  Mean longitude of the Moon minus mean longitude of the Moon's node
      F=sla_DRANGE(AS2R*(335778.877D0+(1342D0*T2AS+295263.137D0
     :                   +(-13.257D0+0.011D0*T)*T)*T))

*  Mean elongation of the Moon from the Sun
      D=sla_DRANGE(AS2R*(1072261.307D0+(1236D0*T2AS+1105601.328D0
     :                   +(-6.891D0+0.019D0*T)*T)*T))

*  Longitude of the mean ascending node of the lunar orbit on the
*   ecliptic, measured from the mean equinox of date
      OM=sla_DRANGE(AS2R*(450160.280D0+(-5D0*T2AS-482890.539D0
     :                   +(7.455D0+0.008D0*T)*T)*T))

*  Multiples of arguments
      EL2=EL+EL
      EL3=EL2+EL
      ELP2=ELP+ELP
      F2=F+F
      F4=F2+F2
      D2=D+D
      D4=D2+D2
      OM2=OM+OM


*
*  SERIES FOR THE NUTATION
*
      DP=0D0
      DE=0D0

*  106
      DP=DP+SIN(ELP+D)
*  105
      DP=DP-SIN(F2+D4+OM2)
*  104
      DP=DP+SIN(EL2+D2)
*  103
      DP=DP-SIN(EL-F2+D2)
*  102
      DP=DP-SIN(EL+ELP-D2+OM)
*  101
      DP=DP-SIN(-ELP+F2+OM)
*  100
      DP=DP-SIN(EL-F2-D2)
*  99
      DP=DP-SIN(ELP+D2)
*  98
      DP=DP-SIN(F2-D+OM2)
*  97
      DP=DP-SIN(-F2+OM)
*  96
      DP=DP+SIN(-EL-ELP+D2+OM)
*  95
      DP=DP+SIN(ELP+F2+OM)
*  94
      DP=DP-SIN(EL+F2-D2)
*  93
      DP=DP+SIN(EL3+F2-D2+OM2)
*  92
      DP=DP+SIN(F4-D2+OM2)
*  91
      DP=DP-SIN(EL+D2+OM)
*  90
      DP=DP-SIN(EL2+F2+D2+OM2)
*  89
      A=EL2+F2-D2+OM
      DP=DP+SIN(A)
      DE=DE-COS(A)
*  88
      DP=DP+SIN(EL-ELP-D2)
*  87
      DP=DP+SIN(-EL+F4+OM2)
*  86
      A=-EL2+F2+D4+OM2
      DP=DP-SIN(A)
      DE=DE+COS(A)
*  85
      A=EL+F2+D2+OM
      DP=DP-SIN(A)
      DE=DE+COS(A)
*  84
      A=EL+ELP+F2-D2+OM2
      DP=DP+SIN(A)
      DE=DE-COS(A)
*  83
      DP=DP-SIN(EL2-D4)
*  82
      A=-EL+F2+D4+OM2
      DP=DP-2D0*SIN(A)
      DE=DE+COS(A)
*  81
      A=-EL2+F2+D2+OM2
      DP=DP+SIN(A)
      DE=DE-COS(A)
*  80
      DP=DP-SIN(EL-D4)
*  79
      A=-EL+OM2
      DP=DP+SIN(A)
      DE=DE-COS(A)
*  78
      A=F2+D+OM2
      DP=DP+2D0*SIN(A)
      DE=DE-COS(A)
*  77
      DP=DP+2D0*SIN(EL3)
*  76
      A=EL+OM2
      DP=DP-2D0*SIN(A)
      DE=DE+COS(A)
*  75
      A=EL2+OM
      DP=DP+2D0*SIN(A)
      DE=DE-COS(A)
*  74
      A=-EL+F2-D2+OM
      DP=DP-2D0*SIN(A)
      DE=DE+COS(A)
*  73
      A=EL+ELP+F2+OM2
      DP=DP+2D0*SIN(A)
      DE=DE-COS(A)
*  72
      A=-ELP+F2+D2+OM2
      DP=DP-3D0*SIN(A)
      DE=DE+COS(A)
*  71
      A=EL3+F2+OM2
      DP=DP-3D0*SIN(A)
      DE=DE+COS(A)
*  70
      A=-EL2+OM
      DP=DP-2D0*SIN(A)
      DE=DE+COS(A)
*  69
      A=-EL-ELP+F2+D2+OM2
      DP=DP-3D0*SIN(A)
      DE=DE+COS(A)
*  68
      A=EL-ELP+F2+OM2
      DP=DP-3D0*SIN(A)
      DE=DE+COS(A)
*  67
      DP=DP+3D0*SIN(EL+F2)
*  66
      DP=DP-3D0*SIN(EL+ELP)
*  65
      DP=DP-4D0*SIN(D)
*  64
      DP=DP+4D0*SIN(EL-F2)
*  63
      DP=DP-4D0*SIN(ELP-D2)
*  62
      A=EL2+F2+OM
      DP=DP-5D0*SIN(A)
      DE=DE+3D0*COS(A)
*  61
      DP=DP+5D0*SIN(EL-ELP)
*  60
      A=-D2+OM
      DP=DP-5D0*SIN(A)
      DE=DE+3D0*COS(A)
*  59
      A=EL+F2-D2+OM
      DP=DP+6D0*SIN(A)
      DE=DE-3D0*COS(A)
*  58
      A=F2+D2+OM
      DP=DP-7D0*SIN(A)
      DE=DE+3D0*COS(A)
*  57
      A=D2+OM
      DP=DP-6D0*SIN(A)
      DE=DE+3D0*COS(A)
*  56
      A=EL2+F2-D2+OM2
      DP=DP+6D0*SIN(A)
      DE=DE-3D0*COS(A)
*  55
      DP=DP+6D0*SIN(EL+D2)
*  54
      A=EL+F2+D2+OM2
      DP=DP-8D0*SIN(A)
      DE=DE+3D0*COS(A)
*  53
      A=-ELP+F2+OM2
      DP=DP-7D0*SIN(A)
      DE=DE+3D0*COS(A)
*  52
      A=ELP+F2+OM2
      DP=DP+7D0*SIN(A)
      DE=DE-3D0*COS(A)
*  51
      DP=DP-7D0*SIN(EL+ELP-D2)
*  50
      A=-EL+F2+D2+OM
      DP=DP-10D0*SIN(A)
      DE=DE+5D0*COS(A)
*  49
      A=EL-D2+OM
      DP=DP-13D0*SIN(A)
      DE=DE+7D0*COS(A)
*  48
      A=-EL+D2+OM
      DP=DP+16D0*SIN(A)
      DE=DE-8D0*COS(A)
*  47
      A=-EL+F2+OM
      DP=DP+21D0*SIN(A)
      DE=DE-10D0*COS(A)
*  46
      DP=DP+26D0*SIN(F2)
      DE=DE-COS(F2)
*  45
      A=EL2+F2+OM2
      DP=DP-31D0*SIN(A)
      DE=DE+13D0*COS(A)
*  44
      A=EL+F2-D2+OM2
      DP=DP+29D0*SIN(A)
      DE=DE-12D0*COS(A)
*  43
      DP=DP+29D0*SIN(EL2)
      DE=DE-COS(EL2)
*  42
      A=F2+D2+OM2
      DP=DP-38D0*SIN(A)
      DE=DE+16D0*COS(A)
*  41
      A=EL+F2+OM
      DP=DP-51D0*SIN(A)
      DE=DE+27D0*COS(A)
*  40
      A=-EL+F2+D2+OM2
      DP=DP-59D0*SIN(A)
      DE=DE+26D0*COS(A)
*  39
      A=-EL+OM
      DP=DP+(-58D0-0.1D0*T)*SIN(A)
      DE=DE+32D0*COS(A)
*  38
      A=EL+OM
      DP=DP+(63D0+0.1D0*T)*SIN(A)
      DE=DE-33D0*COS(A)
*  37
      DP=DP+63D0*SIN(D2)
      DE=DE-2D0*COS(D2)
*  36
      A=-EL+F2+OM2
      DP=DP+123D0*SIN(A)
      DE=DE-53D0*COS(A)
*  35
      A=EL-D2
      DP=DP-158D0*SIN(A)
      DE=DE-COS(A)
*  34
      A=EL+F2+OM2
      DP=DP-301D0*SIN(A)
      DE=DE+(129D0-0.1D0*T)*COS(A)
*  33
      A=F2+OM
      DP=DP+(-386D0-0.4D0*T)*SIN(A)
      DE=DE+200D0*COS(A)
*  32
      DP=DP+(712D0+0.1D0*T)*SIN(EL)
      DE=DE-7D0*COS(EL)
*  31
      A=F2+OM2
      DP=DP+(-2274D0-0.2D0*T)*SIN(A)
      DE=DE+(977D0-0.5D0*T)*COS(A)
*  30
      DP=DP-SIN(ELP+F2-D2)
*  29
      DP=DP+SIN(-EL+D+OM)
*  28
      DP=DP+SIN(ELP+OM2)
*  27
      DP=DP-SIN(ELP-F2+D2)
*  26
      DP=DP+SIN(-F2+D2+OM)
*  25
      DP=DP+SIN(EL2+ELP-D2)
*  24
      DP=DP-4D0*SIN(EL-D)
*  23
      A=ELP+F2-D2+OM
      DP=DP+4D0*SIN(A)
      DE=DE-2D0*COS(A)
*  22
      A=EL2-D2+OM
      DP=DP+4D0*SIN(A)
      DE=DE-2D0*COS(A)
*  21
      A=-ELP+F2-D2+OM
      DP=DP-5D0*SIN(A)
      DE=DE+3D0*COS(A)
*  20
      A=-EL2+D2+OM
      DP=DP-6D0*SIN(A)
      DE=DE+3D0*COS(A)
*  19
      A=-ELP+OM
      DP=DP-12D0*SIN(A)
      DE=DE+6D0*COS(A)
*  18
      A=ELP2+F2-D2+OM2
      DP=DP+(-16D0+0.1D0*T)*SIN(A)
      DE=DE+7D0*COS(A)
*  17
      A=ELP+OM
      DP=DP-15D0*SIN(A)
      DE=DE+9D0*COS(A)
*  16
      DP=DP+(17D0-0.1D0*T)*SIN(ELP2)
*  15
      DP=DP-22D0*SIN(F2-D2)
*  14
      A=EL2-D2
      DP=DP+48D0*SIN(A)
      DE=DE+COS(A)
*  13
      A=F2-D2+OM
      DP=DP+(129D0+0.1D0*T)*SIN(A)
      DE=DE-70D0*COS(A)
*  12
      A=-ELP+F2-D2+OM2
      DP=DP+(217D0-0.5D0*T)*SIN(A)
      DE=DE+(-95D0+0.3D0*T)*COS(A)
*  11
      A=ELP+F2-D2+OM2
      DP=DP+(-517D0+1.2D0*T)*SIN(A)
      DE=DE+(224D0-0.6D0*T)*COS(A)
*  10
      DP=DP+(1426D0-3.4D0*T)*SIN(ELP)
      DE=DE+(54D0-0.1D0*T)*COS(ELP)
*  9
      A=F2-D2+OM2
      DP=DP+(-13187D0-1.6D0*T)*SIN(A)
      DE=DE+(5736D0-3.1D0*T)*COS(A)
*  8
      DP=DP+SIN(EL2-F2+OM)
*  7
      A=-ELP2+F2-D2+OM
      DP=DP-2D0*SIN(A)
      DE=DE+1D0*COS(A)
*  6
      DP=DP-3D0*SIN(EL-ELP-D)
*  5
      A=-EL2+F2+OM2
      DP=DP-3D0*SIN(A)
      DE=DE+1D0*COS(A)
*  4
      DP=DP+11D0*SIN(EL2-F2)
*  3
      A=-EL2+F2+OM
      DP=DP+46D0*SIN(A)
      DE=DE-24D0*COS(A)
*  2
      DP=DP+(2062D0+0.2D0*T)*SIN(OM2)
      DE=DE+(-895D0+0.5D0*T)*COS(OM2)
*  1
      DP=DP+(-171996D0-174.2D0*T)*SIN(OM)
      DE=DE+(92025D0+8.9D0*T)*COS(OM)

*  Convert results to radians
      DPSI=DP*U2R
      DEPS=DE*U2R

*  Mean obliquity
      EPS0=AS2R*(84381.448D0+
     :           (-46.8150D0+
     :           (-0.00059D0+
     :           0.001813D0*T)*T)*T)

      END
