      SUBROUTINE sla_REFCOQ ( TDK, PMB, RH, WL, REFA, REFB )
*+
*     - - - - - - -
*      R E F C O Q
*     - - - - - - -
*
*  Determine the constants A and B in the atmospheric refraction
*  model dZ = A tan Z + B tan**3 Z.  This is a fast alternative
*  to the sla_REFCO routine - see notes.
*
*  Z is the "observed" zenith distance (i.e. affected by refraction)
*  and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
*  zenith distance.
*
*  Given:
*    TDK      d      ambient temperature at the observer (K)
*    PMB      d      pressure at the observer (millibar)
*    RH       d      relative humidity at the observer (range 0-1)
*    WL       d      effective wavelength of the source (micrometre)
*
*  Returned:
*    REFA     d      tan Z coefficient (radian)
*    REFB     d      tan**3 Z coefficient (radian)
*
*  The radio refraction is chosen by specifying WL > 100 micrometres.
*
*  Notes:
*
*  1  The model is an approximation, for moderate zenith distances,
*     to the predictions of the sla_REFRO routine.  The approximation
*     is maintained across a range of conditions, and applies to
*     both optical/IR and radio.
*
*  2  The algorithm is a fast alternative to the sla_REFCO routine.
*     The latter calls the sla_REFRO routine itself:  this involves
*     integrations through a model atmosphere, and is costly in
*     processor time.  However, the model which is produced is precisely
*     correct for two zenith distance (45 degrees and about 76 degrees)
*     and at other zenith distances is limited in accuracy only by the
*     A tan Z + B tan**3 Z formulation itself.  The present routine
*     is not as accurate, though it satisfies most practical
*     requirements.
*
*  3  The model omits the effects of (i) height above sea level (apart
*     from the reduced pressure itself), (ii) latitude (i.e. the
*     flattening of the Earth) and (iii) variations in tropospheric
*     lapse rate.
*
*     The model was tested using the following range of conditions:
*
*       lapse rates 0.0055, 0.0065, 0.0075 K/metre
*       latitudes 0, 25, 50, 75 degrees
*       heights 0, 2500, 5000 metres ASL
*       pressures mean for height -10% to +5% in steps of 5%
*       temperatures -10 deg to +20 deg with respect to 280 deg at SL
*       relative humidity 0, 0.5, 1
*       wavelengths 0.4, 0.6, ... 2 micron, + radio
*       zenith distances 15, 45, 75 degrees
*
*     The accuracy with respect to direct use of the sla_REFRO routine
*     was as follows:
*
*                            worst         RMS
*
*       optical/IR           62 mas       8 mas
*       radio               319 mas      49 mas
*
*     For this particular set of conditions:
*
*       lapse rate 0.0065 K/metre
*       latitude 50 degrees
*       sea level
*       pressure 1005 mb
*       temperature 280.15 K
*       humidity 80%
*       wavelength 5740 Angstroms
*
*     the results were as follows:
*
*       ZD        sla_REFRO   sla_REFCOQ  Saastamoinen
*
*       10         10.27        10.27        10.27
*       20         21.19        21.20        21.19
*       30         33.61        33.61        33.60
*       40         48.82        48.83        48.81
*       45         58.16        58.18        58.16
*       50         69.28        69.30        69.27
*       55         82.97        82.99        82.95
*       60        100.51       100.54       100.50
*       65        124.23       124.26       124.20
*       70        158.63       158.68       158.61
*       72        177.32       177.37       177.31
*       74        200.35       200.38       200.32
*       76        229.45       229.43       229.42
*       78        267.44       267.29       267.41
*       80        319.13       318.55       319.10
*
*      deg        arcsec       arcsec       arcsec
*
*     The values for Saastamoinen's formula (which includes terms
*     up to tan^5) are taken from Hohenkerk and Sinclair (1985).
*
*     The results from the much slower but more accurate sla_REFCO
*     routine have not been included in the tabulation as they are
*     identical to those in the sla_REFRO column to the 0.01 arcsec
*     resolution used.
*
*  4  Outlandish input parameters are silently limited to mathematically
*     safe values.  Zero pressure is permissible, and causes zeroes to
*     be returned.
*
*  5  The algorithm draws on several sources, as follows:
*
*     a) The formula for the saturation vapour pressure of water as
*        a function of temperature and temperature is taken from
*        expressions A4.5-A4.7 of Gill (1982).
*
*     b) The formula for the water vapour pressure, given the
*        saturation pressure and the relative humidity, is from
*        Crane (1976), expression 2.5.5.
*
*     c) The refractivity of air is a function of temperature,
*        total pressure, water-vapour pressure and, in the case
*        of optical/IR but not radio, wavelength.  The formulae
*        for the two cases are developed from Hohenkerk & Sinclair
*        (1985) and Rueger (2002).
*
*     The above three items are as used in the sla_REFRO routine.
*
*     d) The formula for beta, the ratio of the scale height of the
*        atmosphere to the geocentric distance of the observer, is
*        an adaption of expression 9 from Stone (1996).  The
*        adaptations, arrived at empirically, consist of (i) a
*        small adjustment to the coefficient and (ii) a humidity
*        term for the radio case only.
*
*     e) The formulae for the refraction constants as a function of
*        n-1 and beta are from Green (1987), expression 4.31.
*
*  References:
*
*     Crane, R.K., Meeks, M.L. (ed), "Refraction Effects in the Neutral
*     Atmosphere", Methods of Experimental Physics: Astrophysics 12B,
*     Academic Press, 1976.
*
*     Gill, Adrian E., "Atmosphere-Ocean Dynamics", Academic Press, 1982.
*
*     Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987.
*
*     Hohenkerk, C.Y., & Sinclair, A.T., NAO Technical Note No. 63, 1985.
*
*     Rueger, J.M., "Refractive Index Formulae for Electronic Distance
*     Measurement with Radio and Millimetre Waves", in Unisurv Report
*     S-68, School of Surveying and Spatial Information Systems,
*     University of New South Wales, Sydney, Australia, 2002.
*
*     Stone, Ronald C., P.A.S.P. 108 1051-1058, 1996.
*
*  Last revision:   2 December 2005
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

      DOUBLE PRECISION TDK,PMB,RH,WL,REFA,REFB

      LOGICAL OPTIC
      DOUBLE PRECISION T,P,R,W,TDC,PS,PW,WLSQ,GAMMA,BETA



*  Decide whether optical/IR or radio case:  switch at 100 microns.
      OPTIC = WL.LE.100D0

*  Restrict parameters to safe values.
      T = MIN(MAX(TDK,100D0),500D0)
      P = MIN(MAX(PMB,0D0),10000D0)
      R = MIN(MAX(RH,0D0),1D0)
      W = MIN(MAX(WL,0.1D0),1D6)

*  Water vapour pressure at the observer.
      IF (P.GT.0D0) THEN
         TDC = T-273.15D0
         PS = 10D0**((0.7859D0+0.03477D0*TDC)/(1D0+0.00412D0*TDC))*
     :                                    (1D0+P*(4.5D-6+6D-10*TDC*TDC))
         PW = R*PS/(1D0-(1D0-R)*PS/P)
      ELSE
         PW = 0D0
      END IF

*  Refractive index minus 1 at the observer.
      IF (OPTIC) THEN
         WLSQ = W*W
         GAMMA = ((77.53484D-6+(4.39108D-7+3.666D-9/WLSQ)/WLSQ)*P
     :                                                 -11.2684D-6*PW)/T
      ELSE
         GAMMA = (77.6890D-6*P-(6.3938D-6-0.375463D0/T)*PW)/T
      END IF

*  Formula for beta adapted from Stone, with empirical adjustments.
      BETA=4.4474D-6*T
      IF (.NOT.OPTIC) BETA=BETA-0.0074D0*PW*BETA

*  Refraction constants from Green.
      REFA = GAMMA*(1D0-BETA)
      REFB = -GAMMA*(BETA-GAMMA/2D0)

      END
