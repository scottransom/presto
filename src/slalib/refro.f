      SUBROUTINE sla_REFRO ( ZOBS, HM, TDK, PMB, RH, WL, PHI, TLR,
     :                       EPS, REF )
*+
*     - - - - - -
*      R E F R O
*     - - - - - -
*
*  Atmospheric refraction for radio and optical/IR wavelengths.
*
*  Given:
*    ZOBS    d  observed zenith distance of the source (radian)
*    HM      d  height of the observer above sea level (metre)
*    TDK     d  ambient temperature at the observer (K)
*    PMB     d  pressure at the observer (millibar)
*    RH      d  relative humidity at the observer (range 0-1)
*    WL      d  effective wavelength of the source (micrometre)
*    PHI     d  latitude of the observer (radian, astronomical)
*    TLR     d  temperature lapse rate in the troposphere (K/metre)
*    EPS     d  precision required to terminate iteration (radian)
*
*  Returned:
*    REF     d  refraction: in vacuo ZD minus observed ZD (radian)
*
*  Notes:
*
*  1  A suggested value for the TLR argument is 0.0065D0.  The
*     refraction is significantly affected by TLR, and if studies
*     of the local atmosphere have been carried out a better TLR
*     value may be available.  The sign of the supplied TLR value
*     is ignored.
*
*  2  A suggested value for the EPS argument is 1D-8.  The result is
*     usually at least two orders of magnitude more computationally
*     precise than the supplied EPS value.
*
*  3  The routine computes the refraction for zenith distances up
*     to and a little beyond 90 deg using the method of Hohenkerk
*     and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted
*     in the Explanatory Supplement, 1992 edition - see section 3.281).
*
*  4  The code is a development of the optical/IR refraction subroutine
*     AREF of C.Hohenkerk (HMNAO, September 1984), with extensions to
*     support the radio case.  Apart from merely cosmetic changes, the
*     following modifications to the original HMNAO optical/IR refraction
*     code have been made:
*
*     .  The angle arguments have been changed to radians.
*
*     .  Any value of ZOBS is allowed (see note 6, below).
*
*     .  Other argument values have been limited to safe values.
*
*     .  Murray's values for the gas constants have been used
*        (Vectorial Astrometry, Adam Hilger, 1983).
*
*     .  The numerical integration phase has been rearranged for
*        extra clarity.
*
*     .  A better model for Ps(T) has been adopted (taken from
*        Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
*
*     .  More accurate expressions for Pwo have been adopted
*        (again from Gill 1982).
*
*     .  The formula for the water vapour pressure, given the
*        saturation pressure and the relative humidity, is from
*        Crane (1976), expression 2.5.5.
*
*     .  Provision for radio wavelengths has been added using
*        expressions devised by A.T.Sinclair, RGO (private
*        communication 1989).  The refractivity model currently
*        used is from J.M.Rueger, "Refractive Index Formulae for
*        Electronic Distance Measurement with Radio and Millimetre
*        Waves", in Unisurv Report S-68 (2002), School of Surveying
*        and Spatial Information Systems, University of New South
*        Wales, Sydney, Australia.
*
*     .  The optical refractivity for dry air is from Resolution 3 of
*        the International Association of Geodesy adopted at the XXIIth
*        General Assembly in Birmingham, UK, 1999.
*
*     .  Various small changes have been made to gain speed.
*
*  5  The radio refraction is chosen by specifying WL > 100 micrometres.
*     Because the algorithm takes no account of the ionosphere, the
*     accuracy deteriorates at low frequencies, below about 30 MHz.
*
*  6  Before use, the value of ZOBS is expressed in the range +/- pi.
*     If this ranged ZOBS is -ve, the result REF is computed from its
*     absolute value before being made -ve to match.  In addition, if
*     it has an absolute value greater than 93 deg, a fixed REF value
*     equal to the result for ZOBS = 93 deg is returned, appropriately
*     signed.
*
*  7  As in the original Hohenkerk and Sinclair algorithm, fixed values
*     of the water vapour polytrope exponent, the height of the
*     tropopause, and the height at which refraction is negligible are
*     used.
*
*  8  The radio refraction has been tested against work done by
*     Iain Coulson, JACH, (private communication 1995) for the
*     James Clerk Maxwell Telescope, Mauna Kea.  For typical conditions,
*     agreement at the 0.1 arcsec level is achieved for moderate ZD,
*     worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg.  At hot and
*     humid sea-level sites the accuracy will not be as good.
*
*  9  It should be noted that the relative humidity RH is formally
*     defined in terms of "mixing ratio" rather than pressures or
*     densities as is often stated.  It is the mass of water per unit
*     mass of dry air divided by that for saturated air at the same
*     temperature and pressure (see Gill 1982).
*
*  10 The algorithm is designed for observers in the troposphere.  The
*     supplied temperature, pressure and lapse rate are assumed to be
*     for a point in the troposphere and are used to define a model
*     atmosphere with the tropopause at 11km altitude and a constant
*     temperature above that.  However, in practice, the refraction
*     values returned for stratospheric observers, at altitudes up to
*     25km, are quite usable.
*
*  Called:  sla_DRANGE, sla__ATMT, sla__ATMS
*
*  Last revision:   5 December 2005
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

      DOUBLE PRECISION ZOBS,HM,TDK,PMB,RH,WL,PHI,TLR,EPS,REF

*
*  Fixed parameters
*
      DOUBLE PRECISION D93,GCR,DMD,DMW,S,DELTA,HT,HS
      INTEGER ISMAX
*  93 degrees in radians
      PARAMETER (D93=1.623156204D0)
*  Universal gas constant
      PARAMETER (GCR=8314.32D0)
*  Molecular weight of dry air
      PARAMETER (DMD=28.9644D0)
*  Molecular weight of water vapour
      PARAMETER (DMW=18.0152D0)
*  Mean Earth radius (metre)
      PARAMETER (S=6378120D0)
*  Exponent of temperature dependence of water vapour pressure
      PARAMETER (DELTA=18.36D0)
*  Height of tropopause (metre)
      PARAMETER (HT=11000D0)
*  Upper limit for refractive effects (metre)
      PARAMETER (HS=80000D0)
*  Numerical integration: maximum number of strips.
      PARAMETER (ISMAX=16384)

      INTEGER IS,K,N,I,J
      LOGICAL OPTIC,LOOP
      DOUBLE PRECISION ZOBS1,ZOBS2,HMOK,TDKOK,PMBOK,RHOK,WLOK,ALPHA,
     :                 TOL,WLSQ,GB,A,GAMAL,GAMMA,GAMM2,DELM2,
     :                 TDC,PSAT,PWO,W,
     :                 C1,C2,C3,C4,C5,C6,R0,TEMPO,DN0,RDNDR0,SK0,F0,
     :                 RT,TT,DNT,RDNDRT,SINE,ZT,FT,DNTS,RDNDRP,ZTS,FTS,
     :                 RS,DNS,RDNDRS,ZS,FS,REFOLD,Z0,ZRANGE,FB,FF,FO,FE,
     :                 H,R,SZ,RG,DR,TG,DN,RDNDR,T,F,REFP,REFT

      DOUBLE PRECISION sla_DRANGE

*  The refraction integrand
      DOUBLE PRECISION REFI
      REFI(DN,RDNDR) = RDNDR/(DN+RDNDR)



*  Transform ZOBS into the normal range.
      ZOBS1 = sla_DRANGE(ZOBS)
      ZOBS2 = MIN(ABS(ZOBS1),D93)

*  Keep other arguments within safe bounds.
      HMOK = MIN(MAX(HM,-1D3),HS)
      TDKOK = MIN(MAX(TDK,100D0),500D0)
      PMBOK = MIN(MAX(PMB,0D0),10000D0)
      RHOK = MIN(MAX(RH,0D0),1D0)
      WLOK = MAX(WL,0.1D0)
      ALPHA = MIN(MAX(ABS(TLR),0.001D0),0.01D0)

*  Tolerance for iteration.
      TOL = MIN(MAX(ABS(EPS),1D-12),0.1D0)/2D0

*  Decide whether optical/IR or radio case - switch at 100 microns.
      OPTIC = WLOK.LE.100D0

*  Set up model atmosphere parameters defined at the observer.
      WLSQ = WLOK*WLOK
      GB = 9.784D0*(1D0-0.0026D0*COS(PHI+PHI)-0.00000028D0*HMOK)
      IF (OPTIC) THEN
         A = (287.6155D0+(1.62887D0+0.01360D0/WLSQ)/WLSQ)
     :                                              *273.15D-6/1013.25D0
      ELSE
         A = 77.6890D-6
      END IF
      GAMAL = (GB*DMD)/GCR
      GAMMA = GAMAL/ALPHA
      GAMM2 = GAMMA-2D0
      DELM2 = DELTA-2D0
      TDC = TDKOK-273.15D0
      PSAT = 10D0**((0.7859D0+0.03477D0*TDC)/(1D0+0.00412D0*TDC))*
     :                                (1D0+PMBOK*(4.5D-6+6D-10*TDC*TDC))
      IF (PMBOK.GT.0D0) THEN
         PWO = RHOK*PSAT/(1D0-(1D0-RHOK)*PSAT/PMBOK)
      ELSE
         PWO = 0D0
      END IF
      W = PWO*(1D0-DMW/DMD)*GAMMA/(DELTA-GAMMA)
      C1 = A*(PMBOK+W)/TDKOK
      IF (OPTIC) THEN
         C2 = (A*W+11.2684D-6*PWO)/TDKOK
      ELSE
         C2 = (A*W+6.3938D-6*PWO)/TDKOK
      END IF
      C3 = (GAMMA-1D0)*ALPHA*C1/TDKOK
      C4 = (DELTA-1D0)*ALPHA*C2/TDKOK
      IF (OPTIC) THEN
         C5 = 0D0
         C6 = 0D0
      ELSE
         C5 = 375463D-6*PWO/TDKOK
         C6 = C5*DELM2*ALPHA/(TDKOK*TDKOK)
      END IF

*  Conditions at the observer.
      R0 = S+HMOK
      CALL sla__ATMT(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,
     :                                              R0,TEMPO,DN0,RDNDR0)
      SK0 = DN0*R0*SIN(ZOBS2)
      F0 = REFI(DN0,RDNDR0)

*  Conditions in the troposphere at the tropopause.
      RT = S+MAX(HT,HMOK)
      CALL sla__ATMT(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,
     :                                                 RT,TT,DNT,RDNDRT)
      SINE = SK0/(RT*DNT)
      ZT = ATAN2(SINE,SQRT(MAX(1D0-SINE*SINE,0D0)))
      FT = REFI(DNT,RDNDRT)

*  Conditions in the stratosphere at the tropopause.
      CALL sla__ATMS(RT,TT,DNT,GAMAL,RT,DNTS,RDNDRP)
      SINE = SK0/(RT*DNTS)
      ZTS = ATAN2(SINE,SQRT(MAX(1D0-SINE*SINE,0D0)))
      FTS = REFI(DNTS,RDNDRP)

*  Conditions at the stratosphere limit.
      RS = S+HS
      CALL sla__ATMS(RT,TT,DNT,GAMAL,RS,DNS,RDNDRS)
      SINE = SK0/(RS*DNS)
      ZS = ATAN2(SINE,SQRT(MAX(1D0-SINE*SINE,0D0)))
      FS = REFI(DNS,RDNDRS)

*  Variable initialization to avoid compiler warning.
      REFT = 0D0

*  Integrate the refraction integral in two parts;  first in the
*  troposphere (K=1), then in the stratosphere (K=2).

      DO K = 1,2

*     Initialize previous refraction to ensure at least two iterations.
         REFOLD = 1D0

*     Start off with 8 strips.
         IS = 8

*     Start Z, Z range, and start and end values.
         IF (K.EQ.1) THEN
            Z0 = ZOBS2
            ZRANGE = ZT-Z0
            FB = F0
            FF = FT
         ELSE
            Z0 = ZTS
            ZRANGE = ZS-Z0
            FB = FTS
            FF = FS
         END IF

*     Sums of odd and even values.
         FO = 0D0
         FE = 0D0

*     First time through the loop we have to do every point.
         N = 1

*     Start of iteration loop (terminates at specified precision).
         LOOP = .TRUE.
         DO WHILE (LOOP)

*        Strip width.
            H = ZRANGE/DBLE(IS)

*        Initialize distance from Earth centre for quadrature pass.
            IF (K.EQ.1) THEN
               R = R0
            ELSE
               R = RT
            END IF

*        One pass (no need to compute evens after first time).
            DO I=1,IS-1,N

*           Sine of observed zenith distance.
               SZ = SIN(Z0+H*DBLE(I))

*           Find R (to the nearest metre, maximum four iterations).
               IF (SZ.GT.1D-20) THEN
                  W = SK0/SZ
                  RG = R
                  DR = 1D6
                  J = 0
                  DO WHILE (ABS(DR).GT.1D0.AND.J.LT.4)
                     J=J+1
                     IF (K.EQ.1) THEN
                        CALL sla__ATMT(R0,TDKOK,ALPHA,GAMM2,DELM2,
     :                                 C1,C2,C3,C4,C5,C6,RG,TG,DN,RDNDR)
                     ELSE
                        CALL sla__ATMS(RT,TT,DNT,GAMAL,RG,DN,RDNDR)
                     END IF
                     DR = (RG*DN-W)/(DN+RDNDR)
                     RG = RG-DR
                  END DO
                  R = RG
               END IF

*           Find the refractive index and integrand at R.
               IF (K.EQ.1) THEN
                  CALL sla__ATMT(R0,TDKOK,ALPHA,GAMM2,DELM2,
     :                                   C1,C2,C3,C4,C5,C6,R,T,DN,RDNDR)
               ELSE
                  CALL sla__ATMS(RT,TT,DNT,GAMAL,R,DN,RDNDR)
               END IF
               F = REFI(DN,RDNDR)

*           Accumulate odd and (first time only) even values.
               IF (N.EQ.1.AND.MOD(I,2).EQ.0) THEN
                  FE = FE+F
               ELSE
                  FO = FO+F
               END IF
            END DO

*        Evaluate the integrand using Simpson's Rule.
            REFP = H*(FB+4D0*FO+2D0*FE+FF)/3D0

*        Has the required precision been achieved (or can't be)?
            IF (ABS(REFP-REFOLD).GT.TOL.AND.IS.LT.ISMAX) THEN

*           No: prepare for next iteration.

*           Save current value for convergence test.
               REFOLD = REFP

*           Double the number of strips.
               IS = IS+IS

*           Sum of all current values = sum of next pass's even values.
               FE = FE+FO

*           Prepare for new odd values.
               FO = 0D0

*           Skip even values next time.
               N = 2
            ELSE

*           Yes: save troposphere component and terminate the loop.
               IF (K.EQ.1) REFT = REFP
               LOOP = .FALSE.
            END IF
         END DO
      END DO

*  Result.
      REF = REFT+REFP
      IF (ZOBS1.LT.0D0) REF = -REF

      END
