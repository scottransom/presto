      SUBROUTINE sla_REFZ (ZU, REFA, REFB, ZR)
*+
*     - - - - -
*      R E F Z
*     - - - - -
*
*  Adjust an unrefracted zenith distance to include the effect of
*  atmospheric refraction, using the simple A tan Z + B tan**3 Z
*  model (plus special handling for large ZDs).
*
*  Given:
*    ZU    dp    unrefracted zenith distance of the source (radian)
*    REFA  dp    tan Z coefficient (radian)
*    REFB  dp    tan**3 Z coefficient (radian)
*
*  Returned:
*    ZR    dp    refracted zenith distance (radian)
*
*  Notes:
*
*  1  This routine applies the adjustment for refraction in the
*     opposite sense to the usual one - it takes an unrefracted
*     (in vacuo) position and produces an observed (refracted)
*     position, whereas the A tan Z + B tan**3 Z model strictly
*     applies to the case where an observed position is to have the
*     refraction removed.  The unrefracted to refracted case is
*     harder, and requires an inverted form of the text-book
*     refraction models;  the formula used here is based on the
*     Newton-Raphson method.  For the utmost numerical consistency
*     with the refracted to unrefracted model, two iterations are
*     carried out, achieving agreement at the 1D-11 arcseconds level
*     for a ZD of 80 degrees.  The inherent accuracy of the model
*     is, of course, far worse than this - see the documentation for
*     sla_REFCO for more information.
*
*  2  At ZD 83 degrees, the rapidly-worsening A tan Z + B tan^3 Z
*     model is abandoned and an empirical formula takes over.  For
*     optical/IR wavelengths, over a wide range of observer heights and
*     corresponding temperatures and pressures, the following levels of
*     accuracy (arcsec, worst case) are achieved, relative to numerical
*     integration through a model atmosphere:
*
*              ZR    error
*
*              80      0.7
*              81      1.3
*              82      2.4
*              83      4.7
*              84      6.2
*              85      6.4
*              86      8
*              87     10
*              88     15
*              89     30
*              90     60
*              91    150         } relevant only to
*              92    400         } high-elevation sites
*
*     For radio wavelengths the errors are typically 50% larger than
*     the optical figures and by ZD 85 deg are twice as bad, worsening
*     rapidly below that.  To maintain 1 arcsec accuracy down to ZD=85
*     at the Green Bank site, Condon (2004) has suggested amplifying
*     the amount of refraction predicted by sla_REFZ below 10.8 deg
*     elevation by the factor (1+0.00195*(10.8-E_t)), where E_t is the
*     unrefracted elevation in degrees.
*
*     The high-ZD model is scaled to match the normal model at the
*     transition point;  there is no glitch.
*
*  3  Beyond 93 deg zenith distance, the refraction is held at its
*     93 deg value.
*
*  4  See also the routine sla_REFV, which performs the adjustment in
*     Cartesian Az/El coordinates, and with the emphasis on speed
*     rather than numerical accuracy.
*
*  Reference:
*
*     Condon,J.J., Refraction Corrections for the GBT, PTCS/PN/35.2,
*     NRAO Green Bank, 2004.
*
*  P.T.Wallace   Starlink   9 April 2004
*
*  Copyright (C) 2004 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION ZU,REFA,REFB,ZR

*  Radians to degrees
      DOUBLE PRECISION R2D
      PARAMETER (R2D=57.29577951308232D0)

*  Largest usable ZD (deg)
      DOUBLE PRECISION D93
      PARAMETER (D93=93D0)

*  Coefficients for high ZD model (used beyond ZD 83 deg)
      DOUBLE PRECISION C1,C2,C3,C4,C5
      PARAMETER (C1=+0.55445D0,
     :           C2=-0.01133D0,
     :           C3=+0.00202D0,
     :           C4=+0.28385D0,
     :           C5=+0.02390D0)

*  ZD at which one model hands over to the other (radians)
      DOUBLE PRECISION Z83
      PARAMETER (Z83=83D0/R2D)

*  High-ZD-model prediction (deg) for that point
      DOUBLE PRECISION REF83
      PARAMETER (REF83=(C1+C2*7D0+C3*49D0)/(1D0+C4*7D0+C5*49D0))

      DOUBLE PRECISION ZU1,ZL,S,C,T,TSQ,TCU,REF,E,E2



*  Perform calculations for ZU or 83 deg, whichever is smaller
      ZU1 = MIN(ZU,Z83)

*  Functions of ZD
      ZL = ZU1
      S = SIN(ZL)
      C = COS(ZL)
      T = S/C
      TSQ = T*T
      TCU = T*TSQ

*  Refracted ZD (mathematically to better than 1 mas at 70 deg)
      ZL = ZL-(REFA*T+REFB*TCU)/(1D0+(REFA+3D0*REFB*TSQ)/(C*C))

*  Further iteration
      S = SIN(ZL)
      C = COS(ZL)
      T = S/C
      TSQ = T*T
      TCU = T*TSQ
      REF = ZU1-ZL+
     :          (ZL-ZU1+REFA*T+REFB*TCU)/(1D0+(REFA+3D0*REFB*TSQ)/(C*C))

*  Special handling for large ZU
      IF (ZU.GT.ZU1) THEN
         E = 90D0-MIN(D93,ZU*R2D)
         E2 = E*E
         REF = (REF/REF83)*(C1+C2*E+C3*E2)/(1D0+C4*E+C5*E2)
      END IF

*  Return refracted ZD
      ZR = ZU-REF

      END
