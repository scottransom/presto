      REAL FUNCTION sla_RVLSRD (R2000, D2000)
*+
*     - - - - - - -
*      R V L S R D
*     - - - - - - -
*
*  Velocity component in a given direction due to the Sun's motion
*  with respect to the dynamical Local Standard of Rest.
*
*  (single precision)
*
*  Given:
*     R2000,D2000   r    J2000.0 mean RA,Dec (radians)
*
*  Result:
*     Component of "peculiar" solar motion in direction R2000,D2000 (km/s)
*
*  Sign convention:
*     The result is +ve when the Sun is receding from the given point on
*     the sky.
*
*  Note:  The Local Standard of Rest used here is the "dynamical" LSR,
*         a point in the vicinity of the Sun which is in a circular
*         orbit around the Galactic centre.  The Sun's motion with
*         respect to the dynamical LSR is called the "peculiar" solar
*         motion.
*
*         There is another type of LSR, called a "kinematical" LSR.  A
*         kinematical LSR is the mean standard of rest of specified star
*         catalogues or stellar populations, and several slightly
*         different kinematical LSRs are in use.  The Sun's motion with
*         respect to an agreed kinematical LSR is known as the "standard"
*         solar motion.  To obtain a radial velocity correction with
*         respect to an adopted kinematical LSR use the routine sla_RVLSRK.
*
*  Reference:  Delhaye (1965), in "Stars and Stellar Systems", vol 5,
*              p73.
*
*  Called:
*     sla_CS2C, sla_VDV
*
*  P.T.Wallace   Starlink   9 March 1994
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

      REAL R2000,D2000

      REAL VA(3), VB(3)

      REAL sla_VDV

*
*  Peculiar solar motion from Delhaye 1965: in Galactic Cartesian
*  coordinates (+9,+12,+7) km/s.  This corresponds to about 16.6 km/s
*  towards Galactic coordinates L2 = 53 deg, B2 = +25 deg, or RA,Dec
*  17 49 58.7 +28 07 04 J2000.
*
*  The solar motion is expressed here in the form of a J2000.0
*  equatorial Cartesian vector:
*
*      VA(1) = X = -SPEED*COS(RA)*COS(DEC)
*      VA(2) = Y = -SPEED*SIN(RA)*COS(DEC)
*      VA(3) = Z = -SPEED*SIN(DEC)

      DATA VA / +0.63823, +14.58542, -7.80116 /



*  Convert given J2000 RA,Dec to x,y,z
      CALL sla_CS2C(R2000,D2000,VB)

*  Compute dot product with solar motion vector
      sla_RVLSRD=sla_VDV(VA,VB)

      END
