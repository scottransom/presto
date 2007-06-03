      SUBROUTINE sla_UE2EL (U, JFORMR,
     :                      JFORM, EPOCH, ORBINC, ANODE, PERIH,
     :                      AORQ, E, AORL, DM, JSTAT)
*+
*     - - - - - -
*      U E 2 E L
*     - - - - - -
*
*  Transform universal elements into conventional heliocentric
*  osculating elements.
*
*  Given:
*     U         d(13)  universal orbital elements (Note 1)
*
*                 (1)  combined mass (M+m)
*                 (2)  total energy of the orbit (alpha)
*                 (3)  reference (osculating) epoch (t0)
*               (4-6)  position at reference epoch (r0)
*               (7-9)  velocity at reference epoch (v0)
*                (10)  heliocentric distance at reference epoch
*                (11)  r0.v0
*                (12)  date (t)
*                (13)  universal eccentric anomaly (psi) of date, approx
*
*     JFORMR    i      requested element set (1-3; Note 3)
*
*  Returned:
*     JFORM     d      element set actually returned (1-3; Note 4)
*     EPOCH     d      epoch of elements (TT MJD)
*     ORBINC    d      inclination (radians)
*     ANODE     d      longitude of the ascending node (radians)
*     PERIH     d      longitude or argument of perihelion (radians)
*     AORQ      d      mean distance or perihelion distance (AU)
*     E         d      eccentricity
*     AORL      d      mean anomaly or longitude (radians, JFORM=1,2 only)
*     DM        d      daily motion (radians, JFORM=1 only)
*     JSTAT     i      status:  0 = OK
*                              -1 = illegal combined mass
*                              -2 = illegal JFORMR
*                              -3 = position/velocity out of range
*
*  Notes
*
*  1  The "universal" elements are those which define the orbit for the
*     purposes of the method of universal variables (see reference 2).
*     They consist of the combined mass of the two bodies, an epoch,
*     and the position and velocity vectors (arbitrary reference frame)
*     at that epoch.  The parameter set used here includes also various
*     quantities that can, in fact, be derived from the other
*     information.  This approach is taken to avoiding unnecessary
*     computation and loss of accuracy.  The supplementary quantities
*     are (i) alpha, which is proportional to the total energy of the
*     orbit, (ii) the heliocentric distance at epoch, (iii) the
*     outwards component of the velocity at the given epoch, (iv) an
*     estimate of psi, the "universal eccentric anomaly" at a given
*     date and (v) that date.
*
*  2  The universal elements are with respect to the mean equator and
*     equinox of epoch J2000.  The orbital elements produced are with
*     respect to the J2000 ecliptic and mean equinox.
*
*  3  Three different element-format options are supported:
*
*     Option JFORM=1, suitable for the major planets:
*
*     EPOCH  = epoch of elements (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = longitude of perihelion, curly pi (radians)
*     AORQ   = mean distance, a (AU)
*     E      = eccentricity, e
*     AORL   = mean longitude L (radians)
*     DM     = daily motion (radians)
*
*     Option JFORM=2, suitable for minor planets:
*
*     EPOCH  = epoch of elements (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = argument of perihelion, little omega (radians)
*     AORQ   = mean distance, a (AU)
*     E      = eccentricity, e
*     AORL   = mean anomaly M (radians)
*
*     Option JFORM=3, suitable for comets:
*
*     EPOCH  = epoch of perihelion (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = argument of perihelion, little omega (radians)
*     AORQ   = perihelion distance, q (AU)
*     E      = eccentricity, e
*
*  4  It may not be possible to generate elements in the form
*     requested through JFORMR.  The caller is notified of the form
*     of elements actually returned by means of the JFORM argument:
*
*      JFORMR   JFORM     meaning
*
*        1        1       OK - elements are in the requested format
*        1        2       never happens
*        1        3       orbit not elliptical
*
*        2        1       never happens
*        2        2       OK - elements are in the requested format
*        2        3       orbit not elliptical
*
*        3        1       never happens
*        3        2       never happens
*        3        3       OK - elements are in the requested format
*
*  5  The arguments returned for each value of JFORM (cf Note 6: JFORM
*     may not be the same as JFORMR) are as follows:
*
*         JFORM         1              2              3
*         EPOCH         t0             t0             T
*         ORBINC        i              i              i
*         ANODE         Omega          Omega          Omega
*         PERIH         curly pi       omega          omega
*         AORQ          a              a              q
*         E             e              e              e
*         AORL          L              M              -
*         DM            n              -              -
*
*     where:
*
*         t0           is the epoch of the elements (MJD, TT)
*         T              "    epoch of perihelion (MJD, TT)
*         i              "    inclination (radians)
*         Omega          "    longitude of the ascending node (radians)
*         curly pi       "    longitude of perihelion (radians)
*         omega          "    argument of perihelion (radians)
*         a              "    mean distance (AU)
*         q              "    perihelion distance (AU)
*         e              "    eccentricity
*         L              "    longitude (radians, 0-2pi)
*         M              "    mean anomaly (radians, 0-2pi)
*         n              "    daily motion (radians)
*         -             means no value is set
*
*  6  At very small inclinations, the longitude of the ascending node
*     ANODE becomes indeterminate and under some circumstances may be
*     set arbitrarily to zero.  Similarly, if the orbit is close to
*     circular, the true anomaly becomes indeterminate and under some
*     circumstances may be set arbitrarily to zero.  In such cases,
*     the other elements are automatically adjusted to compensate,
*     and so the elements remain a valid description of the orbit.
*
*  References:
*
*     1  Sterne, Theodore E., "An Introduction to Celestial Mechanics",
*        Interscience Publishers Inc., 1960.  Section 6.7, p199.
*
*     2  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
*
*  Called:  sla_PV2EL
*
*  P.T.Wallace   Starlink   18 March 1999
*
*  Copyright (C) 1999 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION U(13)
      INTEGER JFORMR,JFORM
      DOUBLE PRECISION EPOCH,ORBINC,ANODE,PERIH,AORQ,E,AORL,DM
      INTEGER JSTAT

*  Gaussian gravitational constant (exact)
      DOUBLE PRECISION GCON
      PARAMETER (GCON=0.01720209895D0)

*  Canonical days to seconds
      DOUBLE PRECISION CD2S
      PARAMETER (CD2S=GCON/86400D0)

      INTEGER I
      DOUBLE PRECISION PMASS,DATE,PV(6)


*  Unpack the universal elements.
      PMASS = U(1)-1D0
      DATE = U(3)
      DO I=1,3
         PV(I) = U(I+3)
         PV(I+3) = U(I+6)*CD2S
      END DO

*  Convert the position and velocity etc into conventional elements.
      CALL sla_PV2EL(PV,DATE,PMASS,JFORMR,JFORM,EPOCH,ORBINC,ANODE,
     :               PERIH,AORQ,E,AORL,DM,JSTAT)

      END
