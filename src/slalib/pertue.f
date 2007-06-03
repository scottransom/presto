      SUBROUTINE sla_PERTUE (DATE, U, JSTAT)
*+
*     - - - - - - -
*      P E R T U E
*     - - - - - - -
*
*  Update the universal elements of an asteroid or comet by applying
*  planetary perturbations.
*
*  Given:
*     DATE     d       final epoch (TT MJD) for the updated elements
*
*  Given and returned:
*     U        d(13)   universal elements (updated in place)
*
*                (1)   combined mass (M+m)
*                (2)   total energy of the orbit (alpha)
*                (3)   reference (osculating) epoch (t0)
*              (4-6)   position at reference epoch (r0)
*              (7-9)   velocity at reference epoch (v0)
*               (10)   heliocentric distance at reference epoch
*               (11)   r0.v0
*               (12)   date (t)
*               (13)   universal eccentric anomaly (psi) of date, approx
*
*  Returned:
*     JSTAT    i       status:
*                          +102 = warning, distant epoch
*                          +101 = warning, large timespan ( > 100 years)
*                     +1 to +10 = coincident with major planet (Note 5)
*                             0 = OK
*                            -1 = numerical error
*
*  Called:  sla_EPJ, sla_PLANET, sla_PV2UE, sla_UE2PV, sla_EPV,
*           sla_PREC, sla_DMOON, sla_DMXV
*
*  Notes:
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
*  2  The universal elements are with respect to the J2000 equator and
*     equinox.
*
*  3  The epochs DATE, U(3) and U(12) are all Modified Julian Dates
*     (JD-2400000.5).
*
*  4  The algorithm is a simplified form of Encke's method.  It takes as
*     a basis the unperturbed motion of the body, and numerically
*     integrates the perturbing accelerations from the major planets.
*     The expression used is essentially Sterne's 6.7-2 (reference 1).
*     Everhart and Pitkin (reference 2) suggest rectifying the orbit at
*     each integration step by propagating the new perturbed position
*     and velocity as the new universal variables.  In the present
*     routine the orbit is rectified less frequently than this, in order
*     to gain a slight speed advantage.  However, the rectification is
*     done directly in terms of position and velocity, as suggested by
*     Everhart and Pitkin, bypassing the use of conventional orbital
*     elements.
*
*     The f(q) part of the full Encke method is not used.  The purpose
*     of this part is to avoid subtracting two nearly equal quantities
*     when calculating the "indirect member", which takes account of the
*     small change in the Sun's attraction due to the slightly displaced
*     position of the perturbed body.  A simpler, direct calculation in
*     double precision proves to be faster and not significantly less
*     accurate.
*
*     Apart from employing a variable timestep, and occasionally
*     "rectifying the orbit" to keep the indirect member small, the
*     integration is done in a fairly straightforward way.  The
*     acceleration estimated for the middle of the timestep is assumed
*     to apply throughout that timestep;  it is also used in the
*     extrapolation of the perturbations to the middle of the next
*     timestep, to predict the new disturbed position.  There is no
*     iteration within a timestep.
*
*     Measures are taken to reach a compromise between execution time
*     and accuracy.  The starting-point is the goal of achieving
*     arcsecond accuracy for ordinary minor planets over a ten-year
*     timespan.  This goal dictates how large the timesteps can be,
*     which in turn dictates how frequently the unperturbed motion has
*     to be recalculated from the osculating elements.
*
*     Within predetermined limits, the timestep for the numerical
*     integration is varied in length in inverse proportion to the
*     magnitude of the net acceleration on the body from the major
*     planets.
*
*     The numerical integration requires estimates of the major-planet
*     motions.  Approximate positions for the major planets (Pluto
*     alone is omitted) are obtained from the routine sla_PLANET.  Two
*     levels of interpolation are used, to enhance speed without
*     significantly degrading accuracy.  At a low frequency, the routine
*     sla_PLANET is called to generate updated position+velocity "state
*     vectors".  The only task remaining to be carried out at the full
*     frequency (i.e. at each integration step) is to use the state
*     vectors to extrapolate the planetary positions.  In place of a
*     strictly linear extrapolation, some allowance is made for the
*     curvature of the orbit by scaling back the radius vector as the
*     linear extrapolation goes off at a tangent.
*
*     Various other approximations are made.  For example, perturbations
*     by Pluto and the minor planets are neglected and relativistic
*     effects are not taken into account.
*
*     In the interests of simplicity, the background calculations for
*     the major planets are carried out en masse.  The mean elements and
*     state vectors for all the planets are refreshed at the same time,
*     without regard for orbit curvature, mass or proximity.
*
*     The Earth-Moon system is treated as a single body when the body is
*     distant but as separate bodies when closer to the EMB than the
*     parameter RNE, which incurs a time penalty but improves accuracy
*     for near-Earth objects.
*
*  5  This routine is not intended to be used for major planets.
*     However, if major-planet elements are supplied, sensible results
*     will, in fact, be produced.  This happens because the routine
*     checks the separation between the body and each of the planets and
*     interprets a suspiciously small value (0.001 AU) as an attempt to
*     apply the routine to the planet concerned.  If this condition is
*     detected, the contribution from that planet is ignored, and the
*     status is set to the planet number (1-10 = Mercury, Venus, EMB,
*     Mars, Jupiter, Saturn, Uranus, Neptune, Earth, Moon) as a warning.
*
*  References:
*
*     1  Sterne, Theodore E., "An Introduction to Celestial Mechanics",
*        Interscience Publishers Inc., 1960.  Section 6.7, p199.
*
*     2  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
*
*  Last revision:   27 December 2004
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
      DOUBLE PRECISION DATE,U(13)
      INTEGER JSTAT

*  Distance from EMB at which Earth and Moon are treated separately
      DOUBLE PRECISION RNE
      PARAMETER (RNE=1D0)

*  Coincidence with major planet distance
      DOUBLE PRECISION COINC
      PARAMETER (COINC=0.0001D0)

*  Coefficient relating timestep to perturbing force
      DOUBLE PRECISION TSC
      PARAMETER (TSC=1D-4)

*  Minimum and maximum timestep (days)
      DOUBLE PRECISION TSMIN,TSMAX
      PARAMETER (TSMIN=0.01D0,TSMAX=10D0)

*  Age limit for major-planet state vector (days)
      DOUBLE PRECISION AGEPMO
      PARAMETER (AGEPMO=5D0)

*  Age limit for major-planet mean elements (days)
      DOUBLE PRECISION AGEPEL
      PARAMETER (AGEPEL=50D0)

*  Margin for error when deciding whether to renew the planetary data
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1D-6)

*  Age limit for the body's osculating elements (before rectification)
      DOUBLE PRECISION AGEBEL
      PARAMETER (AGEBEL=100D0)

*  Gaussian gravitational constant (exact) and square
      DOUBLE PRECISION GCON,GCON2
      PARAMETER (GCON=0.01720209895D0,GCON2=GCON*GCON)

*  The final epoch
      DOUBLE PRECISION TFINAL

*  The body's current universal elements
      DOUBLE PRECISION UL(13)

*  Current reference epoch
      DOUBLE PRECISION T0

*  Timespan from latest orbit rectification to final epoch (days)
      DOUBLE PRECISION TSPAN

*  Time left to go before integration is complete
      DOUBLE PRECISION TLEFT

*  Time direction flag: +1=forwards, -1=backwards
      DOUBLE PRECISION FB

*  First-time flag
      LOGICAL FIRST

*
*  The current perturbations
*
*  Epoch (days relative to current reference epoch)
      DOUBLE PRECISION RTN
*  Position (AU)
      DOUBLE PRECISION PERP(3)
*  Velocity (AU/d)
      DOUBLE PRECISION PERV(3)
*  Acceleration (AU/d/d)
      DOUBLE PRECISION PERA(3)
*

*  Length of current timestep (days), and half that
      DOUBLE PRECISION TS,HTS

*  Epoch of middle of timestep
      DOUBLE PRECISION T

*  Epoch of planetary mean elements
      DOUBLE PRECISION TPEL

*  Planet number (1=Mercury, 2=Venus, 3=EMB...8=Neptune)
      INTEGER NP

*  Planetary universal orbital elements
      DOUBLE PRECISION UP(13,8)

*  Epoch of planetary state vectors
      DOUBLE PRECISION TPMO

*  State vectors for the major planets (AU,AU/s)
      DOUBLE PRECISION PVIN(6,8)

*  Earth velocity and position vectors (AU,AU/s)
      DOUBLE PRECISION VB(3),PB(3),VH(3),PE(3)

*  Moon geocentric state vector (AU,AU/s) and position part
      DOUBLE PRECISION PVM(6),PM(3)

*  Date to J2000 de-precession matrix
      DOUBLE PRECISION PMAT(3,3)

*
*  Correction terms for extrapolated major planet vectors
*
*  Sun-to-planet distances squared multiplied by 3
      DOUBLE PRECISION R2X3(8)
*  Sunward acceleration terms, G/2R^3
      DOUBLE PRECISION GC(8)
*  Tangential-to-circular correction factor
      DOUBLE PRECISION FC
*  Radial correction factor due to Sunwards acceleration
      DOUBLE PRECISION FG
*

*  The body's unperturbed and perturbed state vectors (AU,AU/s)
      DOUBLE PRECISION PV0(6),PV(6)

*  The body's perturbed and unperturbed heliocentric distances (AU) cubed
      DOUBLE PRECISION R03,R3

*  The perturbating accelerations, indirect and direct
      DOUBLE PRECISION FI(3),FD(3)

*  Sun-to-planet vector, and distance cubed
      DOUBLE PRECISION RHO(3),RHO3

*  Body-to-planet vector, and distance cubed
      DOUBLE PRECISION DELTA(3),DELTA3

*  Miscellaneous
      INTEGER I,J
      DOUBLE PRECISION R2,W,DT,DT2,R,FT
      LOGICAL NE

      DOUBLE PRECISION sla_EPJ

*  Planetary inverse masses, Mercury through Neptune then Earth and Moon
      DOUBLE PRECISION AMAS(10)
      DATA AMAS / 6023600D0, 408523.5D0, 328900.5D0, 3098710D0,
     :            1047.355D0, 3498.5D0, 22869D0, 19314D0,
     :            332946.038D0, 27068709D0 /

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
*----------------------------------------------------------------------*


*  Preset the status to OK.
      JSTAT = 0

*  Copy the final epoch.
      TFINAL = DATE

*  Copy the elements (which will be periodically updated).
      DO I=1,13
         UL(I) = U(I)
      END DO

*  Initialize the working reference epoch.
      T0=UL(3)

*  Total timespan (days) and hence time left.
      TSPAN = TFINAL-T0
      TLEFT = TSPAN

*  Warn if excessive.
      IF (ABS(TSPAN).GT.36525D0) JSTAT=101

*  Time direction: +1 for forwards, -1 for backwards.
      FB = SIGN(1D0,TSPAN)

*  Initialize relative epoch for start of current timestep.
      RTN = 0D0

*  Reset the perturbations (position, velocity, acceleration).
      DO I=1,3
         PERP(I) = 0D0
         PERV(I) = 0D0
         PERA(I) = 0D0
      END DO

*  Set "first iteration" flag.
      FIRST = .TRUE.

*  Step through the time left.
      DO WHILE (FB*TLEFT.GT.0D0)

*     Magnitude of current acceleration due to planetary attractions.
         IF (FIRST) THEN
            TS = TSMIN
         ELSE
            R2 = 0D0
            DO I=1,3
               W = FD(I)
               R2 = R2+W*W
            END DO
            W = SQRT(R2)

*        Use the acceleration to decide how big a timestep can be tolerated.
            IF (W.NE.0D0) THEN
               TS = MIN(TSMAX,MAX(TSMIN,TSC/W))
            ELSE
               TS = TSMAX
            END IF
         END IF
         TS = TS*FB

*     Override if final epoch is imminent.
         TLEFT = TSPAN-RTN
         IF (ABS(TS).GT.ABS(TLEFT)) TS=TLEFT

*     Epoch of middle of timestep.
         HTS = TS/2D0
         T = T0+RTN+HTS

*     Is it time to recompute the major-planet elements?
         IF (FIRST.OR.ABS(T-TPEL)-AGEPEL.GE.TINY) THEN

*        Yes: go forward in time by just under the maximum allowed.
            TPEL = T+FB*AGEPEL

*        Compute the state vector for the new epoch.
            DO NP=1,8
               CALL sla_PLANET(TPEL,NP,PV,J)

*           Warning if remote epoch, abort if error.
               IF (J.EQ.1) THEN
                  JSTAT = 102
               ELSE IF (J.NE.0) THEN
                  GO TO 9010
               END IF

*           Transform the vector into universal elements.
               CALL sla_PV2UE(PV,TPEL,0D0,UP(1,NP),J)
               IF (J.NE.0) GO TO 9010
            END DO
         END IF

*     Is it time to recompute the major-planet motions?
         IF (FIRST.OR.ABS(T-TPMO)-AGEPMO.GE.TINY) THEN

*        Yes: look ahead.
            TPMO = T+FB*AGEPMO

*        Compute the motions of each planet (AU,AU/d).
            DO NP=1,8

*           The planet's position and velocity (AU,AU/s).
               CALL sla_UE2PV(TPMO,UP(1,NP),PVIN(1,NP),J)
               IF (J.NE.0) GO TO 9010

*           Scale velocity to AU/d.
               DO J=4,6
                  PVIN(J,NP) = PVIN(J,NP)*86400D0
               END DO

*           Precompute also the extrapolation correction terms.
               R2 = 0D0
               DO I=1,3
                  W = PVIN(I,NP)
                  R2 = R2+W*W
               END DO
               R2X3(NP) = R2*3D0
               GC(NP) = GCON2/(2D0*R2*SQRT(R2))
            END DO
         END IF

*     Reset the first-time flag.
         FIRST = .FALSE.

*     Unperturbed motion of the body at middle of timestep (AU,AU/s).
         CALL sla_UE2PV(T,UL,PV0,J)
         IF (J.NE.0) GO TO 9010

*     Perturbed position of the body (AU) and heliocentric distance cubed.
         R2 = 0D0
         DO I=1,3
            W = PV0(I)+PERP(I)+(PERV(I)+PERA(I)*HTS/2D0)*HTS
            PV(I) = W
            R2 = R2+W*W
         END DO
         R3 = R2*SQRT(R2)

*     The body's unperturbed heliocentric distance cubed.
         R2 = 0D0
         DO I=1,3
            W = PV0(I)
            R2 = R2+W*W
         END DO
         R03 = R2*SQRT(R2)

*     Compute indirect and initialize direct parts of the perturbation.
         DO I=1,3
            FI(I) = PV0(I)/R03-PV(I)/R3
            FD(I) = 0D0
         END DO

*     Ready to compute the direct planetary effects.

*     Reset the "near-Earth" flag.
         NE = .FALSE.

*     Interval from state-vector epoch to middle of current timestep.
         DT = T-TPMO
         DT2 = DT*DT

*     Planet by planet, including separate Earth and Moon.
         DO NP=1,10

*        Which perturbing body?
            IF (NP.LE.8) THEN

*           Planet: compute the extrapolation in longitude (squared).
               R2 = 0D0
               DO J=4,6
                  W = PVIN(J,NP)*DT
                  R2 = R2+W*W
               END DO

*           Hence the tangential-to-circular correction factor.
               FC = 1D0+R2/R2X3(NP)

*           The radial correction factor due to the inwards acceleration.
               FG = 1D0-GC(NP)*DT2

*           Planet's position.
               DO I=1,3
                  RHO(I) = FG*(PVIN(I,NP)+FC*PVIN(I+3,NP)*DT)
               END DO

            ELSE IF (NE) THEN

*           Near-Earth and either Earth or Moon.

               IF (NP.EQ.9) THEN

*              Earth: position.
                  CALL sla_EPV(T,PE,VH,PB,VB)
                  DO I=1,3
                     RHO(I) = PE(I)
                  END DO

               ELSE

*              Moon: position.
                  CALL sla_PREC(sla_EPJ(T),2000D0,PMAT)
                  CALL sla_DMOON(T,PVM)
                  CALL sla_DMXV(PMAT,PVM,PM)
                  DO I=1,3
                     RHO(I) = PM(I)+PE(I)
                  END DO
               END IF
            END IF

*        Proceed unless Earth or Moon and not the near-Earth case.
            IF (NP.LE.8.OR.NE) THEN

*           Heliocentric distance cubed.
               R2 = 0D0
               DO I=1,3
                  W = RHO(I)
                  R2 = R2+W*W
               END DO
               R = SQRT(R2)
               RHO3 = R2*R

*           Body-to-planet vector, and distance.
               R2 = 0D0
               DO I=1,3
                  W = RHO(I)-PV(I)
                  DELTA(I) = W
                  R2 = R2+W*W
               END DO
               R = SQRT(R2)

*           If this is the EMB, set the near-Earth flag appropriately.
               IF (NP.EQ.3.AND.R.LT.RNE) NE = .TRUE.

*           Proceed unless EMB and this is the near-Earth case.
               IF (.NOT.(NE.AND.NP.EQ.3)) THEN

*              If too close, ignore this planet and set a warning.
                  IF (R.LT.COINC) THEN
                     JSTAT = NP

                  ELSE

*                 Accumulate "direct" part of perturbation acceleration.
                     DELTA3 = R2*R
                     W = AMAS(NP)
                     DO I=1,3
                        FD(I) = FD(I)+(DELTA(I)/DELTA3-RHO(I)/RHO3)/W
                     END DO
                  END IF
               END IF
            END IF
         END DO

*     Update the perturbations to the end of the timestep.
         RTN = RTN+TS
         DO I=1,3
            W = (FI(I)+FD(I))*GCON2
            FT = W*TS
            PERP(I) = PERP(I)+(PERV(I)+FT/2D0)*TS
            PERV(I) = PERV(I)+FT
            PERA(I) = W
         END DO

*     Time still to go.
         TLEFT = TSPAN-RTN

*     Is it either time to rectify the orbit or the last time through?
         IF (ABS(RTN).GE.AGEBEL.OR.FB*TLEFT.LE.0D0) THEN

*        Yes: update to the end of the current timestep.
            T0 = T0+RTN
            RTN = 0D0

*        The body's unperturbed motion (AU,AU/s).
            CALL sla_UE2PV(T0,UL,PV0,J)
            IF (J.NE.0) GO TO 9010

*        Add and re-initialize the perturbations.
            DO I=1,3
               J = I+3
               PV(I) = PV0(I)+PERP(I)
               PV(J) = PV0(J)+PERV(I)/86400D0
               PERP(I) = 0D0
               PERV(I) = 0D0
               PERA(I) = FD(I)*GCON2
            END DO

*        Use the position and velocity to set up new universal elements.
            CALL sla_PV2UE(PV,T0,0D0,UL,J)
            IF (J.NE.0) GO TO 9010

*        Adjust the timespan and time left.
            TSPAN = TFINAL-T0
            TLEFT = TSPAN
         END IF

*     Next timestep.
      END DO

*  Return the updated universal-element set.
      DO I=1,13
         U(I) = UL(I)
      END DO

*  Finished.
      GO TO 9999

*  Miscellaneous numerical error.
 9010 CONTINUE
      JSTAT = -1

 9999 CONTINUE
      END
