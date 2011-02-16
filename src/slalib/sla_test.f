      PROGRAM SLA_TEST
*+
*  - - - - - - - - -
*   S L A _ T E S T
*  - - - - - - - - -
*
*  Validate the slalib library.
*
*  Each slalib function is tested to some useful but in most cases
*  not exhaustive level.  Successful completion is signalled by an
*  absence of output messages.  Failure of a given function or
*  group of functions results in error messages.
*
*  Any messages go to standard output.
*
*  Adapted from original C code by P.T.Wallace.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink and P.T.Wallace.  All rights reserved.
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

      LOGICAL STATUS
      INTEGER EXITSTATUS


*  Preset the status to success.
      STATUS = .TRUE.

*  Test all the slalib functions.
      CALL T_ADDET ( STATUS )
      CALL T_AFIN ( STATUS )
      CALL T_AIRMAS ( STATUS )
      CALL T_ALTAZ ( STATUS )
      CALL T_AMP ( STATUS )
      CALL T_AOP ( STATUS )
      CALL T_BEAR ( STATUS )
      CALL T_CAF2R ( STATUS )
      CALL T_CALDJ ( STATUS )
      CALL T_CALYD ( STATUS )
      CALL T_CC2S ( STATUS )
      CALL T_CC62S ( STATUS )
      CALL T_CD2TF ( STATUS )
      CALL T_CLDJ ( STATUS )
      CALL T_CR2AF ( STATUS )
      CALL T_CR2TF ( STATUS )
      CALL T_CS2C6 ( STATUS )
      CALL T_CTF2D ( STATUS )
      CALL T_CTF2R ( STATUS )
      CALL T_DAT ( STATUS )
      CALL T_DBJIN ( STATUS )
      CALL T_DJCAL ( STATUS )
      CALL T_DMAT ( STATUS )
      CALL T_E2H ( STATUS )
      CALL T_EARTH ( STATUS )
      CALL T_ECLEQ ( STATUS )
      CALL T_ECMAT ( STATUS )
      CALL T_ECOR ( STATUS )
      CALL T_EG50 ( STATUS )
      CALL T_EPB ( STATUS )
      CALL T_EPB2D ( STATUS )
      CALL T_EPCO ( STATUS )
      CALL T_EPJ ( STATUS )
      CALL T_EPJ2D ( STATUS )
      CALL T_EQECL ( STATUS )
      CALL T_EQEQX ( STATUS )
      CALL T_EQGAL ( STATUS )
      CALL T_ETRMS ( STATUS )
      CALL T_EVP ( STATUS )
      CALL T_FITXY ( STATUS )
      CALL T_FK425 ( STATUS )
      CALL T_FK45Z ( STATUS )
      CALL T_FK524 ( STATUS )
      CALL T_FK52H ( STATUS )
      CALL T_FK54Z ( STATUS )
      CALL T_FLOTIN ( STATUS )
      CALL T_GALEQ ( STATUS )
      CALL T_GALSUP ( STATUS )
      CALL T_GE50 ( STATUS )
      CALL T_GMST ( STATUS )
      CALL T_INTIN ( STATUS )
      CALL T_KBJ ( STATUS )
      CALL T_MAP ( STATUS )
      CALL T_MOON ( STATUS )
      CALL T_NUT ( STATUS )
      CALL T_OBS ( STATUS )
      CALL T_PA ( STATUS )
      CALL T_PCD ( STATUS )
      CALL T_PDA2H ( STATUS )
      CALL T_PDQ2H ( STATUS )
      CALL T_PERCOM ( STATUS )
      CALL T_PLANET ( STATUS )
      CALL T_PM ( STATUS )
      CALL T_POLMO ( STATUS )
      CALL T_PREBN ( STATUS )
      CALL T_PREC ( STATUS )
      CALL T_PRECES ( STATUS )
      CALL T_PRENUT ( STATUS )
      CALL T_PVOBS ( STATUS )
      CALL T_RANGE ( STATUS )
      CALL T_RANORM ( STATUS )
      CALL T_RCC ( STATUS )
      CALL T_REF ( STATUS )
      CALL T_RV ( STATUS )
      CALL T_SEP ( STATUS )
      CALL T_SMAT ( STATUS )
      CALL T_SUPGAL ( STATUS )
      CALL T_SVD ( STATUS )
      CALL T_TP ( STATUS )
      CALL T_TPV ( STATUS )
      CALL T_VECMAT ( STATUS )
      CALL T_ZD ( STATUS )

*  Report any errors and set up an appropriate exit status.  Set the
*  EXITSTATUS to 0 on success, 1 on any error -- Unix-style.  The
*  EXIT intrinsic is non-standard but common (which is portable enough
*  for a regression test).

      IF ( STATUS ) THEN
         WRITE (*,'(1X,''SLALIB validation OK!'')')
         EXITSTATUS = 0
      ELSE
         WRITE (*,'(1X,''SLALIB validation failed!'')')
         EXITSTATUS = 1
      ENDIF

      CALL EXIT(EXITSTATUS)

      END

      SUBROUTINE VCS ( S, SOK, FUNC, TEST, STATUS )
*+
*  - - - -
*   V C S
*  - - - -
*
*  Validate a character string result.
*
*  Internal routine used by sla_TEST program.
*
*  Given:
*     S        CHARACTER    string produced by routine under test
*     SOK      CHARACTER    correct value
*     FUNC     CHARACTER    name of routine under test
*     TEST     CHARACTER    name of individual test (or ' ')
*
*  Given and returned:
*     STATUS   LOGICAL      set to .FALSE. if test fails
*
*  Called:  ERR
*
*  Last revision:   25 May 2002
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      CHARACTER*(*) S, SOK, FUNC, TEST
      LOGICAL STATUS


      IF ( S .NE. SOK ) THEN
         CALL ERR ( FUNC, TEST, STATUS )
         WRITE (*,'(1X,''  expected ='',6X,''"'',A,''"'')') SOK
         WRITE (*,'(1X,''  actual =  '',6X,''"'',A,''"'')') S
      END IF

      END

      SUBROUTINE VIV ( IVAL, IVALOK, FUNC, TEST, STATUS )
*+
*  - - - -
*   V I V
*  - - - -
*
*  Validate an integer result.
*
*  Internal routine used by sla_TEST program.
*
*  Given:
*     IVAL     INTEGER      value computed by routine under test
*     IVALOK   INTEGER      correct value
*     FUNC     CHARACTER    name of routine under test
*     TEST     CHARACTER    name of individual test (or ' ')
*
*  Given and returned:
*     STATUS   LOGICAL      set to .FALSE. if test fails
*
*  Called:  ERR
*
*  Last revision:   25 May 2002
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      INTEGER IVAL, IVALOK
      CHARACTER*(*) FUNC, TEST
      LOGICAL STATUS


      IF ( IVAL .NE. IVALOK ) THEN
         CALL ERR ( FUNC, TEST, STATUS )
         WRITE (*,'(1X,''  expected ='',I10)') IVALOK
         WRITE (*,'(1X,''  actual =  '',I10)') IVAL
      END IF

      END

      SUBROUTINE VLV ( IVAL, IVALOK, FUNC, TEST, STATUS )
*+
*  - - - -
*   V L V
*  - - - -
*
*  Validate a long result.
*
*  Internal routine used by sla_TEST program.
*
*  Given:
*     IVAL     INTEGER*4    value computed by routine under test
*     IVALOK   INTEGER*4    correct value
*     FUNC     CHARACTER    name of routine under test
*     TEST     CHARACTER    name of individual test (or ' ')
*
*  Given and returned:
*     STATUS   LOGICAL      set to .FALSE. if test fails
*
*  Called:  ERR
*
*  Last revision:   25 May 2002
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      INTEGER*4 IVAL, IVALOK
      CHARACTER*(*) FUNC, TEST
      LOGICAL STATUS


      IF ( IVAL .NE. IVALOK ) THEN
         CALL ERR ( FUNC, TEST, STATUS )
         WRITE (*,'(1X,''  expected ='',I10)') IVALOK
         WRITE (*,'(1X,''  actual =  '',I10)') IVAL
      END IF

      END

      SUBROUTINE VVD ( VAL, VALOK, DVAL, FUNC, TEST, STATUS )
*+
*  - - - -
*   V V D
*  - - - -
*
*  Validate a double result.
*
*  Internal routine used by sla_TEST program.
*
*  Given:
*     VAL      DOUBLE       value computed by routine under test
*     VALOK    DOUBLE       correct value
*     DVAL     DOUBLE       maximum allowable error
*     FUNC     CHARACTER    name of routine under test
*     TEST     CHARACTER    name of individual test (or ' ')
*
*  Given and returned:
*     STATUS   LOGICAL      set to .FALSE. if test fails
*
*  Called:  ERR
*
*  Last revision:   25 May 2002
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      DOUBLE PRECISION VAL, VALOK, DVAL
      CHARACTER*(*) FUNC, TEST
      LOGICAL STATUS


      IF ( DABS ( VAL - VALOK ) .GT. DVAL ) THEN
         CALL ERR ( FUNC, TEST, STATUS )
         WRITE (*,'(1X,''  expected ='',G30.19)') VALOK
         WRITE (*,'(1X,''  actual =  '',G30.19)') VAL
      END IF

      END

      SUBROUTINE ERR ( FUNC, TEST, STATUS )
*+
*  - - - -
*   E R R
*  - - - -
*
*  Report a failed test.
*
*  Internal routine used by sla_TEST program.
*
*  Given:
*     FUNC     CHARACTER    name of routine under test
*     TEST     CHARACTER    name of individual test (or ' ')
*
*  Given and returned:
*     STATUS   LOGICAL      set to .FALSE.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      CHARACTER*(*) FUNC, TEST
      LOGICAL STATUS


      WRITE (*,'(1X,A,'' test '',A,'' fails:'')') FUNC, TEST
      STATUS = .FALSE.

      END

      SUBROUTINE T_ADDET ( STATUS )
*+
*  - - - - - - - -
*   T _ A D D E T
*  - - - - - - - -
*
*  Test sla_ADDET, sla_SUBET routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_ADDET, VVD, sla_SUBET.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RM, DM, EQ, R1, D1, R2, D2

      RM = 2D0
      DM = -1D0
      EQ = 1975D0

      CALL sla_ADDET ( RM, DM, EQ, R1, D1 )
      CALL VVD ( R1 - RM, 2.983864874295250D-6, 1D-12, 'sla_ADDET',
     :           'R', STATUS )
      CALL VVD ( D1 - DM, 2.379650804185118D-7, 1D-12, 'sla_ADDET',
     :           'D', STATUS )

      CALL sla_SUBET ( R1, D1, EQ, R2, D2 )
      CALL VVD ( R2 - RM, 0D0, 1D-12, 'sla_SUBET', 'R', STATUS )
      CALL VVD ( D2 - DM, 0D0, 1D-12, 'sla_SUBET', 'D', STATUS )

      END

      SUBROUTINE T_AFIN ( STATUS )
*+
*  - - - - - - -
*   T _ A F I N
*  - - - - - - -
*
*  Test sla_AFIN and sla_DAFIN routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_AFIN, VIV, VVD, sla_DAFIN.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I, J
      REAL F
      DOUBLE PRECISION D
      CHARACTER*12 S
      DATA S /'12 34 56.7 |'/


      I = 1
      CALL sla_AFIN (S, I, F, J)
      CALL VIV ( I, 12, 'sla_AFIN', 'I', STATUS )
      CALL VVD ( DBLE( F ), 0.2196045986911432D0, 1D-6, 'sla_AFIN',
     :           'A', STATUS )
      CALL VIV ( J, 0, 'sla_AFIN', 'J', STATUS )

      I = 1
      CALL sla_DAFIN (S, I, D, J)
      CALL VIV ( I, 12, 'sla_DAFIN', 'I', STATUS )
      CALL VVD ( D, 0.2196045986911432D0, 1D-12, 'sla_DAFIN', 'A',
     :           STATUS )
      CALL VIV ( J, 0, 'sla_DAFIN', 'J', STATUS )

      END

      SUBROUTINE T_AIRMAS ( STATUS )
*+
*  - - - - - - - - -
*   T _ A I R M A S
*  - - - - - - - - -
*
*  Test sla_AIRMAS routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called: VVD, sla_AIRMAS.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_AIRMAS


      CALL VVD ( sla_AIRMAS ( 1.2354D0 ), 3.015698990074724D0,
     :           1D-12, 'sla_AIRMAS', ' ', STATUS )

      END

      SUBROUTINE T_ALTAZ ( STATUS )
*+
*  - - - - - - - -
*   T _ A L T A Z
*  - - - - - - - -
*
*  Test sla_ALTAZ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_ALTAZ, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION AZ, AZD, AZDD, EL, ELD, ELDD, PA, PAD, PADD

      CALL sla_ALTAZ ( 0.7D0, -0.7D0, -0.65D0,
     :                 AZ, AZD, AZDD, EL, ELD, ELDD, PA, PAD, PADD )

      CALL VVD ( AZ, 4.400560746660174D0, 1D-12, 'sla_ALTAZ',
     :           'AZ', STATUS )
      CALL VVD ( AZD, -0.2015438937145421D0, 1D-13, 'sla_ALTAZ',
     :           'AZD', STATUS )
      CALL VVD ( AZDD, -0.4381266949668748D0, 1D-13, 'sla_ALTAZ',
     :           'AZDD', STATUS )
      CALL VVD ( EL, 1.026646506651396D0, 1D-12, 'sla_ALTAZ',
     :           'EL', STATUS )
      CALL VVD ( ELD, -0.7576920683826450D0, 1D-13, 'sla_ALTAZ',
     :           'ELD', STATUS )
      CALL VVD ( ELDD, 0.04922465406857453D0, 1D-14, 'sla_ALTAZ',
     :           'ELDD', STATUS )
      CALL VVD ( PA, 1.707639969653937D0, 1D-12, 'sla_ALTAZ',
     :           'PA', STATUS )
      CALL VVD ( PAD, 0.4717832355365627D0, 1D-13, 'sla_ALTAZ',
     :           'PAD', STATUS )
      CALL VVD ( PADD, -0.2957914128185515D0, 1D-13, 'sla_ALTAZ',
     :           'PADD', STATUS )

      END

      SUBROUTINE T_AMP ( STATUS )
*+
*  - - - - - -
*   T _ A M P
*  - - - - - -
*
*  Test sla_AMP, sla_MAPPA, sla_AMPQK routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_AMP, VVD.
*
*  Last revision:   16 November 2001
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RM, DM

      CALL sla_AMP ( 2.345D0, -1.234D0, 50100D0, 1990D0, RM, DM )
      CALL VVD ( RM, 2.344472180027961D0, 1D-11, 'sla_AMP', 'R',
     :           STATUS )
      CALL VVD ( DM, -1.233573099847705D0, 1D-11, 'sla_AMP', 'D',
     :           STATUS )

      END

      SUBROUTINE T_AOP ( STATUS )
*+
*  - - - - - -
*   T _ A O P
*  - - - - - -
*
*  Test sla_AOP, sla_AOPPA, sla_AOPQK, sla_OAP, sla_OAPQK,
*  sla_AOPPAT routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_AOP, VVD, sla_AOPPA, sla_AOPQK, sla_OAP, sla_OAPQK,
*  sla_AOPPAT.
*
*  Defined in slamac.h:  DS2R
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I
      DOUBLE PRECISION DS2R
      DOUBLE PRECISION RAP, DAP, DATE, DUT, ELONGM, PHIM, HM, XP, YP,
     :    TDK, PMB, RH, WL, TLR, AOB, ZOB, HOB, DOB, ROB, AOPRMS(14)

      PARAMETER (DS2R =
     :      7.2722052166430399038487115353692196393452995355905D-5)

      DAP = -0.1234D0
      DATE = 51000.1D0
      DUT = 25D0
      ELONGM = 2.1D0
      PHIM = 0.5D0
      HM = 3000D0
      XP = -0.5D-6
      YP = 1D-6
      TDK = 280D0
      PMB = 550D0
      RH = 0.6D0
      TLR = 0.006D0

      DO I = 1, 3

         IF ( I .EQ. 1 ) THEN
            RAP = 2.7D0
            WL = 0.45D0
            ELSE IF ( I .EQ. 2 ) THEN
                RAP = 2.345D0
            ELSE
                WL = 1D6
         END IF

         CALL sla_AOP ( RAP, DAP, DATE, DUT, ELONGM, PHIM, HM, XP, YP,
     :                  TDK, PMB, RH, WL, TLR, AOB, ZOB, HOB, DOB, ROB )

         IF ( I .EQ. 1 ) THEN
            CALL VVD ( AOB, 1.812817787123283034D0, 1D-10, 'sla_AOP',
     :                 'lo aob', STATUS )
            CALL VVD ( ZOB, 1.393860816635714034D0, 1D-10, 'sla_AOP',
     :                 'lo zob', STATUS )
            CALL VVD ( HOB, -1.297808009092456683D0, 1D-10, 'sla_AOP',
     :                 'lo hob', STATUS )
            CALL VVD ( DOB, -0.122967060534561D0, 1D-10, 'sla_AOP',
     :                 'lo dob', STATUS )
            CALL VVD ( ROB, 2.699270287872084D0, 1D-10, 'sla_AOP',
     :                    'lo rob', STATUS )
         ELSE IF ( I .EQ. 2 ) THEN
            CALL VVD ( AOB, 2.019928026670621442D0, 1D-10, 'sla_AOP',
     :                 'aob/o', STATUS )
            CALL VVD ( ZOB, 1.101316172427482466D0, 1D-10, 'sla_AOP',
     :                 'zob/o', STATUS )
            CALL VVD ( HOB, -0.9432923558497740862D0, 1D-10, 'sla_AOP',
     :                 'hob/o', STATUS )
            CALL VVD ( DOB, -0.1232144708194224D0, 1D-10, 'sla_AOP',
     :                 'dob/o', STATUS )
            CALL VVD ( ROB, 2.344754634629428D0, 1D-10, 'sla_AOP',
     :                 'rob/o', STATUS )
         ELSE
            CALL VVD ( AOB, 2.019928026670621442D0, 1D-10, 'sla_AOP',
     :                 'aob/r', STATUS )
            CALL VVD ( ZOB, 1.101267532198003760D0, 1D-10, 'sla_AOP',
     :                 'zob/r', STATUS )
            CALL VVD ( HOB, -0.9432533138143315937D0, 1D-10, 'sla_AOP',
     :                 'hob/r', STATUS )
            CALL VVD ( DOB, -0.1231850665614878D0, 1D-10, 'sla_AOP',
     :                 'dob/r', STATUS )
            CALL VVD ( ROB, 2.344715592593984D0, 1D-10, 'sla_AOP',
     :                 'rob/r', STATUS )
         END IF
      END DO

      DATE = 48000.3D0
      WL = 0.45D0

      CALL sla_AOPPA ( DATE, DUT, ELONGM, PHIM, HM, XP, YP, TDK,
     :                 PMB, RH, WL, TLR, AOPRMS )
      CALL VVD ( AOPRMS(1), 0.4999993892136306D0, 1D-13, 'sla_AOPPA',
     :           '1', STATUS )
      CALL VVD ( AOPRMS(2), 0.4794250025886467D0, 1D-13, 'sla_AOPPA',
     :           '2', STATUS )
      CALL VVD ( AOPRMS(3), 0.8775828547167932D0, 1D-13, 'sla_AOPPA',
     :           '3', STATUS )
      CALL VVD ( AOPRMS(4), 1.363180872136126D-6, 1D-13, 'sla_AOPPA',
     :           '4', STATUS )
      CALL VVD ( AOPRMS(5), 3000D0, 1D-10, 'sla_AOPPA', '5',
     :           STATUS )
      CALL VVD ( AOPRMS(6), 280D0, 1D-11, 'sla_AOPPA', '6',
     :           STATUS )
      CALL VVD ( AOPRMS(7), 550D0, 1D-11, 'sla_AOPPA', '7',
     :           STATUS )
      CALL VVD ( AOPRMS(8), 0.6D0, 1D-13, 'sla_AOPPA', '8',
     :           STATUS )
      CALL VVD ( AOPRMS(9), 0.45D0, 1D-13, 'sla_AOPPA', '9',
     :           STATUS )
      CALL VVD ( AOPRMS(10), 0.006D0, 1D-15, 'sla_AOPPA', '10',
     :           STATUS )
      CALL VVD ( AOPRMS(11), 0.0001562803328459898D0, 1D-13,
     :           'sla_AOPPA', '11', STATUS )
      CALL VVD ( AOPRMS(12), -1.792293660141D-7, 1D-13,
     :           'sla_AOPPA', '12', STATUS )
      CALL VVD ( AOPRMS(13), 2.101874231495843D0, 1D-13,
     :           'sla_AOPPA', '13', STATUS )
      CALL VVD ( AOPRMS(14), 7.601916802079765D0, 1D-8,
     :           'sla_AOPPA', '14', STATUS )

      CALL sla_OAP ( 'R', 1.6D0, -1.01D0, DATE, DUT, ELONGM, PHIM,
     :               HM, XP, YP, TDK, PMB, RH, WL, TLR, RAP, DAP )
      CALL VVD ( RAP, 1.601197569844787D0, 1D-10, 'sla_OAP',
     :           'Rr', STATUS )
      CALL VVD ( DAP, -1.012528566544262D0, 1D-10, 'sla_OAP',
     :           'Rd', STATUS )
      CALL sla_OAP ( 'H', -1.234D0, 2.34D0, DATE, DUT, ELONGM, PHIM,
     :               HM, XP, YP, TDK, PMB, RH, WL, TLR, RAP, DAP )
      CALL VVD ( RAP, 5.693087688154886463D0, 1D-10, 'sla_OAP',
     :           'Hr', STATUS )
      CALL VVD ( DAP, 0.8010281167405444D0, 1D-10, 'sla_OAP',
     :           'Hd', STATUS )
      CALL sla_OAP ( 'A', 6.1D0, 1.1D0, DATE, DUT, ELONGM, PHIM,
     :               HM, XP, YP, TDK, PMB, RH, WL, TLR, RAP, DAP )
      CALL VVD ( RAP, 5.894305175192448940D0, 1D-10, 'sla_OAP',
     :           'Ar', STATUS )
      CALL VVD ( DAP, 1.406150707974922D0, 1D-10, 'sla_OAP',
     :           'Ad', STATUS )

      CALL sla_OAPQK ( 'R', 2.1D0, -0.345D0, AOPRMS, RAP, DAP )
      CALL VVD ( RAP, 2.10023962776202D0, 1D-10, 'sla_OAPQK',
     :           'Rr', STATUS )
      CALL VVD ( DAP, -0.3452428692888919D0, 1D-10, 'sla_OAPQK',
     :           'Rd', STATUS )
      CALL sla_OAPQK ( 'H', -0.01D0, 1.03D0, AOPRMS, RAP, DAP )
      CALL VVD ( RAP, 1.328731933634564995D0, 1D-10, 'sla_OAPQK',
     :           'Hr', STATUS )
      CALL VVD ( DAP, 1.030091538647746D0, 1D-10, 'sla_OAPQK',
     :           'Hd', STATUS )
      CALL sla_OAPQK ( 'A', 4.321D0, 0.987D0, AOPRMS, RAP, DAP )
      CALL VVD ( RAP, 0.4375507112075065923D0, 1D-10, 'sla_OAPQK',
     :           'Ar', STATUS )
      CALL VVD ( DAP, -0.01520898480744436D0, 1D-10, 'sla_OAPQK',
     :           'Ad', STATUS )

      CALL sla_AOPPAT ( DATE + DS2R, AOPRMS )
      CALL VVD ( AOPRMS(14), 7.602374979243502D0, 1D-8, 'sla_AOPPAT',
     :           ' ', STATUS )

      END

      SUBROUTINE T_BEAR ( STATUS )
*+
*  - - - - - - -
*   T _ B E A R
*  - - - - - - -
*
*  Test sla_BEAR, sla_DBEAR, sla_DPAV, sla_PAV routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  VVD, sla_BEAR, sla_DBEAR,
*           sla_DCS2C, sla_PAV, sla_DPAV.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I
      REAL F1(3), F2(3)
      REAL sla_BEAR, sla_PAV
      DOUBLE PRECISION D1(3), D2(3)
      DOUBLE PRECISION A1, B1, A2, B2
      DOUBLE PRECISION sla_DBEAR, sla_DPAV


      A1 = 1.234D0
      B1 = -0.123D0
      A2 = 2.345D0
      B2 = 0.789D0

      CALL VVD ( DBLE( sla_BEAR ( SNGL( A1 ), SNGL( B1 ), SNGL( A2 ),
     :           SNGL( B2 ) ) ), 0.7045970341781791D0, 1D-6,
     :           'sla_BEAR', ' ', STATUS )
      CALL VVD ( sla_DBEAR ( A1, B1, A2, B2 ), 0.7045970341781791D0,
     :           1D-12, 'sla_DBEAR', ' ', STATUS )
      CALL sla_DCS2C ( A1, B1, D1 )
      CALL sla_DCS2C ( A2, B2, D2 )

      DO I = 1, 3
        F1(I) = SNGL( D1(I) )
        F2(I) = SNGL( D2(I) )
      END DO

      CALL VVD ( DBLE( sla_PAV ( F1, F2 ) ), 0.7045970341781791D0,
     :           1D-6, 'sla_PAV', ' ', STATUS )
      CALL VVD ( sla_DPAV ( D1, D2 ), 0.7045970341781791D0,
     :           1D-12, 'sla_DPAV', ' ', STATUS )

      END

      SUBROUTINE T_CAF2R ( STATUS )
*+
*  - - - - - - - -
*   T _ C A F 2 R
*  - - - - - - - -
*
*  Test sla_CAF2R, sla_DAF2R routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CAF2R, VVD, VIV, sla_DAF2R.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      REAL R
      DOUBLE PRECISION DR

      CALL sla_CAF2R ( 76, 54, 32.1E0, R, J )
      CALL VVD ( DBLE( R ), 1.342313819975276D0, 1D-6, 'sla_CAF2R',
     :           'R', STATUS )
      CALL VIV ( J, 0, 'sla_CAF2R', 'J', STATUS )
      CALL sla_DAF2R ( 76, 54, 32.1D0, DR, J )
      CALL VVD ( DR, 1.342313819975276D0, 1D-12, 'sla_DAF2R',
     :           'R', STATUS )
      CALL VIV ( J, 0, 'sla_CAF2R', 'J', STATUS )

      END

      SUBROUTINE T_CALDJ ( STATUS )
*+
*  - - - - - - - -
*   T _ C A L D J
*  - - - - - - - -
*
*  Test sla_CALDJ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CALDJ, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      DOUBLE PRECISION DJM

      CALL sla_CALDJ ( 1999, 12, 31, DJM, J )
      CALL VVD ( DJM, 51543D0, 0D0, 'sla_CALDJ', ' ', STATUS )

      END

      SUBROUTINE T_CALYD ( STATUS )
*+
*  - - - - - - - -
*   T _ C A L Y D
*  - - - - - - - -
*
*  Test sla_CALYD and sla_CLYD routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CALYD, sla_CLYD, VIV.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER NY, ND, J

      CALL sla_CALYD ( 46, 4, 30, NY, ND, J )
      CALL VIV ( NY, 2046, 'sla_CALYD', 'Y', STATUS )
      CALL VIV ( ND, 120, 'sla_CALYD', 'D', STATUS )
      CALL VIV ( J, 0, 'sla_CALYD', 'J', STATUS )
      CALL sla_CLYD ( -5000, 1, 1, NY, ND, J )
      CALL VIV ( J, 1, 'sla_CLYD', 'illegal year', STATUS )
      CALL sla_CLYD ( 1900, 0, 1, NY, ND, J )
      CALL VIV ( J, 2, 'sla_CLYD', 'illegal month', STATUS )
      CALL sla_CLYD ( 1900, 2, 29, NY, ND, J)
      CALL VIV ( NY, 1900, 'sla_CLYD', 'illegal day (Y)', STATUS )
      CALL VIV ( ND, 61, 'sla_CLYD', 'illegal day (D)', STATUS )
      CALL VIV ( J, 3, 'sla_CLYD', 'illegal day (J)', STATUS )
      CALL sla_CLYD ( 2000, 2, 29, NY, ND, J )
      CALL VIV ( NY, 2000, 'sla_CLYD', 'Y', STATUS )
      CALL VIV ( ND, 60, 'sla_CLYD', 'D', STATUS )
      CALL VIV ( J, 0, 'sla_CLYD', 'J', STATUS )

      END

      SUBROUTINE T_CC2S ( STATUS )
*+
*  - - - - - - -
*   T _ C C 2 S
*  - - - - - - -
*
*  Test sla_CC2S, sla_DCC2S routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CC2S, VVD, sla_DCC2S.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL V(3), A, B
      DOUBLE PRECISION DV(3), DA, DB

      DATA V/100.0, -50.0, 25.0/
      DATA DV/100D0, -50D0, 25D0/

      CALL sla_CC2S ( V, A, B )
      CALL VVD ( DBLE( A), -0.4636476090008061D0, 1D-6, 'sla_CC2S',
     :           'A', STATUS )
      CALL VVD ( DBLE( B ), 0.2199879773954594D0, 1D-6, 'sla_CC2S',
     :           'B', STATUS )

      CALL sla_DCC2S ( DV, DA, DB )
      CALL VVD ( DA, -0.4636476090008061D0, 1D-12, 'sla_DCC2S',
     :           'A', STATUS )
      CALL VVD ( DB, 0.2199879773954594D0, 1D-12, 'sla_DCC2S',
     :           'B', STATUS )

      END

      SUBROUTINE T_CC62S ( STATUS )
*+
*  - - - - - - - -
*   T _ C C 6 2 S
*  - - - - - - - -
*
*  Test sla_CC62S, sla_DC62S routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CC62S, VVD, sla_DC62S.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL V(6), A, B, R, AD, BD, RD
      DOUBLE PRECISION DV(6), DA, DB, DR, DAD, DBD, DRD

      DATA V/100.0, -50.0, 25.0, -0.1, 0.2, 0.7/
      DATA DV/100D0, -50D0, 25D0, -0.1D0, 0.2D0, 0.7D0/

      CALL sla_CC62S ( V, A, B, R, AD, BD, RD )
      CALL VVD ( DBLE( A ), -0.4636476090008061D0, 1D-6, 'sla_CC62S',
     :           'A', STATUS )
      CALL VVD ( DBLE( B ), 0.2199879773954594D0, 1D-6, 'sla_CC62S',
     :           'B', STATUS )
      CALL VVD ( DBLE( R ), 114.564392373896D0, 1D-3, 'sla_CC62S',
     :           'R', STATUS )
      CALL VVD ( DBLE( AD ), 0.001200000000000000D0, 1D-9, 'sla_CC62S',
     :           'AD', STATUS )
      CALL VVD ( DBLE( BD ), 0.006303582107999407D0, 1D-8, 'sla_CC62S',
     :           'BD', STATUS )
      CALL VVD ( DBLE( RD ), -0.02182178902359925D0, 1D-7, 'sla_CC62S',
     :           'RD', STATUS )

      CALL sla_DC62S ( DV, DA, DB, DR, DAD, DBD, DRD )
      CALL VVD ( DA, -0.4636476090008061D0, 1D-6, 'sla_DC62S',
     :           'A', STATUS )
      CALL VVD ( DB, 0.2199879773954594D0, 1D-6, 'sla_DC62S',
     :           'B', STATUS )
      CALL VVD ( DR, 114.564392373896D0, 1D-9, 'sla_DC62S',
     :           'R', STATUS )
      CALL VVD ( DAD, 0.001200000000000000D0, 1D-15, 'sla_DC62S',
     :           'AD', STATUS )
      CALL VVD ( DBD, 0.006303582107999407D0, 1D-14, 'sla_DC62S',
     :           'BD', STATUS )
      CALL VVD ( DRD, -0.02182178902359925D0, 1D-13, 'sla_DC62S',
     :           'RD', STATUS )

      END

      SUBROUTINE T_CD2TF ( STATUS )
*+
*  - - - - - - - -
*   T _ C D 2 T F
*  - - - - - - - -
*
*  Test sla_CD2TF, sla_DD2TF routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CD2TF, VIV, VVD, sla_DD2TF.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER IHMSF(4)
      CHARACTER S

      CALL sla_CD2TF ( 4, -0.987654321E0, S, IHMSF )
      CALL VIV ( ICHAR( S ), ICHAR( '-' ), 'sla_CD2TF', 'S', STATUS )
      CALL VIV ( IHMSF(1), 23, 'sla_CD2TF', '(1)', STATUS )
      CALL VIV ( IHMSF(2), 42, 'sla_CD2TF', '(2)', STATUS )
      CALL VIV ( IHMSF(3), 13, 'sla_CD2TF', '(3)', STATUS )
      CALL VVD ( DFLOAT( IHMSF(4) ), 3333D0, 1000D0, 'sla_CD2TF',
     :           '(4)', STATUS )

      CALL sla_DD2TF ( 4, -0.987654321D0, S, IHMSF )
      CALL VIV ( ICHAR( S ), ICHAR( '-' ), 'sla_DD2TF', 'S', STATUS )
      CALL VIV ( IHMSF(1), 23, 'sla_DD2TF', '(1)', STATUS )
      CALL VIV ( IHMSF(2), 42, 'sla_DD2TF', '(2)', STATUS )
      CALL VIV ( IHMSF(3), 13, 'sla_DD2TF', '(3)', STATUS )
      CALL VIV ( IHMSF(4), 3333, 'sla_DD2TF', '(4)', STATUS )

      END

      SUBROUTINE T_CLDJ ( STATUS )
*+
*  - - - - - - -
*   T _ C L D J
*  - - - - - - -
*
*  Test sla_CLDJ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CLDJ, VVD, VIV.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      DOUBLE PRECISION D

      CALL sla_CLDJ ( 1899, 12, 31, D, J )
      CALL VVD ( D, 15019D0, 0D0, 'sla_CLDJ', 'D', STATUS )
      CALL VIV ( J, 0, 'sla_CLDJ', 'J', STATUS )

      END

      SUBROUTINE T_CR2AF ( STATUS )
*+
*  - - - - - - - -
*   T _ C R 2 A F
*  - - - - - - - -
*
*  Test sla_CR2AF, sla_DR2AF routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CR2AF, VIV, VVD, sla_DR2AF.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER IDMSF(4)
      CHARACTER S

      CALL sla_CR2AF ( 4, 2.345E0, S, IDMSF )
      CALL VIV ( ICHAR( S ), ICHAR( '+' ), 'sla_CR2AF', 'S', STATUS )
      CALL VIV ( IDMSF(1), 134, 'sla_CR2AF', '(1)', STATUS )
      CALL VIV ( IDMSF(2), 21, 'sla_CR2AF', '(2)', STATUS )
      CALL VIV ( IDMSF(3), 30, 'sla_CR2AF', '(3)', STATUS )
      CALL VVD ( DBLE( IDMSF(4) ), 9706D0, 1000D0, 'sla_CR2AF',
     :           '(4)', STATUS )

      CALL sla_DR2AF ( 4, 2.345D0, S, IDMSF )
      CALL VIV ( ICHAR( S ), ICHAR( '+' ), 'sla_DR2AF', 'S', STATUS )
      CALL VIV ( IDMSF(1), 134, 'sla_DR2AF', '(1)', STATUS )
      CALL VIV ( IDMSF(2), 21, 'sla_DR2AF', '(2)', STATUS )
      CALL VIV ( IDMSF(3), 30, 'sla_DR2AF', '(3)', STATUS )
      CALL VIV ( IDMSF(4), 9706, 'sla_DR2AF', '(4)', STATUS )

      END

      SUBROUTINE T_CR2TF ( STATUS )
*+
*  - - - - - - - -
*   T _ C R 2 T F
*  - - - - - - - -
*
*  Test sla_CR2TF, sla_DR2TF routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CR2TF, VIV, VVD, sla_DR2TF.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER IHMSF(4)
      CHARACTER S

      CALL sla_CR2TF ( 4, -3.01234E0, S, IHMSF )
      CALL VIV ( ICHAR( S ), ICHAR( '-' ), 'sla_CR2TF', 'S', STATUS )
      CALL VIV ( IHMSF(1), 11, 'sla_CR2TF', '(1)', STATUS )
      CALL VIV ( IHMSF(2), 30, 'sla_CR2TF', '(2)', STATUS )
      CALL VIV ( IHMSF(3), 22, 'sla_CR2TF', '(3)', STATUS )
      CALL VVD ( DBLE( IHMSF(4) ), 6484D0, 1000D0, 'sla_CR2TF',
     :           '(4)', STATUS )

      CALL sla_DR2TF ( 4, -3.01234D0, S, IHMSF )
      CALL VIV ( ICHAR( S ), ICHAR( '-' ), 'sla_DR2TF', 'S', STATUS )
      CALL VIV ( IHMSF(1), 11, 'sla_DR2TF', '(1)', STATUS )
      CALL VIV ( IHMSF(2), 30, 'sla_DR2TF', '(2)', STATUS )
      CALL VIV ( IHMSF(3), 22, 'sla_DR2TF', '(3)', STATUS )
      CALL VIV ( IHMSF(4), 6484, 'sla_DR2TF', '(4)', STATUS )

      END

      SUBROUTINE T_CS2C6 ( STATUS )
*+
*  - - - - - - - -
*   T _ C S 2 C 6
*  - - - - - - - -
*
*  Test sla_CS2C6, sla_DS2C6 routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CS2C6, VVD, sla_DS2C6.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL V(6)
      DOUBLE PRECISION DV(6)

      CALL sla_CS2C6( -3.21E0, 0.123E0, 0.456E0, -7.8E-6, 9.01E-6,
     :                -1.23E-5, V )
      CALL VVD ( DBLE( V(1) ), -0.4514964673880165D0,
     :           1D-6, 'sla_CS2C6', 'X', STATUS )
      CALL VVD ( DBLE( V(2) ),  0.03093394277342585D0,
     :           1D-6, 'sla_CS2C6', 'Y', STATUS )
      CALL VVD ( DBLE( V(3) ),  0.05594668105108779D0,
     :           1D-6, 'sla_CS2C6', 'Z', STATUS )
      CALL VVD ( DBLE( V(4) ),  1.292270850663260D-5,
     :           1D-6, 'sla_CS2C6', 'XD', STATUS )
      CALL VVD ( DBLE( V(5) ),  2.652814182060692D-6,
     :           1D-6, 'sla_CS2C6', 'YD', STATUS )
      CALL VVD ( DBLE( V(6) ),  2.568431853930293D-6,
     :           1D-6, 'sla_CS2C6', 'ZD', STATUS )

      CALL sla_DS2C6( -3.21D0, 0.123D0, 0.456D0, -7.8D-6, 9.01D-6,
     :                -1.23D-5, DV )
      CALL VVD ( DV(1), -0.4514964673880165D0, 1D-12, 'sla_DS2C6',
     :           'X', STATUS )
      CALL VVD ( DV(2),  0.03093394277342585D0, 1D-12, 'sla_DS2C6',
     :           'Y', STATUS )
      CALL VVD ( DV(3),  0.05594668105108779D0, 1D-12, 'sla_DS2C6',
     :           'Z', STATUS )
      CALL VVD ( DV(4),  1.292270850663260D-5, 1D-12, 'sla_DS2C6',
     :           'XD', STATUS )
      CALL VVD ( DV(5),  2.652814182060692D-6, 1D-12, 'sla_DS2C6',
     :           'YD', STATUS )
      CALL VVD ( DV(6),  2.568431853930293D-6, 1D-12, 'sla_DS2C6',
     :           'ZD', STATUS )

      END

      SUBROUTINE T_CTF2D ( STATUS )
*+
*  - - - - - - - -
*   T _ C T F 2 D
*  - - - - - - - -
*
*  Test sla_CTF2D, sla_DTF2D routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CTF2D, VVD, VIV, sla_DTF2D.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      REAL D
      DOUBLE PRECISION DD

      CALL sla_CTF2D (23, 56, 59.1E0, D, J)
      CALL VVD ( DBLE( D ), 0.99790625D0, 1D-6, 'sla_CTF2D',
     :           'D', STATUS )
      CALL VIV ( J, 0, 'sla_CTF2D', 'J', STATUS )

      CALL sla_DTF2D (23, 56, 59.1D0, DD, J)
      CALL VVD ( DD, 0.99790625D0, 1D-12, 'sla_DTF2D', 'D', STATUS )
      CALL VIV ( J, 0, 'sla_DTF2D', 'J', STATUS )

      END

      SUBROUTINE T_CTF2R ( STATUS )
*+
*  - - - - - - - -
*   T _ C T F 2 R
*  - - - - - - - -
*
*  Test sla_CTF2R, sla_DTF2R routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_CTF2R, VVD, VIV, sla_DTF2R.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      REAL R
      DOUBLE PRECISION DR

      CALL sla_CTF2R (23, 56, 59.1E0, R, J)
      CALL VVD ( DBLE( R ), 6.270029887942679D0, 1D-6, 'sla_CTF2R',
     :           'R', STATUS )
      CALL VIV ( J, 0, 'sla_CTF2R', 'J', STATUS )

      CALL sla_DTF2R (23, 56, 59.1D0, DR, J)
      CALL VVD ( DR, 6.270029887942679D0, 1D-12, 'sla_DTF2R',
     :           'R', STATUS )
      CALL VIV ( J, 0, 'sla_DTF2R', 'J', STATUS )

      END

      SUBROUTINE T_DAT ( STATUS )
*+
*  - - - - - -
*   T _ D A T
*  - - - - - -
*
*  Test sla_DAT, sla_DTT, sla_DT routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_DAT, sla_DTT, sla_DT, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_DAT, sla_DTT, sla_DT


      CALL VVD ( sla_DAT ( 43900D0 ), 18D0, 0D0, 'sla_DAT',
     :           ' ', STATUS )
      CALL VVD ( sla_DTT ( 40404D0 ), 39.709746D0, 1D-12, 'sla_DTT',
     :           ' ', STATUS )
      CALL VVD ( sla_DT ( 500D0 ), 4686.7D0, 1D-10, 'sla_DT',
     :           '500', STATUS )
      CALL VVD ( sla_DT ( 1400D0 ), 408D0, 1D-11, 'sla_DT',
     :           '1400', STATUS )
      CALL VVD ( sla_DT ( 1950D0 ), 27.99145626D0, 1D-12, 'sla_DT',
     :           '1950', STATUS )

      END

      SUBROUTINE T_DBJIN ( STATUS )
*+
*  - - - - - - - -
*   T _ D B J I N
*  - - - - - - - -
*
*  Test sla_DBJIN routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_DBJIN, VVD, VIV.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I, JA, JB
      DOUBLE PRECISION D
      CHARACTER*32 S
      DATA S /'  B1950, , J 2000, B1975 JE     '/

      I = 1
      D = 0D0

      CALL sla_DBJIN ( S, I, D, JA, JB )
      CALL VIV ( I, 9, 'sla_DBJIN', 'I1', STATUS )
      CALL VVD ( D, 1950D0, 0D0, 'sla_DBJIN', 'D1', STATUS )
      CALL VIV ( JA, 0, 'sla_DBJIN', 'JA1', STATUS )
      CALL VIV ( JB, 1, 'sla_DBJIN', 'JB1', STATUS )

      CALL sla_DBJIN ( S, I, D, JA, JB )
      CALL VIV ( I, 11, 'sla_DBJIN', 'I2', STATUS )
      CALL VVD ( D, 1950D0, 0D0, 'sla_DBJIN', 'D2', STATUS )
      CALL VIV ( JA, 1, 'sla_DBJIN', 'JA2', STATUS )
      CALL VIV ( JB, 0, 'sla_DBJIN', 'JB2', STATUS )

      CALL sla_DBJIN ( S, I, D, JA, JB )
      CALL VIV ( I, 19, 'sla_DBJIN', 'I3', STATUS )
      CALL VVD ( D, 2000D0, 0D0, 'sla_DBJIN', 'D3', STATUS )
      CALL VIV ( JA, 0, 'sla_DBJIN', 'JA3', STATUS )
      CALL VIV ( JB, 2, 'sla_DBJIN', 'JB3', STATUS )

      CALL sla_DBJIN ( S, I, D, JA, JB )
      CALL VIV ( I, 26, 'sla_DBJIN', 'I4', STATUS )
      CALL VVD ( D, 1975D0, 0D0, 'sla_DBJIN', 'D4', STATUS )
      CALL VIV ( JA, 0, 'sla_DBJIN', 'JA4', STATUS )
      CALL VIV ( JB, 1, 'sla_DBJIN', 'JB4', STATUS )

      CALL sla_DBJIN ( S, I, D, JA, JB )
      CALL VIV ( I, 26, 'sla_DBJIN', 'I5', STATUS )
      CALL VVD ( D, 1975D0, 0D0, 'sla_DBJIN', 'D5', STATUS )
      CALL VIV ( JA, 1, 'sla_DBJIN', 'JA5', STATUS )
      CALL VIV ( JB, 0, 'sla_DBJIN', 'JB5', STATUS )

      END

      SUBROUTINE T_DJCAL ( STATUS )
*+
*  - - - - - - - -
*   T _ D J C A L
*  - - - - - - - -
*
*  Test sla_DJCAL, sla_DJCL routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_DJCAL, VIV.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER IYDMF(4), J, IY, IM, ID
      DOUBLE PRECISION DJM
      DOUBLE PRECISION F

      DJM = 50123.9999D0

      CALL sla_DJCAL ( 4, DJM, IYDMF, J )
      CALL VIV ( IYDMF(1), 1996, 'sla_DJCAL', 'Y', STATUS )
      CALL VIV ( IYDMF(2), 2, 'sla_DJCAL', 'M', STATUS )
      CALL VIV ( IYDMF(3), 10, 'sla_DJCAL', 'D', STATUS )
      CALL VIV ( IYDMF(4), 9999, 'sla_DJCAL', 'F', STATUS )
      CALL VIV ( J, 0, 'sla_DJCAL', 'J', STATUS )

      CALL sla_DJCL ( DJM, IY, IM, ID, F, J )
      CALL VIV ( IY, 1996, 'sla_DJCL', 'Y', STATUS )
      CALL VIV ( IM, 2, 'sla_DJCL', 'M', STATUS )
      CALL VIV ( ID, 10, 'sla_DJCL', 'D', STATUS )
      CALL VVD ( F, 0.9999D0, 1D-7, 'sla_DJCL', 'F', STATUS )
      CALL VIV ( J, 0, 'sla_DJCL', 'J', STATUS )

      END

      SUBROUTINE T_DMAT ( STATUS )
*+
*  - - - - - - -
*   T _ D M A T
*  - - - - - - -
*
*  Test sla_DMAT routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_DMAT, VVD, VIV.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J, IW(3)
      DOUBLE PRECISION DA(3,3)
      DOUBLE PRECISION DV(3)
      DOUBLE PRECISION DD

      DATA DA/2.22D0,     1.6578D0,     1.380522D0,
     :        1.6578D0,   1.380522D0,   1.22548578D0,
     :        1.380522D0, 1.22548578D0, 1.1356276122D0/
      DATA DV/2.28625D0, 1.7128825D0, 1.429432225D0/

      CALL sla_DMAT ( 3, DA, DV, DD, J, IW )

      CALL VVD ( DA(1,1), 18.02550629769198D0,
     :           1D-10, 'sla_DMAT', 'A(1,1)', STATUS )
      CALL VVD ( DA(1,2), -52.16386644917280607D0,
     :           1D-10, 'sla_DMAT', 'A(1,2)', STATUS )
      CALL VVD ( DA(1,3), 34.37875949717850495D0,
     :           1D-10, 'sla_DMAT', 'A(1,3)', STATUS )
      CALL VVD ( DA(2,1), -52.16386644917280607D0,
     :           1D-10, 'sla_DMAT', 'A(2,1)', STATUS )
      CALL VVD ( DA(2,2), 168.1778099099805627D0,
     :           1D-10, 'sla_DMAT', 'A(2,2)', STATUS )
      CALL VVD ( DA(2,3), -118.0722869694232670D0,
     :           1D-10, 'sla_DMAT', 'A(2,3)', STATUS )
      CALL VVD ( DA(3,1), 34.37875949717850495D0,
     :           1D-10, 'sla_DMAT', 'A(3,1)', STATUS )
      CALL VVD ( DA(3,2), -118.0722869694232670D0,
     :           1D-10, 'sla_DMAT', 'A(3,2)', STATUS )
      CALL VVD ( DA(3,3), 86.50307003740151262D0,
     :           1D-10, 'sla_DMAT', 'A(3,3)', STATUS )
      CALL VVD ( DV(1), 1.002346480763383D0,
     :           1D-12, 'sla_DMAT', 'V(1)', STATUS )
      CALL VVD ( DV(2), 0.03285594016974583489D0,
     :           1D-12, 'sla_DMAT', 'V(2)', STATUS )
      CALL VVD ( DV(3), 0.004760688414885247309D0,
     :           1D-12, 'sla_DMAT', 'V(3)', STATUS )
      CALL VVD ( DD, 0.003658344147359863D0,
     :           1D-12, 'sla_DMAT', 'D', STATUS )
      CALL VIV ( J, 0, 'sla_DMAT', 'J', STATUS )

      END

      SUBROUTINE T_E2H ( STATUS )
*+
*  - - - - - - -
*   T _ E 2 H
*  - - - - - - -
*
*  Test sla_E2H, sla_DE2H, sla_H2E, sla_DH2E routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  All the above plus VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL H, D, P, A, E
      DOUBLE PRECISION DH, DD, DP, DA, DE

      DH = -0.3D0
      DD = -1.1D0
      DP = -0.7D0

      H = SNGL( DH )
      D = SNGL( DD )
      P = SNGL( DP )

      CALL sla_DE2H ( DH, DD, DP, DA, DE )
      CALL VVD ( DA, 2.820087515852369D0, 1D-12, 'sla_DE2H',
     :           'AZ', STATUS )
      CALL VVD ( DE, 1.132711866443304D0, 1D-12, 'sla_DE2H',
     :           'El', STATUS )

      CALL sla_E2H ( H, D, P, A, E )
      CALL VVD ( DBLE( A ), 2.820087515852369D0, 1D-6, 'sla_E2H',
     :           'AZ', STATUS )
      CALL VVD ( DBLE( E ), 1.132711866443304D0, 1D-6, 'sla_E2H',
     :           'El', STATUS )

      CALL sla_DH2E ( DA, DE, DP, DH, DD )
      CALL VVD ( DH, -0.3D0, 1D-12, 'sla_DH2E', 'HA', STATUS )
      CALL VVD ( DD, -1.1D0, 1D-12, 'sla_DH2E', 'DEC', STATUS )

      CALL sla_H2E ( A, E, P, H, D )
      CALL VVD ( DBLE( H ), -0.3D0, 1D-6, 'sla_H2E',
     :           'HA', STATUS )
      CALL VVD ( DBLE( D ), -1.1D0, 1D-6, 'sla_H2E',
     :           'DEC', STATUS )

      END

      SUBROUTINE T_EARTH ( STATUS )
*+
*  - - - - - - - -
*   T _ E A R T H
*  - - - - - - - -
*
*  Test sla_EARTH routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EARTH, VVD.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL PV(6)

      CALL sla_EARTH ( 1978, 174, 0.87E0, PV )

      CALL VVD ( DBLE( PV(1) ), 3.590867086D-2, 1D-6, 'sla_EARTH',
     :           'PV(1)', STATUS )
      CALL VVD ( DBLE( PV(2) ), -9.319285116D-1, 1D-6, 'sla_EARTH',
     :           'PV(2)', STATUS )
      CALL VVD ( DBLE( PV(3) ), -4.041039435D-1, 1D-6, 'sla_EARTH',
     :           'PV(3)', STATUS )
      CALL VVD ( DBLE( PV(4) ), 1.956930055D-7, 1D-13, 'sla_EARTH',
     :           'PV(4)', STATUS )
      CALL VVD ( DBLE( PV(5) ), 5.743797400D-9, 1D-13, 'sla_EARTH',
     :           'PV(5)', STATUS )
      CALL VVD ( DBLE( PV(6) ), 2.512001677D-9, 1D-13, 'sla_EARTH',
     :           'PV(6)', STATUS )

      END

      SUBROUTINE T_ECLEQ ( STATUS )
*+
*  - - - - - - - -
*   T _ E C L E Q
*  - - - - - - - -
*
*  Test sla_ECLEQ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_ECLEQ, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R, D

      CALL sla_ECLEQ ( 1.234D0, -0.123D0, 43210D0, R, D )

      CALL VVD ( R, 1.229910118208851D0, 1D-12, 'sla_ECLEQ',
     :           'RA', STATUS )
      CALL VVD ( D, 0.2638461400411088D0, 1D-12, 'sla_ECLEQ',
     :           'DEC', STATUS )

      END

      SUBROUTINE T_ECMAT ( STATUS )
*+
*  - - - - - - - -
*   T _ E C M A T
*  - - - - - - - -
*
*  Test sla_ECMAT routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_ECMAT, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RM(3,3)

      CALL sla_ECMAT ( 41234D0, RM )

      CALL VVD ( RM(1,1), 1D0, 1D-12, 'sla_ECMAT',
     :           '(1,1)', STATUS )
      CALL VVD ( RM(1,2), 0D0, 1D-12, 'sla_ECMAT',
     :           '(1,2)', STATUS )
      CALL VVD ( RM(1,3), 0D0, 1D-12, 'sla_ECMAT',
     :           '(1,3)', STATUS )
      CALL VVD ( RM(2,1), 0D0, 1D-12, 'sla_ECMAT',
     :           '(2,1)', STATUS )
      CALL VVD ( RM(2,2), 0.917456575085716D0, 1D-12, 'sla_ECMAT',
     :           '(2,2)', STATUS )
      CALL VVD ( RM(2,3), 0.397835937079581D0, 1D-12, 'sla_ECMAT',
     :           '(2,3)', STATUS )
      CALL VVD ( RM(3,1), 0D0, 1D-12, 'sla_ECMAT',
     :           '(3,1)', STATUS )
      CALL VVD ( RM(3,2), -0.397835937079581D0, 1D-12, 'sla_ECMAT',
     :           '(3,2)', STATUS )
      CALL VVD ( RM(3,3), 0.917456575085716D0, 1D-12, 'sla_ECMAT',
     :           '(3,3)', STATUS )

      END

      SUBROUTINE T_ECOR ( STATUS )
*+
*  - - - - - - -
*   T _ E C O R
*  - - - - - - -
*
*  Test sla_ECOR routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_ECOR, VVD.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL RV, Tl

      CALL sla_ECOR ( 2.345E0, -0.567E0, 1995, 306, 0.037E0, RV, Tl )

      CALL VVD ( DBLE( RV ), -19.182460D0, 1D-3, 'sla_ECOR',
     :           'RV', STATUS )
      CALL VVD ( DBLE( Tl ), -120.36632D0, 1D-2, 'sla_ECOR',
     :           'Tl', STATUS )

      END

      SUBROUTINE T_EG50 ( STATUS )
*+
*  - - - - - - -
*   T _ E G 5 0
*  - - - - - - -
*
*  Test sla_EG50 routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EG50, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DL, DB

      CALL sla_EG50 ( 3.012D0, 1.234D0, DL, DB )

      CALL VVD ( DL, 2.305557953813397D0, 1D-12, 'sla_EG50',
     :           'L', STATUS )
      CALL VVD ( DB, 0.7903600886585871D0, 1D-12, 'sla_EG50',
     :           'B', STATUS )

      END

      SUBROUTINE T_EPB ( STATUS )

*+
*  - - - - - -
*   T _ E P B
*  - - - - - -
*
*  Test sla_EPB routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EPB, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_EPB


      CALL VVD ( sla_EPB ( 45123D0 ), 1982.419793168669D0, 1D-8,
     :           'sla_EPB', ' ', STATUS )

      END

      SUBROUTINE T_EPB2D ( STATUS )
*+
*  - - - - - - -
*   T _ E P B 2 D
*  - - - - - - -
*
*  Test sla_EPB2D routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EPB2D, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_EPB2D


      CALL VVD ( sla_EPB2D ( 1975.5D0 ), 42595.5995279655D0, 1D-7,
     :           'sla_EPB2D', ' ', STATUS )

      END

      SUBROUTINE T_EPCO ( STATUS )
*+
*  - - - - - - -
*   T _ E P C O
*  - - - - - - -
*
*  Test sla_EPCO routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EPCO, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_EPCO


      CALL VVD ( sla_EPCO ( 'B', 'J', 2000D0 ), 2000.001277513665D0,
     :           1D-7, 'sla_EPCO', 'BJ', STATUS )
      CALL VVD ( sla_EPCO ( 'J', 'B', 1950D0 ), 1949.999790442300D0,
     :           1D-7, 'sla_EPCO', 'JB', STATUS )
      CALL VVD ( sla_EPCO ( 'J', 'J', 2000D0 ), 2000D0,
     :           1D-7, 'sla_EPCO', 'JJ', STATUS )

      END

      SUBROUTINE T_EPJ ( STATUS )
*+
*  - - - - - -
*   T _ E P J
*  - - - - - -
*
*  Test sla_EPJ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EPJ, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_EPJ


      CALL VVD ( sla_EPJ ( 42999D0 ), 1976.603696098563D0,
     :           1D-7, 'sla_EPJ', ' ', STATUS )

      END

      SUBROUTINE T_EPJ2D ( STATUS )
*+
*  - - - - - - - -
*   T _ E P J 2 D
*  - - - - - - - -
*
*  Test sla_EPJ2D routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EPJ2D, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_EPJ2D


      CALL VVD ( sla_EPJ2D ( 2010.077D0 ), 55225.124250D0,
     :           1D-6, 'sla_EPJ2D', ' ', STATUS )

      END

      SUBROUTINE T_EQECL ( STATUS )
*+
*  - - - - - - - -
*   T _ E Q E C L
*  - - - - - - - -
*
*  Test sla_EQECL routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EQECL, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DL, DB

      CALL sla_EQECL ( 0.789D0, -0.123D0, 46555D0, DL, DB )

      CALL VVD ( DL, 0.7036566430349022D0, 1D-12, 'sla_EQECL',
     :           'L', STATUS )
      CALL VVD ( DB, -0.4036047164116848D0, 1D-12, 'sla_EQECL',
     :           'B', STATUS )

      END

      SUBROUTINE T_EQEQX ( STATUS )
*+
*  - - - - - - - -
*   T _ E Q E Q X
*  - - - - - - - -
*
*  Test sla_EQEQX routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EQEQX, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_EQEQX


      CALL VVD ( sla_EQEQX ( 41234D0 ), 5.376047445838358596D-5,
     :           1D-17, 'sla_EQEQX', ' ', STATUS )

      END

      SUBROUTINE T_EQGAL ( STATUS )
*+
*  - - - - - - - -
*   T _ E Q G A L
*  - - - - - - - -
*
*  Test sla_EQGAL routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EQGAL, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DL, DB

      CALL sla_EQGAL ( 5.67D0, -1.23D0, DL, DB )

      CALL VVD ( DL, 5.612270780904526D0, 1D-12, 'sla_EQGAL',
     :           'DL', STATUS )
      CALL VVD ( DB, -0.6800521449061520D0, 1D-12, 'sla_EQGAL',
     :           'DB', STATUS )

      END

      SUBROUTINE T_ETRMS ( STATUS )
*+
*  - - - - - - - -
*   T _ E T R M S
*  - - - - - - - -
*
*  Test sla_ETRMS routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_ETRMS, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION EV(3)

      CALL sla_ETRMS ( 1976.9D0, EV )

      CALL VVD ( EV(1), -1.621617102537041D-6, 1D-18, 'sla_ETRMS',
     :           'X', STATUS )
      CALL VVD ( EV(2), -3.310070088507914D-7, 1D-18, 'sla_ETRMS',
     :           'Y', STATUS )
      CALL VVD ( EV(3), -1.435296627515719D-7, 1D-18, 'sla_ETRMS',
     :           'Z', STATUS )

      END

      SUBROUTINE T_EVP ( STATUS )
*+
*  - - - - - -
*   T _ E V P
*  - - - - - -
*
*  Test sla_EVP and sla_EPV routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EVP, sla_EPV, VVD.
*
*  Last revision:   21 October 2005
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

      LOGICAL STATUS

      DOUBLE PRECISION DVB(3), DPB(3), DVH(3), DPH(3)

      CALL sla_EVP ( 50100D0, 1990D0, DVB, DPB, DVH, DPH )

      CALL VVD ( DVB(1), -1.807210068604058436D-7, 1D-14, 'sla_EVP',
     :           'DVB(X)', STATUS )
      CALL VVD ( DVB(2), -8.385891022440320D-8, 1D-14, 'sla_EVP',
     :           'DVB(Y)', STATUS )
      CALL VVD ( DVB(3), -3.635846882638055D-8, 1D-14, 'sla_EVP',
     :           'DVB(Z)', STATUS )
      CALL VVD ( DPB(1), -0.4515615297360333D0, 1D-7, 'sla_EVP',
     :           'DPB(X)', STATUS )
      CALL VVD ( DPB(2),  0.8103788166239596D0, 1D-7, 'sla_EVP',
     :           'DPB(Y)', STATUS )
      CALL VVD ( DPB(3),  0.3514505204144827D0, 1D-7, 'sla_EVP',
     :           'DPB(Z)', STATUS )
      CALL VVD ( DVH(1), -1.806354061156890855D-7, 1D-14, 'sla_EVP',
     :           'DVH(X)', STATUS )
      CALL VVD ( DVH(2), -8.383798678086174D-8, 1D-14, 'sla_EVP',
     :           'DVH(Y)', STATUS )
      CALL VVD ( DVH(3), -3.635185843644782D-8, 1D-14, 'sla_EVP',
     :           'DVH(Z)', STATUS )
      CALL VVD ( DPH(1), -0.4478571659918565D0, 1D-7, 'sla_EVP',
     :           'DPH(X)', STATUS )
      CALL VVD ( DPH(2),  0.8036439916076232D0, 1D-7, 'sla_EVP',
     :           'DPH(Y)', STATUS )
      CALL VVD ( DPH(3),  0.3484298459102053D0, 1D-7, 'sla_EVP',
     :           'DPH(Z)', STATUS )

      CALL sla_EPV ( 53411.52501161D0, DPH, DVH, DPB, DVB )

      CALL VVD ( DPH(1), -0.7757238809297653D0, 1D-12, 'sla_EPV',
     :           'DPH(X)', STATUS )
      CALL VVD ( DPH(2), +0.5598052241363390D0, 1D-12, 'sla_EPV',
     :           'DPH(Y)', STATUS )
      CALL VVD ( DPH(3), +0.2426998466481708D0, 1D-12, 'sla_EPV',
     :           'DPH(Z)', STATUS )
      CALL VVD ( DVH(1), -0.0109189182414732D0, 1D-12, 'sla_EPV',
     :           'DVH(X)', STATUS )
      CALL VVD ( DVH(2), -0.0124718726844084D0, 1D-12, 'sla_EPV',
     :           'DVH(Y)', STATUS )
      CALL VVD ( DVH(3), -0.0054075694180650D0, 1D-12, 'sla_EPV',
     :           'DVH(Z)', STATUS )
      CALL VVD ( DPB(1), -0.7714104440491060D0, 1D-12, 'sla_EPV',
     :           'DPB(X)', STATUS )
      CALL VVD ( DPB(2), +0.5598412061824225D0, 1D-12, 'sla_EPV',
     :           'DPB(Y)', STATUS )
      CALL VVD ( DPB(3), +0.2425996277722475D0, 1D-12, 'sla_EPV',
     :           'DPB(Z)', STATUS )
      CALL VVD ( DVB(1), -0.0109187426811683D0, 1D-12, 'sla_EPV',
     :           'DVB(X)', STATUS )
      CALL VVD ( DVB(2), -0.0124652546173285D0, 1D-12, 'sla_EPV',
     :           'DVB(Y)', STATUS )
      CALL VVD ( DVB(3), -0.0054047731809662D0, 1D-12, 'sla_EPV',
     :           'DVB(Z)', STATUS )

      END

      SUBROUTINE T_FITXY ( STATUS )
*+
*  - - - - - - - -
*   T _ F I T X Y
*  - - - - - - - -
*
*  Test sla_FITXY, sla_PXY, sla_INVF, sla_XY2XY, sla_DCMPF routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FITXY, VVD, VIV, sla_PXY, sla_INVF, sla_XY2XY, sla_DCMPF.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J, NPTS

      PARAMETER (NPTS = 8)

      DOUBLE PRECISION XYE(2,NPTS)
      DOUBLE PRECISION XYM(2,NPTS)
      DOUBLE PRECISION COEFFS(6), XYP(2,NPTS), XRMS, YRMS, RRMS,
     :    BKWDS(6), X2, Y2, XZ, YZ, XS, YS, PERP, ORIENT

      DATA XYE/-23.4D0, -12.1D0,   32D0,  -15.3D0,
     :          10.9D0,  23.7D0,   -3D0,   16.1D0,
     :          45D0,  32.5D0,    8.6D0,  -17D0,
     :          15.3D0,  10D0,  121.7D0,   -3.8D0/
      DATA XYM/-23.41D0, 12.12D0,  32.03D0,  15.34D0,
     :          10.93D0,-23.72D0,  -3.01D0, -16.10D0,
     :          44.90D0,-32.46D0,   8.55D0,  17.02D0,
     :          15.31D0,-10.07D0, 120.92D0,   3.81D0/

*  Fit a 4-coeff linear model to relate two sets of (x,y) coordinates.

      CALL sla_FITXY ( 4, NPTS, XYE, XYM, COEFFS, J )
      CALL VVD ( COEFFS(1), -7.938263381515947D-3,
     :           1D-12, 'sla_FITXY', '4/1', STATUS )
      CALL VVD ( COEFFS(2), 1.004640925187200D0,
     :           1D-12, 'sla_FITXY', '4/2', STATUS )
      CALL VVD ( COEFFS(3), 3.976948048238268D-4,
     :           1D-12, 'sla_FITXY', '4/3', STATUS )
      CALL VVD ( COEFFS(4), -2.501031681585021D-2,
     :           1D-12, 'sla_FITXY', '4/4', STATUS )
      CALL VVD ( COEFFS(5), 3.976948048238268D-4,
     :           1D-12, 'sla_FITXY', '4/5', STATUS )
      CALL VVD ( COEFFS(6), -1.004640925187200D0,
     :           1D-12, 'sla_FITXY', '4/6', STATUS )
      CALL VIV ( J, 0, 'sla_FITXY', '4/J', STATUS )

*  Same but 6-coeff.

      CALL sla_FITXY ( 6, NPTS, XYE, XYM, COEFFS, J )
      CALL VVD ( COEFFS(1), -2.617232551841476D-2,
     :           1D-12, 'sla_FITXY', '6/1', STATUS )
      CALL VVD ( COEFFS(2), 1.005634905041421D0,
     :           1D-12, 'sla_FITXY', '6/2', STATUS )
      CALL VVD ( COEFFS(3), 2.133045023329208D-3,
     :           1D-12, 'sla_FITXY', '6/3', STATUS )
      CALL VVD ( COEFFS(4), 3.846993364417779909D-3,
     :           1D-12, 'sla_FITXY', '6/4', STATUS )
      CALL VVD ( COEFFS(5), 1.301671386431460D-4,
     :           1D-12, 'sla_FITXY', '6/5', STATUS )
      CALL VVD ( COEFFS(6), -0.9994827065693964D0,
     :           1D-12, 'sla_FITXY', '6/6', STATUS )
      CALL VIV ( J, 0, 'sla_FITXY', '6/J', STATUS )

*  Compute predicted coordinates and residuals.

      CALL sla_PXY ( NPTS, XYE, XYM, COEFFS, XYP, XRMS, YRMS, RRMS )
      CALL VVD ( XYP(1,1), -23.542232946855340D0,
     :           1D-12, 'sla_PXY', 'X1', STATUS )
      CALL VVD ( XYP(2,1), -12.11293062297230597D0,
     :           1D-12, 'sla_PXY', 'Y1', STATUS )
      CALL VVD ( XYP(1,2), 32.217034593616180D0,
     :           1D-12, 'sla_PXY', 'X2', STATUS )
      CALL VVD ( XYP(2,2), -15.324048471959370D0,
     :           1D-12, 'sla_PXY', 'Y2', STATUS )
      CALL VVD ( XYP(1,3), 10.914821358630950D0,
     :           1D-12, 'sla_PXY', 'X3', STATUS )
      CALL VVD ( XYP(2,3), 23.712999520015880D0,
     :           1D-12, 'sla_PXY', 'Y3', STATUS )
      CALL VVD ( XYP(1,4), -3.087475414568693D0,
     :           1D-12, 'sla_PXY', 'X4', STATUS )
      CALL VVD ( XYP(2,4), 16.09512676604438414D0,
     :           1D-12, 'sla_PXY', 'Y4', STATUS )
      CALL VVD ( XYP(1,5), 45.05759626938414666D0,
     :           1D-12, 'sla_PXY', 'X5', STATUS )
      CALL VVD ( XYP(2,5), 32.45290015313210889D0,
     :           1D-12, 'sla_PXY', 'Y5', STATUS )
      CALL VVD ( XYP(1,6), 8.608310538882801D0,
     :           1D-12, 'sla_PXY', 'X6', STATUS )
      CALL VVD ( XYP(2,6), -17.006235743411300D0,
     :           1D-12, 'sla_PXY', 'Y6', STATUS )
      CALL VVD ( XYP(1,7), 15.348618307280820D0,
     :           1D-12, 'sla_PXY', 'X7', STATUS )
      CALL VVD ( XYP(2,7), 10.07063070741086835D0,
     :           1D-12, 'sla_PXY', 'Y7', STATUS )
      CALL VVD ( XYP(1,8), 121.5833272936291482D0,
     :           1D-12, 'sla_PXY', 'X8', STATUS )
      CALL VVD ( XYP(2,8), -3.788442308260240D0,
     :           1D-12, 'sla_PXY', 'Y8', STATUS )
      CALL VVD ( XRMS ,0.1087247110488075D0,
     :           1D-13, 'sla_PXY', 'XRMS', STATUS )
      CALL VVD ( YRMS, 0.03224481175794666D0,
     :           1D-13, 'sla_PXY', 'YRMS', STATUS )
      CALL VVD ( RRMS, 0.1134054261398109D0,
     :           1D-13, 'sla_PXY', 'RRMS', STATUS )

*  Invert the model.

      CALL sla_INVF ( COEFFS, BKWDS, J )
      CALL VVD ( BKWDS(1), 0.02601750208015891D0,
     :           1D-12, 'sla_INVF', '1', status)
      CALL VVD ( BKWDS(2), 0.9943963945040283D0,
     :           1D-12, 'sla_INVF', '2', status)
      CALL VVD ( BKWDS(3), 0.002122190075497872D0,
     :           1D-12, 'sla_INVF', '3', status)
      CALL VVD ( BKWDS(4), 0.003852372795357474353D0,
     :           1D-12, 'sla_INVF', '4', status)
      CALL VVD ( BKWDS(5), 0.0001295047252932767D0,
     :           1D-12, 'sla_INVF', '5', status)
      CALL VVD ( BKWDS(6), -1.000517284779212D0,
     :           1D-12, 'sla_INVF', '6', status)
      CALL VIV ( J, 0, 'sla_INVF', 'J', STATUS )

*  Transform one x,y.

      CALL sla_XY2XY ( 44.5D0, 32.5D0, COEFFS, X2, Y2 )
      CALL VVD ( X2, 44.793904912083030D0,
     :           1D-11, 'sla_XY2XY', 'X', status)
      CALL VVD ( Y2, -32.473548532471330D0,
     :           1D-11, 'sla_XY2XY', 'Y', status)

*  Decompose the fit into scales etc.

      CALL sla_DCMPF ( COEFFS, XZ, YZ, XS, YS, PERP, ORIENT )
      CALL VVD ( XZ, -0.0260175020801628646D0,
     :           1D-12, 'sla_DCMPF', 'XZ', status)
      CALL VVD ( YZ, -0.003852372795357474353D0,
     :           1D-12, 'sla_DCMPF', 'YZ', status)
      CALL VVD ( XS, -1.00563491346569D0,
     :           1D-12, 'sla_DCMPF', 'XS', status)
      CALL VVD ( YS, 0.999484982684761D0,
     :           1D-12, 'sla_DCMPF', 'YS', status)
      CALL VVD ( PERP,-0.002004707996156263D0,
     :           1D-12, 'sla_DCMPF', 'P', status)
      CALL VVD ( ORIENT, 3.14046086182333D0,
     :           1D-12, 'sla_DCMPF', 'O', status)

      END

      SUBROUTINE T_FK425 ( STATUS )
*+
*  - - - - - - - -
*   T _ F K 4 2 5
*  - - - - - - - -
*
*  Test sla_FK425 routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FK425.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R2000, D2000, DR2000, DD2000, P2000, V2000

      CALL sla_FK425 ( 1.234D0, -0.123D0, -1D-5, 2D-6, 0.5D0,
     :                 20D0, R2000, D2000, DR2000, DD2000, P2000,
     :                 V2000 )

      CALL VVD ( R2000, 1.244117554618727D0, 1D-12, 'sla_FK425',
     :           'R', STATUS )
      CALL VVD ( D2000, -0.1213164254458709D0, 1D-12, 'sla_FK425',
     :           'D', STATUS )
      CALL VVD ( DR2000, -9.964265838268711D-6, 1D-17, 'sla_FK425',
     :           'DR', STATUS )
      CALL VVD ( DD2000, 2.038065265773541D-6, 1D-17, 'sla_FK425',
     :           'DD', STATUS )
      CALL VVD ( P2000, 0.4997443812415410D0, 1D-12, 'sla_FK425',
     :           'P', STATUS )
      CALL VVD ( V2000, 20.010460915421010D0, 1D-11, 'sla_FK425',
     :           'V', STATUS )

      END

      SUBROUTINE T_FK45Z ( STATUS )
*+
*  - - - - - - - -
*   T _ F K 4 5 Z
*  - - - - - - - -
*
*  Test sla_FK45Z routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FK45Z.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R2000, D2000

      CALL sla_FK45Z ( 1.234D0, -0.123D0, 1984D0, R2000, D2000 )

      CALL VVD ( R2000, 1.244616510731691D0, 1D-12, 'sla_FK45Z',
     :           'R', STATUS )
      CALL VVD ( D2000, -0.1214185839586555D0, 1D-12, 'sla_FK45Z',
     :           'D', STATUS )

      END

      SUBROUTINE T_FK524 ( STATUS )
*+
*  - - - - - - - -
*   T _ F K 5 2 4
*  - - - - - - - -
*
*  Test sla_FK524 routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FK524.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R1950, D1950, DR1950, DD1950, P1950, V1950

      CALL sla_FK524 ( 4.567D0, -1.23D0, -3D-5, 8D-6, 0.29D0,
     :                 -35D0, R1950, D1950, DR1950, DD1950, P1950,
     :                 V1950 )

      CALL VVD ( R1950, 4.543778603272084D0, 1D-12, 'sla_FK524',
     :           'R', STATUS )
      CALL VVD ( D1950, -1.229642790187574D0, 1D-12, 'sla_FK524',
     :           'D', STATUS )
      CALL VVD ( DR1950, -2.957873121769244D-5, 1D-17, 'sla_FK524',
     :           'DR', STATUS )
      CALL VVD ( DD1950, 8.117725309659079D-6, 1D-17, 'sla_FK524',
     :           'DD', STATUS )
      CALL VVD ( P1950, 0.2898494999992917D0, 1D-12, 'sla_FK524',
     :           'P', STATUS )
      CALL VVD ( V1950, -35.026862824252680D0, 1D-11, 'sla_FK524',
     :           'V', STATUS )

      END

      SUBROUTINE T_FK52H ( STATUS )
*+
*  - - - - - - - -
*   T _ F K 5 2 H
*  - - - - - - - -
*
*  Test sla_FK52H, sla_H2FK5, sla_FK5HZ, sla_HFK5Z routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FK524, sla_H2FK5.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R5, D5, DR5, DD5, RH, DH, DRH, DDH

      CALL sla_FK52H ( 1.234D0, -0.987D0, 1D-6, -2D-6, RH, DH, DRH,
     :                 DDH )
      CALL VVD ( RH, 1.234000000272122558D0, 1D-13, 'sla_FK52H',
     :           'R', STATUS )
      CALL VVD ( DH, -0.9869999235218543959D0, 1D-13, 'sla_FK52H',
     :           'D', STATUS )
      CALL VVD ( DRH, 0.000000993178295D0, 1D-13, 'sla_FK52H',
     :           'DR', STATUS )
      CALL VVD ( DDH, -0.000001997665915D0, 1D-13, 'sla_FK52H',
     :           'DD', STATUS )
      CALL sla_H2FK5 ( RH, DH, DRH, DDH, r5, D5, DR5, DD5 )
      CALL VVD ( R5, 1.234D0, 1D-13, 'sla_H2FK5', 'R', STATUS )
      CALL VVD ( D5, -0.987D0, 1D-13, 'sla_H2FK5', 'D', STATUS )
      CALL VVD ( DR5, 1D-6, 1D-13, 'sla_H2FK5', 'DR', STATUS )
      CALL VVD ( DD5, -2D-6, 1D-13, 'sla_H2FK5', 'DD', STATUS )
      CALL sla_FK5HZ ( 1.234D0, -0.987D0, 1980D0, RH, DH )
      CALL VVD ( RH, 1.234000136713611301D0, 1D-13, 'sla_FK5HZ',
     :           'R', STATUS )
      CALL VVD ( DH, -0.9869999702020807601D0, 1D-13, 'sla_FK5HZ',
     :           'D', STATUS )
      CALL sla_HFK5Z ( RH, DH, 1980D0, R5, D5, DR5, DD5 )
      CALL VVD ( R5, 1.234D0, 1D-13, 'sla_HFK5Z', 'R', STATUS )
      CALL VVD ( D5, -0.987D0, 1D-13, 'sla_HFK5Z', 'D', STATUS )
      CALL VVD ( DR5, 0.000000006822074D0, 1D-13, 'sla_HFK5Z',
     :           'DR', STATUS )
      CALL VVD ( DD5, -0.000000002334012D0, 1D-13, 'sla_HFK5Z',
     :           'DD', STATUS )

      END

      SUBROUTINE T_FK54Z ( STATUS )
*+
*  - - - - - - - -
*   T _ F K 5 4 Z
*  - - - - - - - -
*
*  Test sla_FK54Z routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FK54Z.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R1950, D1950, DR1950, DD1950

      CALL sla_FK54Z ( 0.001D0, -1.55D0, 1900D0, R1950, D1950,
     :                 DR1950, DD1950 )

      CALL VVD ( R1950, 6.271585543439484D0, 1D-12, 'sla_FK54Z',
     :           'R', STATUS )
      CALL VVD ( D1950, -1.554861715330319D0, 1D-12, 'sla_FK54Z',
     :           'D', STATUS )
      CALL VVD ( DR1950, -4.175410876044916011D-8, 1D-20, 'sla_FK54Z',
     :           'DR', STATUS )
      CALL VVD ( DD1950, 2.118595098308522D-8, 1D-20, 'sla_FK54Z',
     :           'DD', STATUS )

      END

      SUBROUTINE T_FLOTIN ( STATUS )
*+
*  - - - - - - - - -
*   T _ F L O T I N
*  - - - - - - - - -
*
*  Test sla_FLOTIN, sla_DFLTIN routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_FLOTIN, VVD, VIV, sla_DFLTIN.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I, J
      REAL FV
      DOUBLE PRECISION DV
      CHARACTER*33 S
      DATA S /'  12.345, , -0 1E3-4 2000  E     '/

      I = 1
      FV = 0.0

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 10, 'sla_FLOTIN', 'V5', STATUS )
      CALL VVD ( DBLE( FV ), 12.345D0, 1D-5, 'sla_FLOTIN',
     :           'V1', STATUS )
      CALL VIV ( J, 0, 'sla_FLOTIN', 'J1', STATUS )

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 12, 'sla_FLOTIN', 'I2', STATUS )
      CALL VVD ( DBLE( FV ), 12.345D0, 1D-5, 'sla_FLOTIN',
     :           'V2', STATUS )
      CALL VIV ( J, 1, 'sla_FLOTIN', 'J2', STATUS )

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 16, 'sla_FLOTIN', 'I3', STATUS )
      CALL VVD ( DBLE( FV ), 0D0, 0D0, 'sla_FLOTIN', 'V3', STATUS )
      CALL VIV ( J, -1, 'sla_FLOTIN', 'J3', STATUS )

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 19, 'sla_FLOTIN', 'I4', STATUS )
      CALL VVD ( DBLE( FV), 1000D0, 0D0, 'sla_FLOTIN', 'V4', STATUS )
      CALL VIV ( J, 0, 'sla_FLOTIN', 'J4', STATUS )

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 22, 'sla_FLOTIN', 'I5', STATUS )
      CALL VVD ( DBLE( FV ), -4D0, 0D0, 'sla_FLOTIN', 'V5', STATUS )
      CALL VIV ( J, -1, 'sla_FLOTIN', 'J5', STATUS )

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 28, 'sla_FLOTIN', 'I6', STATUS )
      CALL VVD ( DBLE( FV ), 2000D0, 0D0, 'sla_FLOTIN',
     :           'V6', STATUS )
      CALL VIV ( J, 0, 'sla_FLOTIN', 'J6', STATUS )

      CALL sla_FLOTIN ( S, I, FV, J )
      CALL VIV ( I, 34, 'sla_FLOTIN', 'I7', STATUS )
      CALL VVD ( DBLE( FV ), 2000D0, 0D0, 'sla_FLOTIN',
     :           'V7', STATUS )
      CALL VIV ( J, 2, 'sla_FLOTIN', 'J7', STATUS )

      I = 1
      DV = 0D0

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 10, 'sla_DFLTIN', 'I1', STATUS )
      CALL VVD ( DV, 12.345D0, 1D-12, 'sla_DFLTIN', 'V1', STATUS )
      CALL VIV ( J, 0, 'sla_DFLTIN', 'J1', STATUS )

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 12, 'sla_DFLTIN', 'I2', STATUS )
      CALL VVD ( DV, 12.345D0, 1D-12, 'sla_DFLTIN', 'V2', STATUS )
      CALL VIV ( J, 1, 'sla_DFLTIN', 'J2', STATUS )

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 16, 'sla_DFLTIN', 'I3', STATUS )
      CALL VVD ( DV, 0D0, 0D0, 'sla_DFLTIN', 'V3', STATUS )
      CALL VIV ( J, -1, 'sla_DFLTIN', 'J3', STATUS )

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 19, 'sla_DFLTIN', 'I4', STATUS )
      CALL VVD ( DV, 1000D0, 0D0, 'sla_DFLTIN', 'V4', STATUS )
      CALL VIV ( J, 0, 'sla_DFLTIN', 'J4', STATUS )

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 22, 'sla_DFLTIN', 'I5', STATUS )
      CALL VVD ( DV, -4D0, 0D0, 'sla_DFLTIN', 'V5', STATUS )
      CALL VIV ( J, -1, 'sla_DFLTIN', 'J5', STATUS )

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 28, 'sla_DFLTIN', 'I6', STATUS )
      CALL VVD ( DV, 2000D0, 0D0, 'sla_DFLTIN', 'V6', STATUS )
      CALL VIV ( J, 0, 'sla_DFLTIN', 'J6', STATUS )

      CALL sla_DFLTIN ( S, I, DV, J )
      CALL VIV ( I, 34, 'sla_DFLTIN', 'I7', STATUS )
      CALL VVD ( DV, 2000D0, 0D0, 'sla_DFLTIN', 'V7', STATUS )
      CALL VIV ( J, 2, 'sla_DFLTIN', 'J7', STATUS )

      END

      SUBROUTINE T_GALEQ ( STATUS )
*+
*  - - - - - - - -
*   T _ G A L E Q
*  - - - - - - - -
*
*  Test sla_GALEQ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_GALEQ, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DR, DD

      CALL sla_GALEQ ( 5.67D0, -1.23D0, DR, DD )

      CALL VVD ( DR, 0.04729270418071426D0, 1D-12, 'sla_GALEQ',
     :           'DR', STATUS )
      CALL VVD ( DD, -0.7834003666745548D0, 1D-12, 'sla_GALEQ',
     :           'DD', STATUS )

      END

      SUBROUTINE T_GALSUP ( STATUS )
*+
*  - - - - - - - - -
*   T _ G A L S U P
*  - - - - - - - - -
*
*  Test sla_GALSUP routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_GALSUP, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DSL, DSB

      CALL sla_GALSUP ( 6.1D0, -1.4D0, DSL, DSB )

      CALL VVD ( DSL, 4.567933268859171D0, 1D-12, 'sla_GALSUP',
     :           'DSL', STATUS )
      CALL VVD ( DSB, -0.01862369899731829D0, 1D-12, 'sla_GALSUP',
     :           'DSB', STATUS )

      END

      SUBROUTINE T_GE50 ( STATUS )
*+
*  - - - - - - -
*   T _ G E 5 0
*  - - - - - - -
*
*  Test sla_GE50 routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_GE50, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DR, DD

      CALL sla_GE50 ( 6.1D0, -1.55D0, DR, DD )

      CALL VVD ( DR, 0.1966825219934508D0, 1D-12, 'sla_GE50',
     :           'DR', STATUS )
      CALL VVD ( DD, -0.4924752701678960D0, 1D-12, 'sla_GE50',
     :           'DD', STATUS )

      END

      SUBROUTINE T_GMST ( STATUS )
*+
*  - - - - - - -
*   T _ G M S T
*  - - - - - - -
*
*  Test sla_GMST and sla_GMSTA routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_GMST, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_GMST, sla_GMSTA


      CALL VVD ( sla_GMST ( 43999.999D0 ), 3.9074971356487318D0,
     :           1D-9, 'sla_GMST', ' ', STATUS )
      CALL VVD ( sla_GMSTA ( 43999D0, 0.999D0 ),
     :           3.9074971356487318D0, 1D-12, 'sla_GMSTA', ' ', STATUS )

      END

      SUBROUTINE T_INTIN ( STATUS )
*+
*  - - - - - - - -
*   T _ I N T I N
*  - - - - - - - -
*
*  Test sla_INTIN routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_INTIN, VIV.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER*4 N
      INTEGER I, J
      CHARACTER*28 S
      DATA S /'  -12345, , -0  2000  +     '/

      I = 1
      N = 0

      CALL sla_INTIN ( S, I, N, J )
      CALL VIV ( I, 10, 'sla_INTIN', 'I1', STATUS )
      CALL VLV ( N, -12345, 'sla_INTIN', 'V1', STATUS )
      CALL VIV ( J, -1, 'sla_INTIN', 'J1', STATUS )

      CALL sla_INTIN ( S, I, N, J )
      CALL VIV ( I, 12, 'sla_INTIN', 'I2', STATUS )
      CALL VLV ( N, -12345, 'sla_INTIN', 'V2', STATUS )
      CALL VIV ( J, 1, 'sla_INTIN', 'J2', STATUS )

      CALL sla_INTIN ( S, I, N, J )
      CALL VIV ( I, 17, 'sla_INTIN', 'I3', STATUS )
      CALL VLV ( N, 0, 'sla_INTIN', 'V3', STATUS )
      CALL VIV ( J, -1, 'sla_INTIN', 'J3', STATUS )

      CALL sla_INTIN ( S, I, N, J )
      CALL VIV ( I, 23, 'sla_INTIN', 'I4', STATUS )
      CALL VLV ( N, 2000, 'sla_INTIN', 'V4', STATUS )
      CALL VIV ( J, 0, 'sla_INTIN', 'J4', STATUS )

      CALL sla_INTIN ( S, I, N, J )
      CALL VIV ( I, 29, 'sla_INTIN', 'I5', STATUS )
      CALL VLV ( N, 2000, 'sla_INTIN', 'V5', STATUS )
      CALL VIV ( J, 2, 'sla_INTIN', 'J5', STATUS )

      END

      SUBROUTINE T_KBJ ( STATUS )
*+
*  - - - - - -
*   T _ K B J
*  - - - - - -
*
*  Test sla_KBJ routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_KBJ, VCS, VIV.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      DOUBLE PRECISION E
      CHARACTER K
      DATA K /'?'/

      E = 1950D0
      CALL sla_KBJ ( -1, E, K, J )
      CALL VCS ( K, ' ', 'sla_KBJ', 'JB1', STATUS )
      CALL VIV ( J, 1, 'sla_KBJ', 'J1', STATUS )
      CALL sla_KBJ ( 0, E, K, J )
      CALL VCS ( K, 'B', 'sla_KBJ', 'JB2', STATUS )
      CALL VIV ( J, 0, 'sla_KBJ', 'J2', STATUS )
      CALL sla_KBJ ( 1, E, K, J )
      CALL VCS ( K, 'B', 'sla_KBJ', 'JB3', STATUS )
      CALL VIV ( J, 0, 'sla_KBJ', 'J3', STATUS )
      CALL sla_KBJ ( 2, E, K, J )
      CALL VCS ( K, 'J', 'sla_KBJ', 'JB4', STATUS )
      CALL VIV ( J, 0, 'sla_KBJ', 'J4', STATUS )
      CALL sla_KBJ ( 3, E, K, J )
      CALL VCS ( K, ' ', 'sla_KBJ', 'JB5', STATUS )
      CALL VIV ( J, 1, 'sla_KBJ', 'J5', STATUS )

      E = 2000D0
      CALL sla_KBJ ( 0, E, K, J )
      CALL VCS ( K, 'J', 'sla_KBJ', 'JB6', STATUS )
      CALL VIV ( J, 0, 'sla_KBJ', 'J6', STATUS )
      CALL sla_KBJ ( 1, E, K, J )
      CALL VCS ( K, 'B', 'sla_KBJ', 'jB7', STATUS )
      CALL VIV ( J, 0, 'sla_KBJ', 'J7', STATUS )
      CALL sla_KBJ ( 2, E, K, J )
      CALL VCS ( K, 'J', 'sla_KBJ', 'JB8', STATUS )
      CALL VIV ( J, 0, 'sla_KBJ', 'J8', STATUS )

      END

      SUBROUTINE T_MAP ( STATUS )
*+
*  - - - - - -
*   T _ M A P
*  - - - - - -
*
*  Test sla_MAP, sla_MAPPA, sla_MAPQK, sla_MAPQKZ routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_MAP, sla_MAPPA, sla_MAPQK, sla_MAPQKZ, VVD.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RA, DA, AMPRMS(21)

      CALL sla_MAP ( 6.123D0, -0.999D0, 1.23D-5, -0.987D-5,
     :               0.123D0, 32.1D0, 1999D0, 43210.9D0, RA, DA )

      CALL VVD ( RA, 6.117130429775647D0, 1D-12, 'sla_MAP',
     :           'RA', STATUS )
      CALL VVD ( DA, -1.000880769038632D0, 1D-12, 'sla_MAP',
     :           'DA', STATUS )

      CALL sla_MAPPA ( 2020D0, 45012.3D0, AMPRMS )

      CALL VVD ( AMPRMS(1), -37.884188911704310D0,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(1)', STATUS )
      CALL VVD ( AMPRMS(2),  -0.7888341859486424D0,
     :           1D-7, 'sla_MAPPA', 'AMPRMS(2)', STATUS )
      CALL VVD ( AMPRMS(3),   0.5405321789059870D0,
     :           1D-7, 'sla_MAPPA', 'AMPRMS(3)', STATUS )
      CALL VVD ( AMPRMS(4),   0.2340784267119091D0,
     :           1D-7, 'sla_MAPPA', 'AMPRMS(4)', STATUS )
      CALL VVD ( AMPRMS(5),  -0.8067807553217332071D0,
     :           1D-7, 'sla_MAPPA', 'AMPRMS(5)', STATUS )
      CALL VVD ( AMPRMS(6),   0.5420884771236513880D0,
     :           1D-7, 'sla_MAPPA', 'AMPRMS(6)', STATUS )
      CALL VVD ( AMPRMS(7),   0.2350423277034460899D0,
     :           1D-7, 'sla_MAPPA', 'AMPRMS(7)', STATUS )
      CALL VVD ( AMPRMS(8),   1.999729469227807D-8,
     :           1D-12, 'sla_MAPPA', 'AMPRMS(8)', STATUS )
      CALL VVD ( AMPRMS(9),  -6.035531043691568494D-5,
     :           1D-12, 'sla_MAPPA', 'AMPRMS(9)', STATUS )
      CALL VVD ( AMPRMS(10), -7.381891582591552377D-5,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(10)', STATUS )
      CALL VVD ( AMPRMS(11), -3.200897749853207412D-5,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(11)', STATUS )
      CALL VVD ( AMPRMS(12),  0.9999999949417148D0,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(12)', STATUS )
      CALL VVD ( AMPRMS(13),  0.9999566751478850D0,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(13)', STATUS )
      CALL VVD ( AMPRMS(14), -8.537361890149777D-3,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(14)', STATUS )
      CALL VVD ( AMPRMS(15), -3.709619811228171D-3,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(15)', STATUS )
      CALL VVD ( AMPRMS(16), 8.537308717676752D-3,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(16)', STATUS )
      CALL VVD ( AMPRMS(17),  0.9999635560607690D0,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(17)', STATUS )
      CALL VVD ( AMPRMS(18), -3.016886324169151D-5,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(18)', STATUS )
      CALL VVD ( AMPRMS(19),  3.709742180572510D-3,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(19)', STATUS )
      CALL VVD ( AMPRMS(20), -1.502613373498668D-6,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(20)', STATUS )
      CALL VVD ( AMPRMS(21),  0.9999931188816729D0,
     :           1D-11, 'sla_MAPPA', 'AMPRMS(21)', STATUS )

      CALL sla_MAPQK ( 1.234D0, -0.987D0, -1.2D-5, -0.99D0,
     :                 0.75D0, -23.4D0, AMPRMS, RA, DA )

      CALL VVD ( RA, 1.223337584930993D0, 1D-11, 'sla_MAPQK',
     :           'RA', STATUS )
      CALL VVD ( DA, 0.5558838650379129D0, 1D-11, 'sla_MAPQK',
     :           'DA', STATUS )

      CALL sla_MAPQKZ ( 6.012D0, 1.234D0, AMPRMS, RA, DA )

      CALL VVD ( RA, 6.006091119756597D0, 1D-11, 'sla_MAPQKZ',
     :           'RA', STATUS )
      CALL VVD ( DA, 1.23045846622498D0, 1D-11, 'sla_MAPQKZ',
     :           'DA', STATUS )

      END

      SUBROUTINE T_MOON ( STATUS )
*+
*  - - - - - - -
*   T _ M O O N
*  - - - - - - -
*
*  Test sla_MOON and sla_DMOON routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_MOON, sla_DMOON, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL PV(6)

      CALL sla_MOON ( 1999, 365, 0.9E0, PV )

      CALL VVD ( DBLE( PV(1) ), -2.155729505970773D-3, 1D-6,
     :           'sla_MOON', '(1)', STATUS )
      CALL VVD ( DBLE( PV(2) ), -1.538107758633427D-3, 1D-6,
     :           'sla_MOON', '(2)', STATUS )
      CALL VVD ( DBLE( PV(3) ), -4.003940552689305D-4, 1D-6 ,
     :           'sla_MOON', '(3)', STATUS )
      CALL VVD ( DBLE( PV(4) ),  3.629209419071314D-9, 1D-12,
     :           'sla_MOON', '(4)', STATUS )
      CALL VVD ( DBLE( PV(5) ), -4.989667166259157D-9, 1D-12,
     :           'sla_MOON', '(5)', STATUS )
      CALL VVD ( DBLE( PV(6) ), -2.160752457288307D-9, 1D-12,
     :           'sla_MOON', '(6)', STATUS )

      END

      SUBROUTINE T_NUT ( STATUS )
*+
*  - - - - - -
*   T _ N U T
*  - - - - - -
*
*  Test sla_NUT, sla_NUTC, sla_NUTC80 routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_NUT, sla_NUTC, VVD.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RMATN(3,3), DPSI, DEPS, EPS0

      CALL sla_NUT ( 46012.34D0, RMATN )

      CALL VVD ( RMATN(1,1),  9.999999969492166D-1, 1D-12,
     :           'sla_NUT', '(1,1)', STATUS )
      CALL VVD ( RMATN(1,2),  7.166577986249302D-5, 1D-12,
     :           'sla_NUT', '(1,2)', STATUS )
      CALL VVD ( RMATN(1,3),  3.107382973077677D-5, 1D-12,
     :           'sla_NUT', '(1,3)', STATUS )
      CALL VVD ( RMATN(2,1), -7.166503970900504D-5, 1D-12,
     :           'sla_NUT', '(2,1)', STATUS )
      CALL VVD ( RMATN(2,2),  9.999999971483732D-1, 1D-12,
     :           'sla_NUT', '(2,2)', STATUS )
      CALL VVD ( RMATN(2,3), -2.381965032461830D-5, 1D-12,
     :           'sla_NUT', '(2,3)', STATUS )
      CALL VVD ( RMATN(3,1), -3.107553669598237D-5, 1D-12,
     :           'sla_NUT', '(3,1)', STATUS )
      CALL VVD ( RMATN(3,2),  2.381742334472628D-5, 1D-12,
     :           'sla_NUT', '(3,2)', STATUS )
      CALL VVD ( RMATN(3,3),  9.999999992335206818D-1, 1D-12,
     :           'sla_NUT', '(3,3)', STATUS )

      CALL sla_NUTC ( 50123.4D0, DPSI, DEPS, EPS0 )

      CALL VVD ( DPSI, 3.523550954747999709D-5, 1D-17, 'sla_NUTC',
     :           'DPSI', STATUS )
      CALL VVD ( DEPS, -4.143371566683342D-5, 1D-17, 'sla_NUTC',
     :           'DEPS', STATUS )
      CALL VVD ( EPS0, 0.4091014592901651D0, 1D-12, 'sla_NUTC',
     :           'EPS0', STATUS )

      CALL sla_NUTC80 ( 50123.4D0, DPSI, DEPS, EPS0 )

      CALL VVD ( DPSI, 3.537714281665945321D-5, 1D-17, 'sla_NUTC80',
     :           'DPSI', STATUS )
      CALL VVD ( DEPS, -4.140590085987148317D-5, 1D-17, 'sla_NUTC80',
     :           'DEPS', STATUS )
      CALL VVD ( EPS0, 0.4091016349007751D0, 1D-12, 'sla_NUTC80',
     :           'EPS0', STATUS )

      END

      SUBROUTINE T_OBS ( STATUS )
*+
*  - - - - - -
*   T _ O B S
*  - - - - - -
*
*  Test sla_OBS routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_OBS, err, VVD.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER N
      DOUBLE PRECISION W, P, H
      CHARACTER*10 C
      CHARACTER*40 NAME

      N = 0
      C = 'MMT'
      CALL sla_OBS ( N, C, NAME, W, P, H )
      CALL VCS ( C, 'MMT', 'sla_OBS', '1/C', STATUS )
      CALL VCS ( NAME, 'MMT 6.5m, Mt Hopkins', 'sla_OBS', '1/NAME',
     :           STATUS )
      CALL VVD ( W, 1.935300584055477D0, 1D-8, 'sla_OBS',
     :           '1/W', STATUS )
      CALL VVD ( P, 0.5530735081550342238D0, 1D-10, 'sla_OBS',
     :           '1/P', STATUS )
      CALL VVD ( H, 2608D0, 1D-10, 'sla_OBS',
     :           '1/H', STATUS )

      N = 61
      CALL sla_OBS ( N, C, NAME, W, P, H )
      CALL VCS ( C, 'KECK1', 'sla_OBS', '2/C', STATUS )
      CALL VCS ( NAME, 'Keck 10m Telescope #1', 'sla_OBS',
     :           '2/NAME', STATUS )
      CALL VVD ( W, 2.713545757918895D0, 1D-8, 'sla_OBS',
     :           '2/W', STATUS )
      CALL VVD ( P, 0.3460280563536619D0, 1D-8, 'sla_OBS',
     :           '2/P', STATUS )
      CALL VVD ( H, 4160D0, 1D-10, 'sla_OBS',
     :           '2/H', STATUS )

      N = 83
      CALL sla_OBS ( N, C, NAME, W, P, H )
      CALL VCS ( C, 'MAGELLAN2', 'sla_OBS', '3/C', STATUS )
      CALL VCS ( NAME, 'Magellan 2, 6.5m, Las Campanas',
     :           'sla_OBS', '3/NAME', STATUS )
      CALL VVD ( W, 1.233819305534497D0, 1D-8, 'sla_OBS',
     :           '3/W', STATUS )
      CALL VVD ( P, -0.506389344359954D0, 1D-8, 'sla_OBS',
     :           '3/P', STATUS )
      CALL VVD ( H, 2408D0, 1D-10, 'sla_OBS',
     :           '3/H', STATUS )

      N = 85
      CALL sla_OBS ( N, C, NAME, W, P, H )
      CALL VCS ( NAME, '?', 'sla_OBS', '4/NAME', STATUS )

      END

      SUBROUTINE T_PA ( STATUS )
*+
*  - - - - -
*   T _ P A
*  - - - - -
*
*  Test sla_PA routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PA, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_PA


      CALL VVD ( sla_PA ( -1.567D0, 1.5123D0, 0.987D0 ),
     :           -1.486288540423851D0, 1D-12, 'sla_PA', ' ', STATUS )
      CALL VVD ( sla_PA ( 0D0, 0.789D0, 0.789D0 ),
     :           0D0, 0D0, 'sla_PA', 'zenith', STATUS )

      END

      SUBROUTINE T_PCD ( STATUS )
*+
*  - - - - - -
*   T _ P C D
*  - - - - - -
*
*  Test sla_PCD, sla_UNPCD routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PCD, VVD, sla_UNPCD.
*
*  Last revision:   4 September 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DISCO, X, Y

      DISCO = 178.585D0
      X = 0.0123D0
      Y = -0.00987D0

      CALL sla_PCD ( DISCO, X, Y )
      CALL VVD ( X, 0.01284630845735895D0, 1D-14, 'sla_PCD',
     :           'X', STATUS )
      CALL VVD ( Y, -0.01030837922553926D0, 1D-14, 'sla_PCD',
     :           'Y', STATUS )

      CALL sla_UNPCD ( DISCO, X, Y )
      CALL VVD ( X, 0.0123D0, 1D-14, 'sla_UNPCD',
     :           'X', STATUS )
      CALL VVD ( Y, -0.00987D0, 1D-14, 'sla_UNPCD',
     :           'Y', STATUS )

      END

      SUBROUTINE T_PDA2H ( STATUS )
*+
*  - - - - - - - -
*   T _ P D A 2 H
*  - - - - - - - -
*
*  Test sla_PDA2H routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PDA2H, VVD.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J1, J2
      DOUBLE PRECISION H1, H2

      CALL sla_PDA2H ( -0.51D0, -1.31D0, 3.1D0, H1, J1, H2, J2 )
      CALL VVD ( H1, -0.1161784556585304927D0, 1D-14, 'sla_PDA2H',
     :           'H1', STATUS )
      CALL VIV ( J1, 0, 'sla_PDA2H', 'J1', STATUS )
      CALL VVD ( H2, -2.984787179226459D0, 1D-13, 'sla_PDA2H',
     :           'H2', STATUS )
      CALL VIV ( J2, 0, 'sla_PDA2H', 'J2', STATUS )

      END

      SUBROUTINE T_PDQ2H ( STATUS )
*+
*  - - - - - - - -
*   T _ P D Q 2 H
*  - - - - - - - -
*
*  Test sla_PDQ2H routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PDQ2H, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J1, J2
      DOUBLE PRECISION H1, H2

      CALL sla_PDQ2H ( 0.9D0, 0.2D0, 0.1D0, H1, J1, H2, J2 )
      CALL VVD ( H1, 0.1042809894435257D0, 1D-14, 'sla_PDQ2H',
     :           'H1', STATUS )
      CALL VIV ( J1, 0, 'sla_PDQ2H', 'J1', STATUS )
      CALL VVD ( H2, 2.997450098818439D0, 1D-13, 'sla_PDQ2H',
     :           'H2', STATUS )
      CALL VIV ( J2, 0, 'sla_PDQ2H', 'J2', STATUS )

      END

      SUBROUTINE T_PERCOM ( STATUS )
*+
*  - - - - - - - - -
*   T _ P E R C O M
*  - - - - - - - - -
*
*  Test sla_COMBN, sla_PERMUT routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_COMBN, VIV, sla_PERMUT.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER LIST(3), I, J, ISTATE(4), IORDER(4)

      LIST(1) = 0

      DO I = 1, 11
         CALL sla_COMBN ( 3, 5, LIST, J )
      END DO

      CALL VIV ( J, 1, 'sla_COMBN', 'J', STATUS )
      CALL VIV ( LIST(1), 1, 'sla_COMBN', 'LIST(1)', STATUS )
      CALL VIV ( LIST(2), 2, 'sla_COMBN', 'LIST(2)', STATUS )
      CALL VIV ( LIST(3), 3, 'sla_COMBN', 'LIST(3)', STATUS )

      ISTATE(1) = -1

      DO I = 1, 25
         CALL sla_PERMUT ( 4, ISTATE, IORDER, J )
      END DO

      CALL VIV ( J, 1, 'sla_PERMUT', 'J', STATUS )
      CALL VIV ( IORDER(1), 4, 'sla_PERMUT', 'IORDER(1)', STATUS )
      CALL VIV ( IORDER(2), 3, 'sla_PERMUT', 'IORDER(2)', STATUS )
      CALL VIV ( IORDER(3), 2, 'sla_PERMUT', 'IORDER(3)', STATUS )
      CALL VIV ( IORDER(4), 1, 'sla_PERMUT', 'IORDER(4)', STATUS )

      END

      SUBROUTINE T_PLANET ( STATUS )
*+
*  - - - - - - - - -
*   T _ P L A N E T
*  - - - - - - - - -
*
*  Test sla_EL2UE, sla_PERTEL, sla_PERTUE, sla_PLANEL, sla_PLANET,
*  sla_PLANTE, sla_PLANTU, sla_PV2EL, sla_PV2UE, sla_RDPLAN, sla_UE2EL
*  and sla_UE2PV routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_EL2UE, sla_PERTEL, sla_PERTUE, sla_PLANEL, sla_PLANET,
*           sla_PLANTE, sla_PLANTU, sla_PV2EL, sla_PV2UE, sla_RDPLAN,
*           sla_UE2EL, sla_UE2PV, VIV, VVD.
*
*  Last revision:   22 October 2005
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

      LOGICAL STATUS

      INTEGER J, JFORM
      DOUBLE PRECISION U(13), PV(6), RA, DEC, R, DIAM, EPOCH, ORBINC,
     :    ANODE, PERIH, AORQ, E, AORL, DM


      CALL sla_EL2UE ( 50000D0, 1, 49000D0, 0.1D0, 2D0, 0.2D0,
     :                 3D0, 0.05D0, 3D0, 0.003312D0, U, J )
      CALL VVD ( U(1), 1.000878908362435284D0, 1D-12, 'sla_EL2UE',
     :           'U(1)', STATUS )
      CALL VVD ( U(2), -0.3336263027874777288D0, 1D-12, 'sla_EL2UE',
     :           'U(2)', STATUS )
      CALL VVD ( U(3), 50000D0, 1D-12, 'sla_EL2UE',
     :           'U(3)', STATUS )
      CALL VVD ( U(4), 2.840425801310305210D0, 1D-12, 'sla_EL2UE',
     :           'U(4)', STATUS )
      CALL VVD ( U(5), 0.1264380368035014224D0, 1D-12, 'sla_EL2UE',
     :           'U(5)', STATUS )
      CALL VVD ( U(6), -0.2287711835229143197D0, 1D-12, 'sla_EL2UE',
     :           'U(6)', STATUS )
      CALL VVD ( U(7), -0.01301062595106185195D0, 1D-12, 'sla_EL2UE',
     :           'U(7)', STATUS )
      CALL VVD ( U(8), 0.5657102158104651697D0, 1D-12, 'sla_EL2UE',
     :           'U(8)', STATUS )
      CALL VVD ( U(9), 0.2189745287281794885D0, 1D-12, 'sla_EL2UE',
     :           'U(9)', STATUS )
      CALL VVD ( U(10), 2.852427310959998500D0, 1D-12, 'sla_EL2UE',
     :           'U(10)', STATUS )
      CALL VVD ( U(11), -0.01552349065435120900D0, 1D-12, 'sla_EL2UE',
     :           'U(11)', STATUS )
      CALL VVD ( U(12), 50000D0, 1D-12, 'sla_EL2UE',
     :           'U(12)', STATUS )
      CALL VVD ( U(13), 0D0, 1D-12, 'sla_EL2UE',
     :           'U(13)', STATUS )
      CALL VIV ( J, 0, 'sla_EL2UE', 'J', STATUS )

      CALL sla_PERTEL ( 2, 43000D0, 43200D0, 43000D0,
     :                  0.2D0, 3D0, 4D0, 5D0, 0.02D0, 6D0,
     :                  EPOCH, ORBINC, ANODE, PERIH, AORQ, E, AORL, J )
      CALL VVD ( EPOCH, 43200D0, 1D-10, 'sla_PERTEL',
     :           'EPOCH', STATUS )
      CALL VVD ( ORBINC, 0.1995661466545422381D0, 1D-7, 'sla_PERTEL',
     :           'ORBINC', STATUS )
      CALL VVD ( ANODE, 2.998052737821591215D0, 1D-7, 'sla_PERTEL',
     :           'ANODE', STATUS )
      CALL VVD ( PERIH, 4.009516448441143636D0, 1D-6, 'sla_PERTEL',
     :           'PERIH', STATUS )
      CALL VVD ( AORQ, 5.014216294790922323D0, 1D-7, 'sla_PERTEL',
     :           'AORQ', STATUS )
      CALL VVD ( E, 0.02281386258309823607D0, 1D-7, 'sla_PERTEL',
     :           'E', STATUS )
      CALL VVD ( AORL, 0.01735248648779583748D0, 1D-6, 'sla_PERTEL',
     :           'AORL', STATUS )
      CALL VIV ( J, 0, 'sla_PERTEL', 'J', STATUS )

      CALL sla_PERTUE ( 50100D0, U, J )
      CALL VVD ( U(1), 1.000000000000000D0, 1D-12, 'sla_PERTUE',
     :           'U(1)', STATUS )
      CALL VVD ( U(2), -0.3329769417028020949D0, 1D-11, 'sla_PERTUE',
     :           'U(2)', STATUS )
      CALL VVD ( U(3), 50100D0, 1D-12, 'sla_PERTUE',
     :           'U(3)', STATUS )
      CALL VVD ( U(4), 2.638884303608524597D0, 1D-11, 'sla_PERTUE',
     :           'U(4)', STATUS )
      CALL VVD ( U(5), 1.070994304747824305D0, 1D-11, 'sla_PERTUE',
     :           'U(5)', STATUS )
      CALL VVD ( U(6), 0.1544112080167568589D0, 1D-11, 'sla_PERTUE',
     :           'U(6)', STATUS )
      CALL VVD ( U(7), -0.2188240619161439344D0, 1D-11, 'sla_PERTUE',
     :           'U(7)', STATUS )
      CALL VVD ( U(8), 0.5207557453451906385D0, 1D-11, 'sla_PERTUE',
     :           'U(8)', STATUS )
      CALL VVD ( U(9), 0.2217782439275216936D0, 1D-11, 'sla_PERTUE',
     :           'U(9)', STATUS )
      CALL VVD ( U(10), 2.852118859689216658D0, 1D-11, 'sla_PERTUE',
     :           'U(10)', STATUS )
      CALL VVD ( U(11), 0.01452010174371893229D0, 1D-11, 'sla_PERTUE',
     :           'U(11)', STATUS )
      CALL VVD ( U(12), 50100D0, 1D-12, 'sla_PERTUE',
     :           'U(12)', STATUS )
      CALL VVD ( U(13), 0D0, 1D-12, 'sla_PERTUE',
     :           'U(13)', STATUS )
      CALL VIV ( J, 0, 'sla_PERTUE', 'J', STATUS )

      CALL sla_PLANEL ( 50600D0, 2, 50500D0, 0.1D0, 3D0, 5D0,
     :                  2D0, 0.3D0, 4D0, 0D0, PV, J )
      CALL VVD ( PV(1), 1.947628959288897677D0, 1D-12, 'sla_PLANEL',
     :           'PV(1)', STATUS )
      CALL VVD ( PV(2), -1.013736058752235271D0, 1D-12, 'sla_PLANEL',
     :           'PV(2)', STATUS )
      CALL VVD ( PV(3), -0.3536409947732733647D0, 1D-12, 'sla_PLANEL',
     :           'PV(3)', STATUS )
      CALL VVD ( PV(4), 2.742247411571786194D-8, 1D-19, 'sla_PLANEL',
     :           'PV(4)', STATUS )
      CALL VVD ( PV(5), 1.170467244079075911D-7, 1D-19, 'sla_PLANEL',
     :           'PV(5)', STATUS )
      CALL VVD ( PV(6), 3.709878268217564005D-8, 1D-19, 'sla_PLANEL',
     :           'PV(6)', STATUS )
      CALL VIV ( J, 0, 'sla_PLANEL', 'J', STATUS )

      CALL sla_PLANET ( 1D6, 0, PV, J )
      CALL VVD ( PV(1), 0D0, 0D0, 'sla_PLANET',
     :           'PV(1) 1', STATUS )
      CALL VVD ( PV(2), 0D0, 0D0, 'sla_PLANET',
     :           'PV(2) 1', STATUS )
      CALL VVD ( PV(3), 0D0, 0D0, 'sla_PLANET',
     :           'PV(3) 1', STATUS )
      CALL VVD ( PV(4), 0D0, 0D0, 'sla_PLANET',
     :           'PV(4) 1', STATUS )
      CALL VVD ( PV(5), 0D0, 0D0, 'sla_PLANET',
     :           'PV(5) 1', STATUS )
      CALL VVD ( PV(6), 0D0, 0D0, 'sla_PLANET',
     :           'PV(6) 1', STATUS )
      CALL VIV ( J, -1, 'sla_PLANET', 'J 1', STATUS )

      CALL sla_PLANET ( 1D6, 10, PV, J )
      CALL VIV ( J, -1, 'sla_PLANET', 'J 2', STATUS )

      CALL sla_PLANET ( -320000D0, 3, PV, J )
      CALL VVD ( PV(1), 0.9308038666827242603D0, 1D-11, 'sla_PLANET',
     :           'PV(1) 3', STATUS )
      CALL VVD ( PV(2), 0.3258319040252137618D0, 1D-11, 'sla_PLANET',
     :           'PV(2) 3', STATUS )
      CALL VVD ( PV(3), 0.1422794544477122021D0, 1D-11, 'sla_PLANET',
     :           'PV(3) 3', STATUS )
      CALL VVD ( PV(4), -7.441503423889371696D-8, 1D-17, 'sla_PLANET',
     :           'PV(4) 3', STATUS )
      CALL VVD ( PV(5), 1.699734557528650689D-7, 1D-17, 'sla_PLANET',
     :           'PV(5) 3', STATUS )
      CALL VVD ( PV(6), 7.415505123001430864D-8, 1D-17, 'sla_PLANET',
     :           'PV(6) 3', STATUS )
      CALL VIV ( J, 1, 'sla_PLANET', 'J 3', STATUS )

      CALL sla_PLANET ( 43999.9D0, 1, PV, J )
      CALL VVD ( PV(1), 0.2945293959257422246D0, 1D-11, 'sla_PLANET',
     :           'PV(1) 4', STATUS )
      CALL VVD ( PV(2), -0.2452204176601052181D0, 1D-11, 'sla_PLANET',
     :           'PV(2) 4', STATUS )
      CALL VVD ( PV(3), -0.1615427700571978643D0, 1D-11, 'sla_PLANET',
     :           'PV(3) 4', STATUS )
      CALL VVD ( PV(4), 1.636421147459047057D-7, 1D-18, 'sla_PLANET',
     :           'PV(4) 4', STATUS )
      CALL VVD ( PV(5), 2.252949422574889753D-7, 1D-18, 'sla_PLANET',
     :           'PV(5) 4', STATUS )
      CALL VVD ( PV(6), 1.033542799062371839D-7, 1D-18, 'sla_PLANET',
     :           'PV(6) 4', STATUS )
      CALL VIV ( J, 0, 'sla_PLANET', 'J 4', STATUS )

      CALL sla_PLANTE ( 50600D0, -1.23D0, 0.456D0, 2, 50500D0,
     :                  0.1D0, 3D0, 5D0, 2D0, 0.3D0, 4D0,
     :                  0D0, RA, DEC, R, J )
      CALL VVD ( RA, 6.222958101333794007D0, 1D-10, 'sla_PLANTE',
     :           'RA', STATUS )
      CALL VVD ( DEC, 0.01142220305739771601D0, 1D-10, 'sla_PLANTE',
     :           'DEC', STATUS )
      CALL VVD ( R, 2.288902494080167624D0, 1D-8, 'sla_PLANTE',
     :           'R', STATUS )
      CALL VIV ( J, 0, 'sla_PLANTE', 'J', STATUS )

      U(1) = 1.0005D0
      U(2) = -0.3D0
      U(3) = 55000D0
      U(4) = 2.8D0
      U(5) = 0.1D0
      U(6) = -0.2D0
      U(7) = -0.01D0
      U(8) = 0.5D0
      U(9) = 0.22D0
      U(10) = 2.8D0
      U(11) = -0.015D0
      U(12) = 55001D0
      U(13) = 0D0

      CALL sla_PLANTU ( 55001D0, -1.23D0, 0.456D0, U, RA, DEC, R, J )
      CALL VVD ( RA, 0.3531814831241686647D0, 1D-9, 'sla_PLANTU',
     :           'RA', STATUS )
      CALL VVD ( DEC, 0.06940344580567131328D0, 1D-9, 'sla_PLANTU',
     :           'DEC', STATUS )
      CALL VVD ( R, 3.031687170873274464D0, 1D-8, 'sla_PLANTU',
     :           'R', STATUS )
      CALL VIV ( J, 0, 'sla_PLANTU', 'J', STATUS )

      PV(1) = 0.3D0
      PV(2) = -0.2D0
      PV(3) = 0.1D0
      PV(4) = -0.9D-7
      PV(5) = 0.8D-7
      PV(6) = -0.7D-7
      CALL sla_PV2EL ( PV, 50000D0, 0.00006D0, 1,
     :                 JFORM, EPOCH, ORBINC, ANODE, PERIH,
     :                 AORQ, E, AORL, DM, J )
      CALL VIV ( JFORM, 1, 'sla_PV2EL', 'JFORM', STATUS )
      CALL VVD ( EPOCH, 50000D0, 1D-10, 'sla_PV2EL',
     :           'EPOCH', STATUS )
      CALL VVD ( ORBINC, 1.52099895268912D0, 1D-12, 'sla_PV2EL',
     :           'ORBINC', STATUS )
      CALL VVD ( ANODE, 2.720503180538650D0, 1D-12, 'sla_PV2EL',
     :           'ANODE', STATUS )
      CALL VVD ( PERIH, 2.194081512031836D0, 1D-12, 'sla_PV2EL',
     :           'PERIH', STATUS )
      CALL VVD ( AORQ, 0.2059371035373771D0, 1D-12, 'sla_PV2EL',
     :           'AORQ', STATUS )
      CALL VVD ( E, 0.9866822985810528D0, 1D-12, 'sla_PV2EL',
     :           'E', STATUS )
      CALL VVD ( AORL, 0.2012758344836794D0, 1D-12, 'sla_PV2EL',
     :           'AORL', STATUS )
      CALL VVD ( DM, 0.1840740507951820D0, 1D-12, 'sla_PV2EL',
     :           'DM', STATUS )
      CALL VIV ( J, 0, 'sla_PV2EL', 'J', STATUS )

      CALL sla_PV2UE ( PV, 50000D0, 0.00006D0, U, J )
      CALL VVD ( U(1), 1.00006D0, 1D-12, 'sla_PV2UE',
     :           'U(1)', STATUS )
      CALL VVD ( U(2), -4.856142884511782D0, 1D-12, 'sla_PV2UE',
     :           'U(2)', STATUS )
      CALL VVD ( U(3), 50000D0, 1D-12, 'sla_PV2UE',
     :           'U(3)', STATUS )
      CALL VVD ( U(4), 0.3D0, 1D-12, 'sla_PV2UE',
     :           'U(4)', STATUS )
      CALL VVD ( U(5), -0.2D0, 1D-12, 'sla_PV2UE',
     :           'U(5)', STATUS )
      CALL VVD ( U(6), 0.1D0, 1D-12, 'sla_PV2UE',
     :           'U(6)', STATUS )
      CALL VVD ( U(7), -0.4520378601821727D0, 1D-12, 'sla_PV2UE',
     :           'U(7)', STATUS )
      CALL VVD ( U(8), 0.4018114312730424D0, 1D-12, 'sla_PV2UE',
     :           'U(8)', STATUS )
      CALL VVD ( U(9), -.3515850023639121D0, 1D-12, 'sla_PV2UE',
     :           'U(9)', STATUS )
      CALL VVD ( U(10), 0.3741657386773941D0, 1D-12, 'sla_PV2UE',
     :           'U(10)', STATUS )
      CALL VVD ( U(11), -0.2511321445456515D0, 1D-12, 'sla_PV2UE',
     :           'U(11)', STATUS )
      CALL VVD ( U(12), 50000D0, 1D-12, 'sla_PV2UE',
     :           'U(12)', STATUS )
      CALL VVD ( U(13), 0D0, 1D-12, 'sla_PV2UE',
     :           'U(13)', STATUS )
      CALL VIV ( J, 0, 'sla_PV2UE', 'J', STATUS )

      CALL sla_RDPLAN ( 40999.9D0, 0, 0.1D0, -0.9D0, RA, DEC, DIAM )
      CALL VVD ( RA, 5.772270359389275837D0, 1D-7, 'sla_RDPLAN',
     :           'RA 0', STATUS )
      CALL VVD ( DEC, -0.2089207338795416192D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 0', STATUS )
      CALL VVD ( DIAM, 9.415338935229717875D-3, 1D-14, 'sla_RDPLAN',
     :           'DIAM 0', STATUS )
      CALL sla_RDPLAN ( 41999.9D0, 1, 1.1D0, -0.9D0, RA, DEC, DIAM )
      CALL VVD ( RA, 3.866363420052936653D0, 1D-7, 'sla_RDPLAN',
     :           'RA 1', STATUS )
      CALL VVD ( DEC, -0.2594430577550113130D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 1', STATUS )
      CALL VVD ( DIAM, 4.638468996795023071D-5, 1D-14, 'sla_RDPLAN',
     :           'DIAM 1', STATUS )
      CALL sla_RDPLAN ( 42999.9D0, 2, 2.1D0, 0.9D0, RA, DEC, DIAM )
      CALL VVD ( RA, 2.695383203184077378D0, 1D-7, 'sla_RDPLAN',
     :           'RA 2', STATUS )
      CALL VVD ( DEC, 0.2124044506294805126D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 2', STATUS )
      CALL VVD ( DIAM, 4.892222838681000389D-5, 1D-14, 'sla_RDPLAN',
     :           'DIAM 2', STATUS )
      CALL sla_RDPLAN ( 43999.9D0, 3, 3.1D0, 0.9D0, RA, DEC, DIAM )
      CALL VVD ( RA, 2.908326678461540165D0, 1D-7, 'sla_RDPLAN',
     :           'RA 3', STATUS )
      CALL VVD ( DEC, 0.08729783126905579385D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 3', STATUS )
      CALL VVD ( DIAM, 8.581305866034962476D-3, 1D-14, 'sla_RDPLAN',
     :           'DIAM 3', STATUS )
      CALL sla_RDPLAN ( 44999.9D0, 4, -0.1D0, 1.1D0, RA, DEC, DIAM )
      CALL VVD ( RA, 3.429840787472851721D0, 1D-7, 'sla_RDPLAN',
     :           'RA 4', STATUS )
      CALL VVD ( DEC, -0.06979851055261161013D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 4', STATUS )
      CALL VVD ( DIAM, 4.540536678439300199D-5, 1D-14, 'sla_RDPLAN',
     :           'DIAM 4', STATUS )
      CALL sla_RDPLAN ( 45999.9D0, 5, -1.1D0, 0.1D0, RA, DEC, DIAM )
      CALL VVD ( RA, 4.864669466449422548D0, 1D-7, 'sla_RDPLAN',
     :           'RA 5', STATUS )
      CALL VVD ( DEC, -0.4077714497908953354D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 5', STATUS )
      CALL VVD ( DIAM, 1.727945579027815576D-4, 1D-14, 'sla_RDPLAN',
     :           'DIAM 5', STATUS )
      CALL sla_RDPLAN ( 46999.9D0, 6, -2.1D0, -0.1D0, RA, DEC, DIAM )
      CALL VVD ( RA, 4.432929829176388766D0, 1D-7, 'sla_RDPLAN',
     :           'RA 6', STATUS )
      CALL VVD ( DEC, -0.3682820877854730530D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 6', STATUS )
      CALL VVD ( DIAM, 8.670829016099083311D-5, 1D-14, 'sla_RDPLAN',
     :           'DIAM 6', STATUS )
      CALL sla_RDPLAN ( 47999.9D0, 7, -3.1D0, -1.1D0, RA, DEC, DIAM )
      CALL VVD ( RA, 4.894972492286818487D0, 1D-7, 'sla_RDPLAN',
     :           'RA 7', STATUS )
      CALL VVD ( DEC, -0.4084068901053653125D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 7', STATUS )
      CALL VVD ( DIAM, 1.793916783975974163D-5, 1D-14, 'sla_RDPLAN',
     :           'DIAM 7', STATUS )
      CALL sla_RDPLAN ( 48999.9D0, 8, 0D0, 0D0, RA, DEC, DIAM )
      CALL VVD ( RA, 5.066050284760144000D0, 1D-7, 'sla_RDPLAN',
     :           'RA 8', STATUS )
      CALL VVD ( DEC, -0.3744690779683850609D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 8', STATUS )
      CALL VVD ( DIAM, 1.062210086082700563D-5, 1D-14, 'sla_RDPLAN',
     :           'DIAM 8', STATUS )
      CALL sla_RDPLAN ( 49999.9D0, 9, 0D0, 0D0, RA, DEC, DIAM )
      CALL VVD ( RA, 4.179543143097200945D0, 1D-7, 'sla_RDPLAN',
     :           'RA 9', STATUS )
      CALL VVD ( DEC, -0.1258021632894033300D0, 1D-7, 'sla_RDPLAN',
     :           'DEC 9', STATUS )
      CALL VVD ( DIAM, 5.034057475664904352D-7, 1D-14, 'sla_RDPLAN',
     :           'DIAM 9', STATUS )

      CALL sla_UE2EL ( U, 1, JFORM, EPOCH, ORBINC, ANODE, PERIH,
     :                 AORQ, E, AORL, DM, J )
      CALL VIV ( JFORM, 1, 'sla_UE2EL', 'JFORM', STATUS )
      CALL VVD ( EPOCH, 50000.00000000000D0, 1D-10, 'sla_PV2EL',
     :           'EPOCH', STATUS )
      CALL VVD ( ORBINC, 1.520998952689120D0, 1D-12, 'sla_UE2EL',
     :           'ORBINC', STATUS )
      CALL VVD ( ANODE, 2.720503180538650D0, 1D-12, 'sla_UE2EL',
     :           'ANODE', STATUS )
      CALL VVD ( PERIH, 2.194081512031836D0, 1D-12, 'sla_UE2EL',
     :           'PERIH', STATUS )
      CALL VVD ( AORQ, 0.2059371035373771D0, 1D-12, 'sla_UE2EL',
     :           'AORQ', STATUS )
      CALL VVD ( E, 0.9866822985810528D0, 1D-12, 'sla_UE2EL',
     :           'E', STATUS )
      CALL VVD ( AORL, 0.2012758344836794D0, 1D-12, 'sla_UE2EL',
     :           'AORL', STATUS )
      CALL VIV ( J, 0, 'sla_UE2EL', 'J', STATUS )

      CALL sla_UE2PV ( 50010D0, U, PV, J )
      CALL VVD ( U(1), 1.00006D0, 1D-12, 'sla_UE2PV',
     :           'U(1)', STATUS )
      CALL VVD ( U(2), -4.856142884511782111D0, 1D-12, 'sla_UE2PV',
     :           'U(2)', STATUS )
      CALL VVD ( U(3), 50000D0, 1D-12, 'sla_UE2PV',
     :           'U(3)', STATUS )
      CALL VVD ( U(4), 0.3D0, 1D-12, 'sla_UE2PV',
     :           'U(4)', STATUS )
      CALL VVD ( U(5), -0.2D0, 1D-12, 'sla_UE2PV',
     :           'U(5)', STATUS )
      CALL VVD ( U(6), 0.1D0, 1D-12, 'sla_UE2PV',
     :           'U(6)', STATUS )
      CALL VVD ( U(7), -0.4520378601821727110D0, 1D-12, 'sla_UE2PV',
     :           'U(7)', STATUS )
      CALL VVD ( U(8), 0.4018114312730424097D0, 1D-12, 'sla_UE2PV',
     :           'U(8)', STATUS )
      CALL VVD ( U(9), -0.3515850023639121085D0, 1D-12, 'sla_UE2PV',
     :           'U(9)', STATUS )
      CALL VVD ( U(10), 0.3741657386773941386D0, 1D-12, 'sla_UE2PV',
     :           'U(10)', STATUS )
      CALL VVD ( U(11), -0.2511321445456515061D0, 1D-12, 'sla_UE2PV',
     :           'U(11)', STATUS )
      CALL VVD ( U(12), 50010.00000000000D0, 1D-12, 'sla_UE2PV',
     :           'U(12)', STATUS )
      CALL VVD ( U(13), 0.7194308220038886856D0, 1D-12, 'sla_UE2PV',
     :           'U(13)', STATUS )
      CALL VVD ( PV(1), 0.07944764084631667011D0, 1D-12, 'sla_UE2PV',
     :           'PV(1)', STATUS )
      CALL VVD ( PV(2), -0.04118141077419014775D0, 1D-12, 'sla_UE2PV',
     :           'PV(2)', STATUS )
      CALL VVD ( PV(3), 0.002915180702063625400D0, 1D-12, 'sla_UE2PV',
     :           'PV(3)', STATUS )
      CALL VVD ( PV(4), -0.6890132370721108608D-6, 1D-18,'sla_UE2PV',
     :           'PV(4)', STATUS )
      CALL VVD ( PV(5), 0.4326690733487621457D-6, 1D-18, 'sla_UE2PV',
     :           'PV(5)', STATUS )
      CALL VVD ( PV(6), -0.1763249096254134306D-6, 1D-18, 'sla_UE2PV',
     :           'PV(6)', STATUS )
      CALL VIV ( J, 0, 'sla_UE2PV', 'J', STATUS )

      END

      SUBROUTINE T_PM ( STATUS )
*+
*  - - - - -
*   T _ P M
*  - - - - -
*
*  Test sla_PM routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PM, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION R1, D1

      CALL sla_PM ( 5.43D0, -0.87D0, -0.33D-5, 0.77D-5, 0.7D0,
     :              50.3D0*365.2422D0/365.25D0, 1899D0, 1943D0,
     :              R1, D1 )
      CALL VVD ( R1, 5.429855087793875D0, 1D-12, 'sla_PM',
     :           'R', STATUS )
      CALL VVD ( D1, -0.8696617307805072D0, 1D-12, 'sla_PM',
     :           'D', STATUS )

      END

      SUBROUTINE T_POLMO ( STATUS )
*+
*  - - - - - - - -
*   T _ P O L M O
*  - - - - - - - -
*
*  Test sla_POLMO routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_POLMO, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION ELONG, PHI, DAZ

      CALL sla_POLMO ( 0.7D0, -0.5D0, 1D-6, -2D-6, ELONG, PHI, DAZ )

      CALL VVD ( ELONG,  0.7000004837322044D0, 1D-12, 'sla_POLMO',
     :           'ELONG', STATUS )
      CALL VVD ( PHI, -0.4999979467222241D0, 1D-12, 'sla_POLMO',
     :           'PHI', STATUS )
      CALL VVD ( DAZ,  1.008982781275728D-6, 1D-12, 'sla_POLMO',
     :           'DAZ', STATUS )

      END

      SUBROUTINE T_PREBN ( STATUS )
*+
*  - - - - - - - -
*   T _ P R E B N
*  - - - - - - - -
*
*  Test sla_PREBN routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PREBN, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RMATP(3,3)

      CALL sla_PREBN ( 1925D0, 1975D0, RMATP )

      CALL VVD ( RMATP(1,1),  9.999257613786738D-1, 1D-12,
     :           'sla_PREBN', '(1,1)', STATUS )
      CALL VVD ( RMATP(1,2), -1.117444640880939D-2, 1D-12,
     :           'sla_PREBN', '(1,2)', STATUS )
      CALL VVD ( RMATP(1,3), -4.858341150654265D-3, 1D-12,
     :           'sla_PREBN', '(1,3)', STATUS )
      CALL VVD ( RMATP(2,1),  1.117444639746558D-2, 1D-12,
     :           'sla_PREBN', '(2,1)', STATUS )
      CALL VVD ( RMATP(2,2),  9.999375635561940D-1, 1D-12,
     :           'sla_PREBN', '(2,2)', STATUS )
      CALL VVD ( RMATP(2,3), -2.714797892626396D-5, 1D-12,
     :           'sla_PREBN', '(2,3)', STATUS )
      CALL VVD ( RMATP(3,1),  4.858341176745641D-3, 1D-12,
     :           'sla_PREBN', '(3,1)', STATUS )
      CALL VVD ( RMATP(3,2), -2.714330927085065D-5, 1D-12,
     :           'sla_PREBN', '(3,2)', STATUS )
      CALL VVD ( RMATP(3,3),  9.999881978224798D-1, 1D-12,
     :           'sla_PREBN', '(3,3)', STATUS )

      END

      SUBROUTINE T_PREC ( STATUS )
*+
*  - - - - - - -
*   T _ P R E C
*  - - - - - - -
*
*  Test sla_PREC and sla_PRECL routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PREC, sla_PRECL, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RMATP(3,3)

      CALL sla_PREC ( 1925D0, 1975D0, RMATP )

      CALL VVD ( RMATP(1,1),  9.999257249850045D-1, 1D-12,
     :           'sla_PREC', '(1,1)', STATUS )
      CALL VVD ( RMATP(1,2), -1.117719859160180D-2, 1D-12,
     :           'sla_PREC', '(1,2)', STATUS )
      CALL VVD ( RMATP(1,3), -4.859500474027002D-3, 1D-12,
     :           'sla_PREC', '(1,3)', STATUS )
      CALL VVD ( RMATP(2,1),  1.117719858025860D-2, 1D-12,
     :           'sla_PREC', '(2,1)', STATUS )
      CALL VVD ( RMATP(2,2),  9.999375327960091D-1, 1D-12,
     :           'sla_PREC', '(2,2)', STATUS )
      CALL VVD ( RMATP(2,3), -2.716114374174549D-5, 1D-12,
     :           'sla_PREC', '(2,3)', STATUS )
      CALL VVD ( RMATP(3,1),  4.859500500117173D-3, 1D-12,
     :           'sla_PREC', '(3,1)', STATUS )
      CALL VVD ( RMATP(3,2), -2.715647545167383D-5, 1D-12,
     :           'sla_PREC', '(3,2)', STATUS )
      CALL VVD ( RMATP(3,3),  9.999881921889954D-1, 1D-12,
     :           'sla_PREC', '(3,3)', STATUS )

      CALL sla_PRECL ( 1925D0, 1975D0, RMATP )

      CALL VVD ( RMATP(1,1),  9.999257331781050D-1, 1D-12,
     :           'sla_PREC', '(1,1)', STATUS )
      CALL VVD ( RMATP(1,2), -1.117658038434041D-2, 1D-12,
     :           'sla_PREC', '(1,2)', STATUS )
      CALL VVD ( RMATP(1,3), -4.859236477249598D-3, 1D-12,
     :           'sla_PREC', '(1,3)', STATUS )
      CALL VVD ( RMATP(2,1),  1.117658037299592D-2, 1D-12,
     :           'sla_PREC', '(2,1)', STATUS )
      CALL VVD ( RMATP(2,2),  9.999375397061558D-1, 1D-12,
     :           'sla_PREC', '(2,2)', STATUS )
      CALL VVD ( RMATP(2,3), -2.715816653174189D-5, 1D-12,
     :           'sla_PREC', '(2,3)', STATUS )
      CALL VVD ( RMATP(3,1),  4.859236503342703D-3, 1D-12,
     :           'sla_PREC', '(3,1)', STATUS )
      CALL VVD ( RMATP(3,2), -2.715349745834860D-5, 1D-12,
     :           'sla_PREC', '(3,2)', STATUS )
      CALL VVD ( RMATP(3,3),  9.999881934719490D-1, 1D-12,
     :           'sla_PREC', '(3,3)', STATUS )

      END

      SUBROUTINE T_PRECES ( STATUS )
*+
*  - - - - - - - - -
*   T _ P R E C E S
*  - - - - - - - - -
*
*  Test sla_PRECES routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PRECES, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RA, DC

      RA = 6.28D0
      DC = -1.123D0
      CALL sla_PRECES ( 'FK4', 1925D0, 1950D0, RA, DC )
      CALL VVD ( RA,  0.002403604864728447D0, 1D-12, 'sla_PRECES',
     :           'R', STATUS )
      CALL VVD ( DC, -1.120570643322045D0, 1D-12, 'sla_PRECES',
     :           'D', STATUS )

      RA = 0.0123D0
      DC = 1.0987D0
      CALL sla_PRECES ( 'FK5', 2050D0, 1990D0, RA, DC )
      CALL VVD ( RA, 6.282003602708382D0, 1D-12, 'sla_PRECES',
     :           'R', STATUS )
      CALL VVD ( DC, 1.092870326188383D0, 1D-12, 'sla_PRECES',
     :           'D', STATUS )

      END

      SUBROUTINE T_PRENUT ( STATUS )
*+
*  - - - - - - - - -
*   P R E N U T
*  - - - - - - - - -
*
*  Test sla_PRENUT routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PRENUT, VVD.
*
*  Last revision:   16 November 2001
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION RMATPN(3,3)

      CALL sla_PRENUT ( 1985D0, 50123.4567D0, RMATPN )

      CALL VVD ( RMATPN(1,1),  9.999962358680738D-1, 1D-12,
     :           'sla_PRENUT', '(1,1)', STATUS )
      CALL VVD ( RMATPN(1,2), -2.516417057665452D-3, 1D-12,
     :           'sla_PRENUT', '(1,2)', STATUS )
      CALL VVD ( RMATPN(1,3), -1.093569785342370D-3, 1D-12,
     :           'sla_PRENUT', '(1,3)', STATUS )
      CALL VVD ( RMATPN(2,1),  2.516462370370876D-3, 1D-12,
     :           'sla_PRENUT', '(2,1)', STATUS )
      CALL VVD ( RMATPN(2,2),  9.999968329010883D-1, 1D-12,
     :           'sla_PRENUT', '(2,2)', STATUS )
      CALL VVD ( RMATPN(2,3),  4.006159587358310D-5, 1D-12,
     :           'sla_PRENUT', '(2,3)', STATUS )
      CALL VVD ( RMATPN(3,1),  1.093465510215479D-3, 1D-12,
     :           'sla_PRENUT', '(3,1)', STATUS )
      CALL VVD ( RMATPN(3,2), -4.281337229063151D-5, 1D-12,
     :           'sla_PRENUT', '(3,2)', STATUS )
      CALL VVD ( RMATPN(3,3),  9.999994012499173D-1, 1D-12,
     :           'sla_PRENUT', '(3,3)', STATUS )

      END

      SUBROUTINE T_PVOBS ( STATUS )
*+
*  - - - - - - - -
*   T _ P V O B S
*  - - - - - - - -
*
*  Test sla_PVOBS routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_PVOBS, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION PV(6)

      CALL sla_PVOBS ( 0.5123D0, 3001D0, -0.567D0, PV )

      CALL VVD ( PV(1), 0.3138647803054939D-4, 1D-16, 'sla_PVOBS',
     :           '(1)', STATUS )
      CALL VVD ( PV(2),-0.1998515596527082D-4, 1D-16, 'sla_PVOBS',
     :           '(2)', STATUS )
      CALL VVD ( PV(3), 0.2078572043443275D-4, 1D-16, 'sla_PVOBS',
     :           '(3)', STATUS )
      CALL VVD ( PV(4), 0.1457340726851264D-8, 1D-20, 'sla_PVOBS',
     :           '(4)', STATUS )
      CALL VVD ( PV(5), 0.2288738340888011D-8, 1D-20, 'sla_PVOBS',
     :           '(5)', STATUS )
      CALL VVD ( PV(6), 0D0, 0D0, 'sla_PVOBS',
     :           '(6)', STATUS )

      END

      SUBROUTINE T_RANGE ( STATUS )
*+
*  - - - - - - - -
*   T _ R A N G E
*  - - - - - - - -
*
*  Test sla_RANGE, sla_DRANGE routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_RANGE, VVD, sla_DRANGE.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL sla_RANGE
      DOUBLE PRECISION sla_DRANGE


      CALL VVD ( DBLE( sla_RANGE ( -4.0 ) ), 2.283185307179586D0,
     :           1D-6, 'sla_RANGE', ' ', STATUS )
      CALL VVD ( sla_DRANGE ( -4D0 ), 2.283185307179586D0,
     :           1D-12, 'sla_DRANGE', ' ', STATUS )

      END

      SUBROUTINE T_RANORM ( STATUS )
*+
*  - - - - - - - - -
*   T _ R A N O R M
*  - - - - - - - - -
*
*  Test sla_RANORM, sla_DRANRM routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_RANORM, VVD, sla_DRANRM.
*
*  Last revision:   22 October 2006
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL sla_RANORM
      DOUBLE PRECISION sla_DRANRM


      CALL VVD ( DBLE( sla_RANORM ( -0.1E0 ) ), 6.183185307179587D0,
     :           1D-5, 'sla_RANORM', '1', STATUS )
      CALL VVD ( sla_DRANRM ( -0.1D0 ), 6.183185307179587D0,
     :           1D-12, 'sla_DRANRM', '2', STATUS )

      END

      SUBROUTINE T_RCC ( STATUS )
*+
*  - - - - - -
*   T _ R C C
*  - - - - - -
*
*  Test sla_RCC routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_RCC, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_RCC


      CALL VVD ( sla_RCC ( 48939.123D0, 0.76543D0, 5.0123D0,
     :                     5525.242D0, 3190D0 ),
     :           -1.280131613589158D-3, 1D-15, 'sla_RCC', ' ', STATUS )

      END

      SUBROUTINE T_REF ( STATUS )
*+
*  - - - - - -
*   T _ R E F
*  - - - - - -
*
*  Test sla_REFRO, sla_REFCO, sla_ATMDSP, sla_REFV, sla_REFZ routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_REFRO, VVD, sla_REFCO, sla_REFCOQ, sla_ATMDSP,
*           sla_DCS2C, sla_REFV, sla_REFZ.
*
*  Last revision:   17 January 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION REF, REFA, REFB, REFA2, REFB2, VU(3), VR(3), ZR

      CALL sla_REFRO ( 1.4D0, 3456.7D0, 280D0, 678.9D0, 0.9D0, 0.55D0,
     :                 -0.3D0, 0.006D0, 1D-9, REF )
      CALL VVD ( REF, 0.00106715763018568D0, 1D-12, 'sla_REFRO',
     :           'O', STATUS )

      CALL sla_REFRO ( 1.4D0, 3456.7D0, 280D0, 678.9D0, 0.9D0, 1000D0,
     :                 -0.3D0, 0.006D0, 1D-9, REF )
      CALL VVD ( REF, 0.001296416185295403D0, 1D-12, 'sla_REFRO',
     :           'R', STATUS )

      CALL sla_REFCOQ ( 275.9D0, 709.3D0, 0.9D0, 101D0, REFA, REFB )
      CALL VVD ( REFA, 2.324736903790639D-4, 1D-12, 'sla_REFCOQ',
     :           'A/R', STATUS )
      CALL VVD ( REFB, -2.442884551059D-7, 1D-15, 'sla_REFCOQ',
     :           'B/R', STATUS )

      CALL sla_REFCO ( 2111.1D0, 275.9D0, 709.3D0, 0.9D0, 101D0,
     :                 -1.03D0, 0.0067D0, 1D-12, REFA, REFB )
      CALL VVD ( REFA, 2.324673985217244D-4, 1D-12, 'sla_REFCO',
     :           'A/R', STATUS )
      CALL VVD ( REFB, -2.265040682496D-7, 1D-15, 'sla_REFCO',
     :           'B/R', STATUS )

      CALL sla_REFCOQ ( 275.9D0, 709.3D0, 0.9D0, 0.77D0, REFA, REFB )
      CALL VVD ( REFA, 2.007406521596588D-4, 1D-12, 'sla_REFCOQ',
     :           'A', STATUS )
      CALL VVD ( REFB, -2.264210092590D-7, 1D-15, 'sla_REFCOQ',
     :           'B', STATUS )

      CALL sla_REFCO ( 2111.1D0, 275.9D0, 709.3D0, 0.9D0, 0.77D0,
     :                 -1.03D0, 0.0067D0, 1D-12, REFA, REFB )
      CALL VVD ( REFA, 2.007202720084551D-4, 1D-12, 'sla_REFCO',
     :           'A', STATUS )
      CALL VVD ( REFB, -2.223037748876D-7, 1D-15, 'sla_REFCO',
     :           'B', STATUS )

      CALL sla_ATMDSP ( 275.9D0, 709.3D0, 0.9D0, 0.77D0,
     :                  REFA, REFB, 0.5D0, REFA2, REFB2 )
      CALL VVD ( REFA2, 2.034523658888048D-4, 1D-12, 'sla_ATMDSP',
     :           'A', STATUS )
      CALL VVD ( REFB2, -2.250855362179D-7, 1D-15, 'sla_ATMDSP',
     :           'B', STATUS )

      CALL sla_DCS2C ( 0.345D0, 0.456D0, VU )
      CALL sla_REFV ( VU, REFA, REFB, VR )
      CALL VVD ( VR(1), 0.8447487047790478D0, 1D-12, 'sla_REFV',
     :           'X1', STATUS )
      CALL VVD ( VR(2), 0.3035794890562339D0, 1D-12, 'sla_REFV',
     :           'Y1', STATUS )
      CALL VVD ( VR(3), 0.4407256738589851D0, 1D-12, 'sla_REFV',
     :           'Z1', STATUS )

      CALL sla_DCS2C ( 3.7D0, 0.03D0, VU )
      CALL sla_REFV ( VU, REFA, REFB, VR )
      CALL VVD ( VR(1), -0.8476187691681673D0, 1D-12, 'sla_REFV',
     :           'X2', STATUS )
      CALL VVD ( VR(2), -0.5295354802804889D0, 1D-12, 'sla_REFV',
     :           'Y2', STATUS )
      CALL VVD ( VR(3), 0.0322914582168426D0, 1D-12, 'sla_REFV',
     :           'Z2', STATUS )

      CALL sla_REFZ ( 0.567D0, REFA, REFB, ZR )
      CALL VVD ( ZR, 0.566872285910534D0, 1D-12, 'sla_REFZ',
     :           'hi el', STATUS )

      CALL sla_REFZ ( 1.55D0, REFA, REFB, ZR )
      CALL VVD ( ZR, 1.545697350690958D0, 1D-12, 'sla_REFZ',
     :           'lo el', STATUS )

      END

      SUBROUTINE T_RV ( STATUS )
*+
*  - - - - -
*   T _ R V
*  - - - - -
*
*  Test sla_RVEROT, sla_RVGALC, sla_RVLG, sla_RVLSRD, sla_RVLSRK routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called: VVD, sla_RVEROT, sla_RVGALC, sla_RVLG, sla_RVLSRD, sla_RVLSRK.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      REAL sla_RVEROT, sla_RVGALC, sla_RVLG, sla_RVLSRD, sla_RVLSRK


      CALL VVD ( DBLE( sla_RVEROT ( -0.777E0, 5.67E0, -0.3E0,
     :           3.19E0 ) ), -0.1948098355075913D0, 1D-6,
     :           'sla_RVEROT', ' ', STATUS )
      CALL VVD ( DBLE( sla_RVGALC ( 1.11E0, -0.99E0 ) ),
     :           158.9630759840254D0, 1D-3, 'sla_RVGALC', ' ', STATUS )
      CALL VVD ( DBLE( sla_RVLG ( 3.97E0, 1.09E0 ) ),
     :           -197.818762175363D0, 1D-3, 'sla_RVLG', ' ', STATUS )
      CALL VVD ( DBLE( sla_RVLSRD ( 6.01E0, 0.1E0 ) ),
     :           -4.082811335150567D0, 1D-4, 'sla_RVLSRD', ' ', STATUS )
      CALL VVD ( DBLE( sla_RVLSRK ( 6.01E0, 0.1E0 ) ),
     :           -5.925180579830265D0, 1D-4, 'sla_RVLSRK', ' ', STATUS )

      END

      SUBROUTINE T_SEP ( STATUS )
*+
*  - - - - - - -
*   T _ S E P
*  - - - - - - -
*
*  Test sla_DSEP, sla_DSEPV, sla_SEP, sla_SEPV routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_DSEP, sla_SEP, VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I
      REAL sla_SEP, sla_SEPV
      REAL R1(3), R2(3), AR1, BR1, AR2, BR2
      DOUBLE PRECISION sla_DSEP, sla_DSEPV
      DOUBLE PRECISION D1(3), D2(3), AD1, BD1, AD2, BD2


      R1(1) = 1.0
      R1(2) = 0.1
      R1(3) = 0.2
      R2(1) = -3.0
      R2(2) = 1E-3
      R2(3) = 0.2

      DO I = 1, 3
         D1(I) = DBLE( R1(I) )
         D2(I) = DBLE( R2(I) )
      END DO

      CALL sla_DCC2S ( D1, AD1, BD1 )
      CALL sla_DCC2S ( D2, AD2, BD2 )

      AR1 = SNGL( AD1 )
      BR1 = SNGL( BD1 )
      AR2 = SNGL( AD2 )
      BR2 = SNGL( BD2 )

      CALL VVD ( sla_DSEP ( AD1, BD1, AD2, BD2 ),
     :           2.8603919190246608D0, 1D-7, 'sla_DSEP', ' ', STATUS )
      CALL VVD ( DBLE( sla_SEP ( AR1, BR1, AR2, BR2 ) ),
     :           2.8603919190246608D0, 1D-4, 'sla_SEP', ' ', STATUS )
      CALL VVD ( sla_DSEPV ( D1, D2 ),
     :           2.8603919190246608D0, 1D-7, 'sla_DSEPV', ' ', STATUS )
      CALL VVD ( DBLE( sla_SEPV ( R1, R2 ) ),
     :           2.8603919190246608D0, 1D-4, 'sla_SEPV', ' ', STATUS )

      END

      SUBROUTINE T_SMAT ( STATUS )
*+
*  - - - - - - -
*   T _ S M A T
*  - - - - - - -
*
*  Test sla_SMAT routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_SMAT, VVD, VIV.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J, IW(3)
      REAL A(3,3)
      REAL V(3)
      REAL D

      DATA A/2.22E0,     1.6578E0,     1.380522E0,
     :       1.6578E0,   1.380522E0,   1.22548578E0,
     :       1.380522E0, 1.22548578E0, 1.1356276122E0/
      DATA V/2.28625E0,  1.7128825E0,  1.429432225E0/


      CALL sla_SMAT ( 3, A, V, D, J, IW )

      CALL VVD ( DBLE( A(1,1) ), 18.02550629769198D0,
     :           1D-2, 'sla_SMAT', 'A(0,0)', STATUS )
      CALL VVD ( DBLE( A(1,2) ), -52.16386644917481D0,
     :           1D-2, 'sla_SMAT', 'A(0,1)', STATUS )
      CALL VVD ( DBLE( A(1,3) ), 34.37875949717994D0,
     :           1D-2, 'sla_SMAT', 'A(0,2)', STATUS )
      CALL VVD ( DBLE( A(2,1) ), -52.16386644917477D0,
     :           1D-2, 'sla_SMAT', 'A(1,0)', STATUS )
      CALL VVD ( DBLE( A(2,2) ), 168.1778099099869D0,
     :           1D-1, 'sla_SMAT', 'A(1,1)', STATUS )
      CALL VVD ( DBLE( A(2,3) ), -118.0722869694278D0,
     :           1D-2, 'sla_SMAT', 'A(1,2)', STATUS )
      CALL VVD ( DBLE( A(3,1) ), 34.37875949717988D0,
     :           1D-2, 'sla_SMAT', 'A(2,0)', STATUS )
      CALL VVD ( DBLE( A(3,2) ), -118.07228696942770D0,
     :           1D-2, 'sla_SMAT', 'A(2,1)', STATUS )
      CALL VVD ( DBLE( A(3,3) ), 86.50307003740468D0,
     :           1D-2, 'sla_SMAT', 'A(2,2)', STATUS )
      CALL VVD ( DBLE( V(1) ), 1.002346480763383D0,
     :           1D-4, 'sla_SMAT', 'V(1)', STATUS )
      CALL VVD ( DBLE( V(2) ), 0.0328559401697292D0,
     :           1D-4, 'sla_SMAT', 'V(2)', STATUS )
      CALL VVD ( DBLE( V(3) ), 0.004760688414898454D0,
     :           1D-4, 'sla_SMAT', 'V(3)', STATUS )
      CALL VVD ( DBLE( D ), 0.003658344147359863D0,
     :           1D-4, 'sla_SMAT', 'D', STATUS )
      CALL VIV ( J, 0, 'sla_SMAT', 'J', STATUS )

      END

      SUBROUTINE T_SUPGAL ( STATUS )
*+
*  - - - - - - - - -
*   T _ S U P G A L
*  - - - - - - - - -
*
*  Test sla_SUPGAL routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_SUPGAL, VVD.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION DL, DB

      CALL sla_SUPGAL ( 6.1D0, -1.4D0, Dl, DB )

      CALL VVD ( DL, 3.798775860769474D0, 1D-12, 'sla_SUPGAL',
     :           'DL', STATUS )
      CALL VVD ( DB, -0.1397070490669407D0, 1D-12, 'sla_SUPGAL',
     :           'DB', STATUS )

      END

      SUBROUTINE T_SVD ( STATUS )
*+
*  - - - - - -
*   T _ S V D
*  - - - - - -
*
*  Test sla_SVD, sla_SVDSOL, sla_SVDCOV routines.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  sla_SVD, VVD, sla_SVDSOL, sla_SVDCOV.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER M, N
      INTEGER I, J
      INTEGER MP, NP, NC

      PARAMETER (MP = 10)
      PARAMETER (NP = 6)
      PARAMETER (NC = 7)

      DOUBLE PRECISION A(MP,NP), W(NP), V(NP,NP), WORK(NP),
     :                  B(MP), X(NP), C(NC,NC)
      DOUBLE PRECISION VAL

      M = 5
      N = 4

      DO I = 1, M
         VAL = DFLOAT( ( I ) ) / 2D0
         B(I) = 23D0 - 3D0 * VAL - 11D0 * DSIN ( VAL ) +
     :          13D0 * DCOS ( VAL )
         A(I,1) = 1D0
         A(I,2) = VAL
         A(I,3) = DSIN ( VAL )
         A(I,4) = DCOS ( VAL )
      END DO

      CALL sla_SVD ( M, N, MP, NP, A, W, V, WORK, J )

*  Allow U and V to have reversed signs.
      IF (A(1,1) .GT. 0D0) THEN
         DO I = 1, M
            DO J = 1, N
               A(I,J) = - A(I,J)
               V(I,J) = - V(I,J)
            END DO
         END DO
      END IF

      CALL VVD ( A(1,1), -0.21532492989299D0, 1D-12, 'sla_SVD',
     :           'A(1,1)', STATUS )
      CALL VVD ( A(1,2),  0.67675050651267D0, 1D-12, 'sla_SVD',
     :           'A(1,2)', STATUS )
      CALL VVD ( A(1,3), -0.37267876361644D0, 1D-12, 'sla_SVD',
     :           'A(1,3)', STATUS )
      CALL VVD ( A(1,4),  0.58330405917160D0, 1D-12, 'sla_SVD',
     :           'A(1,4)', STATUS )
      CALL VVD ( A(2,1), -0.33693420368121D0, 1D-12, 'sla_SVD',
     :           'A(2,1)', STATUS )
      CALL VVD ( A(2,2),  0.48011695963936D0, 1D-12, 'sla_SVD',
     :           'A(2,2)', STATUS )
      CALL VVD ( A(2,3),  0.62656568539705D0, 1D-12, 'sla_SVD',
     :           'A(2,3)', STATUS )
      CALL VVD ( A(2,4), -0.17479918328198D0, 1D-12, 'sla_SVD',
     :           'A(2,4)', STATUS )
      CALL VVD ( A(3,1), -0.44396825906047D0, 1D-12, 'sla_SVD',
     :           'A(3,1)', STATUS )
      CALL VVD ( A(3,2),  0.18255923809825D0, 1D-12, 'sla_SVD',
     :           'A(3,2)', STATUS )
      CALL VVD ( A(3,3),  0.02228154115994D0, 1D-12, 'sla_SVD',
     :           'A(3,3)', STATUS )
      CALL VVD ( A(3,4), -0.51743308030238D0, 1D-12, 'sla_SVD',
     :           'A(3,4)', STATUS )
      CALL VVD ( A(4,1), -0.53172583816951D0, 1D-12, 'sla_SVD',
     :           'A(4,1)', STATUS )
      CALL VVD ( A(4,2), -0.16537863535943D0, 1D-12, 'sla_SVD',
     :           'A(4,2)', STATUS )
      CALL VVD ( A(4,3), -0.61134201569990D0, 1D-12, 'sla_SVD',
     :           'A(4,3)', STATUS )
      CALL VVD ( A(4,4), -0.28871221824912D0, 1D-12, 'sla_SVD',
     :           'A(4,4)', STATUS )
      CALL VVD ( A(5,1), -0.60022523682867D0, 1D-12, 'sla_SVD',
     :           'A(5,1)', STATUS )
      CALL VVD ( A(5,2), -0.50081781972404D0, 1D-12, 'sla_SVD',
     :           'A(5,2)', STATUS )
      CALL VVD ( A(5,3),  0.30706750690326D0, 1D-12, 'sla_SVD',
     :           'A(5,3)', STATUS )
      CALL VVD ( A(5,4),  0.52736124480318D0, 1D-12, 'sla_SVD',
     :           'A(5,4)', STATUS )

      CALL VVD ( W(1), 4.57362714220621D0, 1D-12, 'sla_SVD',
     :           'W(1)', STATUS )
      CALL VVD ( W(2), 1.64056393111226D0, 1D-12, 'sla_SVD',
     :           'W(2)', STATUS )
      CALL VVD ( W(3), 0.03999179717447D0, 1D-12, 'sla_SVD',
     :           'W(3)', STATUS )
      CALL VVD ( W(4), 0.37267332634218D0, 1D-12, 'sla_SVD',
     :           'W(4)', STATUS )

      CALL VVD ( V(1,1), -0.46531525230679D0, 1D-12, 'sla_SVD',
     :           'V(1,1)', STATUS )
      CALL VVD ( V(1,2),  0.41036514115630D0, 1D-12, 'sla_SVD',
     :           'V(1,2)', STATUS )
      CALL VVD ( V(1,3), -0.70279526907678D0, 1D-12, 'sla_SVD',
     :           'V(1,3)', STATUS )
      CALL VVD ( V(1,4),  0.34808185338758D0, 1D-12, 'sla_SVD',
     :           'V(1,4)', STATUS )
      CALL VVD ( V(2,1), -0.80342444002914D0, 1D-12, 'sla_SVD',
     :           'V(2,1)', STATUS )
      CALL VVD ( V(2,2), -0.29896472833787D0, 1D-12, 'sla_SVD',
     :           'V(2,2)', STATUS )
      CALL VVD ( V(2,3),  0.46592932810178D0, 1D-12, 'sla_SVD',
     :           'V(2,3)', STATUS )
      CALL VVD ( V(2,4),  0.21917828721921D0, 1D-12, 'sla_SVD',
     :           'V(2,4)', STATUS )
      CALL VVD ( V(3,1), -0.36564497020801D0, 1D-12, 'sla_SVD',
     :           'V(3,1)', STATUS )
      CALL VVD ( V(3,2),  0.28066812941896D0, 1D-12, 'sla_SVD',
     :           'V(3,2)', STATUS )
      CALL VVD ( V(3,3), -0.03324480702665D0, 1D-12, 'sla_SVD',
     :           'V(3,3)', STATUS )
      CALL VVD ( V(3,4), -0.88680546891402D0, 1D-12, 'sla_SVD',
     :           'V(3,4)', STATUS )
      CALL VVD ( V(4,1),  0.06553350971918D0, 1D-12, 'sla_SVD',
     :           'V(4,1)', STATUS )
      CALL VVD ( V(4,2),  0.81452191085452D0, 1D-12, 'sla_SVD',
     :           'V(4,2)', STATUS )
      CALL VVD ( V(4,3),  0.53654771808636D0, 1D-12, 'sla_SVD',
     :           'V(4,3)', STATUS )
      CALL VVD ( V(4,4),  0.21065602782287D0, 1D-12, 'sla_SVD',
     :           'V(4,4)', STATUS )

      CALL sla_SVDSOL ( M, N, MP, NP, B, A, W, V, WORK, X )

      CALL VVD ( X(1),  23D0, 1D-12, 'sla_SVDSOL', 'X(1)', STATUS )
      CALL VVD ( X(2),  -3D0, 1D-12, 'sla_SVDSOL', 'X(2)', STATUS )
      CALL VVD ( X(3), -11D0, 1D-12, 'sla_SVDSOL', 'X(3)', STATUS )
      CALL VVD ( X(4),  13D0, 1D-12, 'sla_SVDSOL', 'X(4)', STATUS )

      CALL sla_SVDCOV ( N, NP, NC, W, V, WORK, C )

      CALL VVD ( C(1,1),  309.77269378273270D0, 1D-10,
     :           'sla_SVDCOV', 'C(1,1)', STATUS )
      CALL VVD ( C(1,2), -204.22043941662150D0, 1D-10,
     :           'sla_SVDCOV', 'C(1,2)', STATUS )
      CALL VVD ( C(1,3),   12.43704316907477D0, 1D-10,
     :           'sla_SVDCOV', 'C(1,3)', STATUS )
      CALL VVD ( C(1,4), -235.12299986206710D0, 1D-10,
     :           'sla_SVDCOV', 'C(1,4)', STATUS )
      CALL VVD ( C(2,1), -204.22043941662150D0, 1D-10,
     :           'sla_SVDCOV', 'C(2,1)', STATUS )
      CALL VVD ( C(2,2),  136.14695961108110D0, 1D-10,
     :           'sla_SVDCOV', 'C(2,2)', STATUS )
      CALL VVD ( C(2,3),  -11.10167446246327D0, 1D-10,
     :           'sla_SVDCOV', 'C(2,3)', STATUS )
      CALL VVD ( C(2,4),  156.54937371198730D0, 1D-10,
     :           'sla_SVDCOV', 'C(2,4)', STATUS )
      CALL VVD ( C(3,1),   12.43704316907477D0, 1D-10,
     :           'sla_SVDCOV', 'C(3,1)', STATUS )
      CALL VVD ( C(3,2),  -11.10167446246327D0, 1D-10,
     :           'sla_SVDCOV', 'C(3,2)', STATUS )
      CALL VVD ( C(3,3),    6.38909830090602D0, 1D-10,
     :           'sla_SVDCOV', 'C(3,3)', STATUS )
      CALL VVD ( C(3,4),  -12.41424302586736D0, 1D-10,
     :           'sla_SVDCOV', 'C(3,4)', STATUS )
      CALL VVD ( C(4,1), -235.12299986206710D0, 1D-10,
     :           'sla_SVDCOV', 'C(4,1)', STATUS )
      CALL VVD ( C(4,2),  156.54937371198730D0, 1D-10,
     :           'sla_SVDCOV', 'C(4,2)', STATUS )
      CALL VVD ( C(4,3),  -12.41424302586736D0, 1D-10,
     :           'sla_SVDCOV', 'C(4,3)', STATUS )
      CALL VVD ( C(4,4),  180.56719842359560D0, 1D-10,
     :           'sla_SVDCOV', 'C(4,4)', STATUS )

      END

      SUBROUTINE T_TP ( STATUS )
*+
*  - - - - - -
*   T _ T P
*  - - - - - -
*
*  Test spherical tangent-planD-projection routines:
*
*       sla_S2TP      sla_DS2TP     sla_DTPS2C
*       sla_TP2S      sla_DTP2S     sla_TPS2C
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  all the above, plus VVD and VIV.
*
*  Last revision:   10 July 2000
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      REAL R0, D0, R1, D1, X, Y, R2, D2, R01, D01, R02, D02
      DOUBLE PRECISION DR0, DD0, DR1, DD1, DX, DY, DR2, DD2, DR01,
     :                 DD01, DR02, DD02

      R0 = 3.1E0
      D0 = -0.9E0
      R1 = R0 + 0.2E0
      D1 = D0 - 0.1E0
      CALL sla_S2TP ( R1, D1, R0, D0, X, Y, J )
      CALL VVD ( DBLE( X ), 0.1086112301590404D0, 1D-6, 'sla_S2TP',
     :           'X', STATUS )
      CALL VVD ( DBLE( Y ), -0.1095506200711452D0, 1D-6, 'sla_S2TP',
     :           'Y', STATUS )
      CALL VIV ( J, 0, 'sla_S2TP', 'J', STATUS )
      CALL sla_TP2S ( X, Y, R0, D0, R2, D2 )
      CALL VVD ( DBLE( ( R2 - R1 ) ), 0D0, 1D-6, 'sla_TP2S',
     :           'R', STATUS )
      CALL VVD ( DBLE( ( D2 - D1 ) ), 0D0, 1D-6, 'sla_TP2S',
     :           'D', STATUS )
      CALL sla_TPS2C ( X, Y, R2, D2, R01, D01, R02, D02, J )
      CALL VVD ( DBLE( R01 ),  3.1D0, 1D-6, 'sla_TPS2C',
     :           'R1', STATUS )
      CALL VVD ( DBLE( D01 ), -0.9D0, 1D-6, 'sla_TPS2C',
     :           'D1', STATUS )
      CALL VVD ( DBLE( R02 ), 0.3584073464102072D0, 1D-6, 'sla_TPS2C',
     :           'R2', STATUS )
      CALL VVD ( DBLE( D02 ), -2.023361658234722D0, 1D-6, 'sla_TPS2C',
     :           'D2', STATUS )
      CALL VIV ( J, 1, 'sla_TPS2C', 'N', STATUS )

      DR0 = 3.1D0
      DD0 = -0.9D0
      DR1 = DR0 + 0.2D0
      DD1 = DD0 - 0.1D0
      CALL sla_DS2TP ( DR1, DD1, DR0, DD0, DX, DY, J )
      CALL VVD ( DX, 0.1086112301590404D0, 1D-12, 'sla_DS2TP',
     :           'X', STATUS )
      CALL VVD ( DY, -0.1095506200711452D0, 1D-12, 'sla_DS2TP',
     :           'Y', STATUS )
      CALL VIV ( J, 0, 'sla_DS2TP', 'J', STATUS )
      CALL sla_DTP2S ( DX, DY, DR0, DD0, DR2, DD2 )
      CALL VVD ( DR2 - DR1, 0D0, 1D-12, 'sla_DTP2S', 'R', STATUS )
      CALL VVD ( DD2 - DD1, 0D0, 1D-12, 'sla_DTP2S', 'D', STATUS )
      CALL sla_DTPS2C ( DX, DY, DR2, DD2, DR01, DD01, DR02, DD02, J )
      CALL VVD ( DR01,  3.1D0, 1D-12, 'sla_DTPS2C', 'R1', STATUS )
      CALL VVD ( DD01, -0.9D0, 1D-12, 'sla_DTPS2C', 'D1', STATUS )
      CALL VVD ( DR02, 0.3584073464102072D0, 1D-12, 'sla_DTPS2C',
     :           'R2', STATUS )
      CALL VVD ( DD02, -2.023361658234722D0, 1D-12, 'sla_DTPS2C',
     :           'D2', STATUS )
      CALL VIV ( J, 1, 'sla_DTPS2C', 'N', STATUS )

      END

      SUBROUTINE T_TPV ( STATUS )
*+
*  - - - - - -
*   T _ T P V
*  - - - - - -
*
*  Test Cartesian tangent-planD-projection routines:
*
*       sla_TP2V      sla_V2TP      sla_TPV2C
*       sla_DTP2V     sla_DV2TP     sla_DTPV2C
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  all the above, plus VVD and VIV.
*
*  Last revision:   21 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER J
      REAL RXI, RETA, RV(3), RV0(3), RTXI, RTETA, RTV(3),
     :     RTV01(3), RTV02(3)
      DOUBLE PRECISION XI, ETA, X, Y, Z, V(3), V0(3), TXI, TETA,
     :                 TV(3), TV01(3), TV02(3)

      XI = -0.1D0
      ETA = 0.055D0
      RXI = SNGL( XI )
      RETA = SNGL( ETA )

      X = -0.7D0
      Y = -0.13D0
      Z = DSQRT ( 1D0 - X * X - Y * Y )
      RV(1) = SNGL( X )
      RV(2) = SNGL( Y )
      RV(3) = SNGL( Z )
      V(1) = X
      V(2) = Y
      V(3) = Z

      X = -0.72D0
      Y = -0.16D0
      Z = DSQRT ( 1D0 - X * X - Y * Y )
      RV0(1) = SNGL( X )
      RV0(2) = SNGL( Y )
      RV0(3) = SNGL( Z )
      V0(1) = X
      V0(2) = Y
      V0(3) = Z

      CALL sla_TP2V ( RXI, RETA, RV0, RTV )
      CALL VVD ( DBLE( RTV(1) ), -0.700887428128D0, 1D-6, 'sla_TP2V',
     :           'V(1)', STATUS )
      CALL VVD ( DBLE( RTV(2) ), -0.05397407D0, 1D-6, 'sla_TP2V',
     :           'V(2)', STATUS )
      CALL VVD ( DBLE( RTV(3) ),  0.711226836562D0, 1D-6, 'sla_TP2V',
     :           'V(3)', STATUS )

      CALL sla_DTP2V ( XI, ETA, V0, TV )
      CALL VVD ( TV(1), -0.7008874281280771D0, 1D-13, 'sla_DTP2V',
     :           'V(1)', STATUS )
      CALL VVD ( TV(2), -0.05397406827952735D0, 1D-13, 'sla_DTP2V',
     :           'V(2)', STATUS )
      CALL VVD ( TV(3), 0.7112268365615617D0, 1D-13, 'sla_DTP2V',
     :           'V(3)', STATUS )

      CALL sla_V2TP ( RV, RV0, RTXI, RTETA, J)
      CALL VVD ( DBLE( RTXI ), -0.02497229197D0, 1D-6, 'sla_V2TP',
     :           'XI', STATUS )
      CALL VVD ( DBLE( RTETA ), 0.03748140764D0, 1D-6, 'sla_V2TP',
     :           'ETA', STATUS )
      CALL VIV ( J, 0, 'sla_V2TP', 'J', STATUS )

      CALL sla_DV2TP ( V, V0, TXI, TETA, J )
      CALL VVD ( TXI, -0.02497229197023852D0, 1D-13, 'sla_DV2TP',
     :           'XI', STATUS )
      CALL VVD ( TETA, 0.03748140764224765D0, 1D-13, 'sla_DV2TP',
     :           'ETA', STATUS )
      CALL VIV ( J, 0, 'sla_DV2TP', 'J', STATUS )

      CALL sla_TPV2C ( RXI, RETA, RV, RTV01, RTV02, J )
      CALL VVD ( DBLE( RTV01(1) ), -0.7074573732537283D0, 1D-6,
     :           'sla_TPV2C', 'V01(1)', STATUS )
      CALL VVD ( DBLE( RTV01(2) ), -0.2372965765309941D0, 1D-6,
     :           'sla_TPV2C', 'V01(2)', STATUS )
      CALL VVD ( DBLE( RTV01(3) ), 0.6657284730245545D0, 1D-6,
     :           'sla_TPV2C', 'V01(3)', STATUS )
      CALL VVD ( DBLE( RTV02(1) ), -0.6680480104758149D0, 1D-6,
     :           'sla_TPV2C', 'V02(1)', STATUS )
      CALL VVD ( DBLE( RTV02(2) ), -0.02915588494045333D0, 1D-6,
     :           'sla_TPV2C', 'V02(2)', STATUS )
      CALL VVD ( DBLE( RTV02(3) ), 0.7435467638774610D0, 1D-6,
     :           'sla_TPV2C', 'V02(3)', STATUS )
      CALL VIV ( J, 1, 'sla_TPV2C', 'N', STATUS )

      CALL sla_DTPV2C ( XI, ETA, V, TV01, TV02, J )
      CALL VVD ( TV01(1), -0.7074573732537283D0, 1D-13, 'sla_DTPV2C',
     :           'V01(1)', STATUS )
      CALL VVD ( TV01(2), -0.2372965765309941D0, 1D-13, 'sla_DTPV2C',
     :           'V01(2)', STATUS )
      CALL VVD ( TV01(3), 0.6657284730245545D0, 1D-13, 'sla_DTPV2C',
     :           'V01(3)', STATUS )
      CALL VVD ( TV02(1), -0.6680480104758149D0, 1D-13, 'sla_DTPV2C',
     :           'V02(1)', STATUS )
      CALL VVD ( TV02(2), -0.02915588494045333D0, 1D-13, 'sla_DTPV2C',
     :           'V02(2)', STATUS )
      CALL VVD ( TV02(3), 0.7435467638774610D0, 1D-13, 'sla_DTPV2C',
     :           'V02(3)', STATUS )
      CALL VIV ( J, 1, 'sla_DTPV2C', 'N', STATUS )

      END

      SUBROUTINE T_VECMAT ( STATUS )
*+
*  - - - - - - - - -
*   T _ V E C M A T
*  - - - - - - - - -
*
*  Test all the 3-vector and 3x3 matrix routines:
*
*          sla_AV2M     sla_DAV2M
*          sla_CC2S     sla_DCC2S
*          sla_CS2C     sla_DCS2C
*          sla_EULER    sla_DEULER
*          sla_IMXV     sla_DIMXV
*          sla_M2AV     sla_DM2AV
*          sla_MXM      sla_DMXM
*          sla_MXV      sla_DMXV
*          sla_VDV      sla_DVDV
*          sla_VN       sla_DVN
*          sla_VXV      sla_DVXV
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called:  all the above, plus VVD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      INTEGER I
      REAL sla_VDV
      REAL AV(3), RM1(3,3), RM2(3,3), RM(3,3), V1(3), V2(3),
     :   V3(3), V4(3), V5(3), VM, V6(3), V7(3)
      DOUBLE PRECISION sla_DVDV
      DOUBLE PRECISION DAV(3), DRM1(3,3), DRM2(3,3), DRM(3,3),
     :                 DV1(3), DV2(3), DV3(3), DV4(3), DV5(3),
     :                 DVM, DV6(3), DV7(3)


*  Make a rotation matrix.
      AV(1) = -0.123E0
      AV(2) = 0.0987E0
      AV(3) = 0.0654E0
      CALL sla_AV2M ( AV, RM1 )
      CALL VVD ( DBLE( RM1(1,1) ), 0.9930075842721269D0,
     :           1D-6, 'sla_AV2M', '11', STATUS )
      CALL VVD ( DBLE( RM1(1,2) ), 0.05902743090199868D0,
     :           1D-6, 'sla_AV2M', '12', STATUS )
      CALL VVD ( DBLE( RM1(1,3) ), -0.1022335560329612D0,
     :           1D-6, 'sla_AV2M', '13', STATUS )
      CALL VVD ( DBLE( RM1(2,1) ), -0.07113807138648245D0,
     :           1D-6, 'sla_AV2M', '21', STATUS )
      CALL VVD ( DBLE( RM1(2,2) ), 0.9903204657727545D0,
     :           1D-6, 'sla_AV2M', '22', STATUS )
      CALL VVD ( DBLE( RM1(2,3) ), -0.1191836812279541D0,
     :           1D-6, 'sla_AV2M', '23', STATUS )
      CALL VVD ( DBLE( RM1(3,1) ), 0.09420887631983825D0,
     :           1D-6, 'sla_AV2M', '31', STATUS )
      CALL VVD ( DBLE( RM1(3,2) ), 0.1256229973879967D0,
     :           1D-6, 'sla_AV2M', '32', STATUS )
      CALL VVD ( DBLE( RM1(3,3) ), 0.9875948309655174D0,
     :           1D-6, 'sla_AV2M', '33', STATUS )

*  Make another.
      CALL sla_EULER ( 'YZY', 2.345E0, -0.333E0, 2.222E0, RM2 )
      CALL VVD ( DBLE( RM2(1,1) ), -0.1681574770810878D0,
     :           1D-6, 'sla_EULER', '11', STATUS )
      CALL VVD ( DBLE( RM2(1,2) ), 0.1981362273264315D0,
     :           1D-6, 'sla_EULER', '12', STATUS )
      CALL VVD ( DBLE( RM2(1,3) ), 0.9656423242187410D0,
     :           1D-6, 'sla_EULER', '13', STATUS )
      CALL VVD ( DBLE( RM2(2,1) ), -0.2285369373983370D0,
     :           1D-6, 'sla_EULER', '21', STATUS )
      CALL VVD ( DBLE( RM2(2,2) ), 0.9450659587140423D0,
     :           1D-6, 'sla_EULER', '22', STATUS )
      CALL VVD ( DBLE( RM2(2,3) ), -0.2337117924378156D0,
     :           1D-6, 'sla_EULER', '23', STATUS )
      CALL VVD ( DBLE( RM2(3,1) ), -0.9589024617479674D0,
     :           1D-6, 'sla_EULER', '31', STATUS )
      CALL VVD ( DBLE( RM2(3,2) ), -0.2599853247796050D0,
     :           1D-6, 'sla_EULER', '32', STATUS )
      CALL VVD ( DBLE( RM2(3,3) ), -0.1136384607117296D0,
     :           1D-6, 'sla_EULER', '33', STATUS )

*  Combine them.
      CALL sla_MXM ( RM2, RM1, RM )
      CALL VVD ( DBLE( RM(1,1) ), -0.09010460088585805D0,
     :           1D-6, 'sla_MXM', '11', STATUS )
      CALL VVD ( DBLE( RM(1,2) ), 0.3075993402463796D0,
     :           1D-6, 'sla_MXM', '12', STATUS )
      CALL VVD ( DBLE( RM(1,3) ), 0.9472400998581048D0,
     :           1D-6, 'sla_MXM', '13', STATUS )
      CALL VVD ( DBLE( RM(2,1) ), -0.3161868071070688D0,
     :           1D-6, 'sla_MXM', '21', STATUS )
      CALL VVD ( DBLE( RM(2,2) ), 0.8930686362478707D0,
     :           1D-6, 'sla_MXM', '22', STATUS )
      CALL VVD ( DBLE( RM(2,3) ),-0.3200848543149236D0,
     :           1D-6, 'sla_MXM', '23', STATUS )
      CALL VVD ( DBLE( RM(3,1) ),-0.9444083141897035D0,
     :           1D-6, 'sla_MXM', '31', STATUS )
      CALL VVD ( DBLE( RM(3,2) ),-0.3283459407855694D0,
     :           1D-6, 'sla_MXM', '32', STATUS )
      CALL VVD ( DBLE( RM(3,3) ), 0.01678926022795169D0,
     :           1D-6, 'sla_MXM', '33', STATUS )

*  Create a vector.
      CALL sla_CS2C ( 3.0123E0, -0.999E0, V1 )
      CALL VVD ( DBLE( V1(1) ), -0.5366267667260525D0,
     :           1D-6, 'sla_CS2C', 'X', STATUS )
      CALL VVD ( DBLE( V1(2) ), 0.06977111097651444D0,
     :           1D-6, 'sla_CS2C', 'Y', STATUS )
      CALL VVD ( DBLE( V1(3) ), -0.8409302618566215D0,
     :           1D-6, 'sla_CS2C', 'Z', STATUS )

*  Rotate it using the two matrices sequentially.
      CALL sla_MXV ( RM1, V1, V2 )
      CALL sla_MXV ( RM2, V2, V3 )
      CALL VVD ( DBLE( V3(1) ), -0.7267487768696160D0,
     :           1D-6, 'sla_MXV', 'X', STATUS )
      CALL VVD ( DBLE( V3(2) ), 0.5011537352639822D0,
     :           1D-6, 'sla_MXV', 'Y', STATUS )
      CALL VVD ( DBLE( V3(3) ), 0.4697671220397141D0,
     :           1D-6, 'sla_MXV', 'Z', STATUS )

*  Derotate it using the combined matrix.
      CALL sla_IMXV ( RM, V3, V4 )
      CALL VVD ( DBLE( V4(1) ), -0.5366267667260526D0,
     :           1D-6, 'sla_IMXV', 'X', STATUS )
      CALL VVD ( DBLE( V4(2) ), 0.06977111097651445D0,
     :           1D-6, 'sla_IMXV', 'Y', STATUS )
      CALL VVD ( DBLE( V4(3) ), -0.8409302618566215D0,
     :           1D-6, 'sla_IMXV', 'Z', STATUS )

*  Convert the combined matrix into an axial vector.
      CALL sla_M2AV ( RM, V5 )
      CALL VVD ( DBLE( V5(1) ), 0.006889040510209034D0,
     :           1D-6, 'sla_M2AV', 'X', STATUS )
      CALL VVD ( DBLE( V5(2) ), -1.577473205461961D0,
     :           1D-6, 'sla_M2AV', 'Y', STATUS )
      CALL VVD ( DBLE( V5(3) ), 0.5201843672856759D0,
     :           1D-6, 'sla_M2AV', 'Z', STATUS )

*  Multiply it by a scalar and then normalize.
      DO I = 1, 3
         V5(I) = V5(I) * 1000.0
      END DO

      CALL sla_VN ( V5, V6, VM )
      CALL VVD ( DBLE( V6(1) ), 0.004147420704640065D0,
     :           1D-6, 'sla_VN', 'X', STATUS )
      CALL VVD ( DBLE( V6(2) ), -0.9496888606842218D0,
     :           1D-6, 'sla_VN', 'Y', STATUS )
      CALL VVD ( DBLE( V6(3) ), 0.3131674740355448D0,
     :           1D-6, 'sla_VN', 'Z', STATUS )
      CALL VVD ( DBLE( VM ), 1661.042127339937D0,
     :           1D-3, 'sla_VN', 'M', STATUS )

*  Dot product with the original vector.
      CALL VVD ( DBLE( sla_VDV ( V6, V1 ) ),
     :           -0.3318384698006295D0, 1D-6, 'sla_VN', ' ', STATUS )

*  Cross product with the original vector.
      CALL sla_VXV (V6, V1, V7 )
      CALL VVD ( DBLE( V7(1) ), 0.7767720597123304D0,
     :           1D-6, 'sla_VXV', 'X', STATUS )
      CALL VVD ( DBLE( V7(2) ), -0.1645663574562769D0,
     :           1D-6, 'sla_VXV', 'Y', STATUS )
      CALL VVD ( DBLE( V7(3) ), -0.5093390925544726D0,
     :           1D-6, 'sla_VXV', 'Z', STATUS )

*  Same in double precision.

      DAV(1) = -0.123D0
      DAV(2) = 0.0987D0
      DAV(3) = 0.0654D0
      CALL sla_DAV2M ( DAV, DRM1 )
      CALL VVD ( DRM1(1,1), 0.9930075842721269D0, 1D-12,
     :           'sla_DAV2M', '11', STATUS )
      CALL VVD ( DRM1(1,2), 0.05902743090199868D0, 1D-12,
     :           'sla_DAV2M', '12', STATUS )
      CALL VVD ( DRM1(1,3), -0.1022335560329612D0, 1D-12,
     :           'sla_DAV2M', '13', STATUS )
      CALL VVD ( DRM1(2,1), -0.07113807138648245D0, 1D-12,
     :           'sla_DAV2M', '21', STATUS )
      CALL VVD ( DRM1(2,2), 0.9903204657727545D0, 1D-12,
     :           'sla_DAV2M', '22', STATUS )
      CALL VVD ( DRM1(2,3), -0.1191836812279541D0, 1D-12,
     :           'sla_DAV2M', '23', STATUS )
      CALL VVD ( DRM1(3,1), 0.09420887631983825D0, 1D-12,
     :           'sla_DAV2M', '31', STATUS )
      CALL VVD ( DRM1(3,2), 0.1256229973879967D0, 1D-12,
     :           'sla_DAV2M', '32', STATUS )
      CALL VVD ( DRM1(3,3), 0.9875948309655174D0, 1D-12,
     :           'sla_DAV2M', '33', STATUS )

      CALL sla_DEULER ( 'YZY', 2.345D0, -0.333D0, 2.222D0, DRM2 )
      CALL VVD ( DRM2(1,1), -0.1681574770810878D0, 1D-12,
     :           'sla_DEULER', '11', STATUS )
      CALL VVD ( DRM2(1,2), 0.1981362273264315D0, 1D-12,
     :           'sla_DEULER', '12', STATUS )
      CALL VVD ( DRM2(1,3), 0.9656423242187410D0, 1D-12,
     :           'sla_DEULER', '13', STATUS )
      CALL VVD ( DRM2(2,1), -0.2285369373983370D0, 1D-12,
     :           'sla_DEULER', '21', STATUS )
      CALL VVD ( DRM2(2,2), 0.9450659587140423D0, 1D-12,
     :           'sla_DEULER', '22', STATUS )
      CALL VVD ( DRM2(2,3), -0.2337117924378156D0, 1D-12,
     :           'sla_DEULER', '23', STATUS )
      CALL VVD ( DRM2(3,1), -0.9589024617479674D0, 1D-12,
     :           'sla_DEULER', '31', STATUS )
      CALL VVD ( DRM2(3,2), -0.2599853247796050D0, 1D-12,
     :           'sla_DEULER', '32', STATUS )
      CALL VVD ( DRM2(3,3), -0.1136384607117296D0, 1D-12,
     :           'sla_DEULER', '33', STATUS )

      CALL sla_DMXM ( DRM2, DRM1, DRM )
      CALL VVD ( DRM(1,1), -0.09010460088585805D0, 1D-12,
     :           'sla_DMXM', '11', STATUS )
      CALL VVD ( DRM(1,2), 0.3075993402463796D0, 1D-12,
     :           'sla_DMXM', '12', STATUS )
      CALL VVD ( DRM(1,3), 0.9472400998581048D0, 1D-12,
     :           'sla_DMXM', '13', STATUS )
      CALL VVD ( DRM(2,1), -0.3161868071070688D0, 1D-12,
     :           'sla_DMXM', '21', STATUS )
      CALL VVD ( DRM(2,2), 0.8930686362478707D0, 1D-12,
     :           'sla_DMXM', '22', STATUS )
      CALL VVD ( DRM(2,3), -0.3200848543149236D0, 1D-12,
     :           'sla_DMXM', '23', STATUS )
      CALL VVD ( DRM(3,1), -0.9444083141897035D0, 1D-12,
     :           'sla_DMXM', '31', STATUS )
      CALL VVD ( DRM(3,2), -0.3283459407855694D0, 1D-12,
     :           'sla_DMXM', '32', STATUS )
      CALL VVD ( DRM(3,3), 0.01678926022795169D0, 1D-12,
     :           'sla_DMXM', '33', STATUS )

      CALL sla_DCS2C ( 3.0123D0, -0.999D0, DV1 )
      CALL VVD ( DV1(1), -0.5366267667260525D0, 1D-12,
     :           'sla_DCS2C', 'X', STATUS )
      CALL VVD ( DV1(2), 0.06977111097651444D0, 1D-12,
     :           'sla_DCS2C', 'Y', STATUS )
      CALL VVD ( DV1(3), -0.8409302618566215D0, 1D-12,
     :           'sla_DCS2C', 'Z', STATUS )

      CALL sla_DMXV ( DRM1, DV1, DV2 )
      CALL sla_DMXV ( DRM2, DV2, DV3 )
      CALL VVD ( DV3(1), -0.7267487768696160D0, 1D-12,
     :           'sla_DMXV', 'X', STATUS )
      CALL VVD ( DV3(2), 0.5011537352639822D0, 1D-12,
     :           'sla_DMXV', 'Y', STATUS )
      CALL VVD ( DV3(3), 0.4697671220397141D0, 1D-12,
     :           'sla_DMXV', 'Z', STATUS )

      CALL sla_DIMXV ( DRM, DV3, DV4 )
      CALL VVD ( DV4(1), -0.5366267667260526D0, 1D-12,
     :           'sla_DIMXV', 'X', STATUS )
      CALL VVD ( DV4(2), 0.06977111097651445D0, 1D-12,
     :           'sla_DIMXV', 'Y', STATUS )
      CALL VVD ( DV4(3), -0.8409302618566215D0, 1D-12,
     :           'sla_DIMXV', 'Z', STATUS )

      CALL sla_DM2AV ( DRM, DV5 )
      CALL VVD ( DV5(1), 0.006889040510209034D0, 1D-12,
     :           'sla_DM2AV', 'X', STATUS )
      CALL VVD ( DV5(2), -1.577473205461961D0, 1D-12,
     :           'sla_DM2AV', 'Y', STATUS )
      CALL VVD ( DV5(3), 0.5201843672856759D0, 1D-12,
     :           'sla_DM2AV', 'Z', STATUS )

      DO I = 1, 3
         DV5(I) = DV5(I) * 1000D0
      END DO

      CALL sla_DVN ( DV5, DV6, DVM )
      CALL VVD ( DV6(1), 0.004147420704640065D0, 1D-12,
     :           'sla_DVN', 'X', STATUS )
      CALL VVD ( DV6(2), -0.9496888606842218D0, 1D-12,
     :           'sla_DVN', 'Y', STATUS )
      CALL VVD ( DV6(3), 0.3131674740355448D0, 1D-12,
     :           'sla_DVN', 'Z', STATUS )
      CALL VVD ( DVM, 1661.042127339937D0, 1D-9, 'sla_DVN',
     :           'M', STATUS )

      CALL VVD ( sla_DVDV ( DV6, DV1 ), -0.3318384698006295D0,
     :           1D-12, 'sla_DVN', ' ', STATUS )

      CALL sla_DVXV (DV6, DV1, DV7 )
      CALL VVD ( DV7(1), 0.7767720597123304D0, 1D-12,
     :           'sla_DVXV', 'X', STATUS )
      CALL VVD ( DV7(2), -0.1645663574562769D0, 1D-12,
     :           'sla_DVXV', 'Y', STATUS )
      CALL VVD ( DV7(3), -0.5093390925544726D0, 1D-12,
     :           'sla_DVXV', 'Z', STATUS )

      END

      SUBROUTINE T_ZD ( STATUS )
*+
*  - - - - -
*   T _ Z D
*  - - - - -
*
*  Test sla_ZD routine.
*
*  Returned:
*     STATUS    LOGICAL     .TRUE. = success, .FALSE. = fail
*
*  Called: VVD, sla_ZD.
*
*  Last revision:   22 October 2005
*
*  Copyright CLRC/Starlink.  All rights reserved.
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

      LOGICAL STATUS

      DOUBLE PRECISION sla_ZD


      CALL VVD ( sla_ZD ( -1.023D0, -0.876D0, -0.432D0 ),
     :           0.8963914139430839D0, 1D-12, 'sla_ZD', ' ', STATUS )

      END
