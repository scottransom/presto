      SUBROUTINE sla_OBS (N, C, NAME, W, P, H)
*+
*     - - - -
*      O B S
*     - - - -
*
*  Parameters of selected groundbased observing stations
*
*  Given:
*     N       int     number specifying observing station
*
*  Either given or returned
*     C       c*(*)   identifier specifying observing station
*
*  Returned:
*     NAME    c*(*)   name of specified observing station
*     W       dp      longitude (radians, West +ve)
*     P       dp      geodetic latitude (radians, North +ve)
*     H       dp      height above sea level (metres)
*
*  Notes:
*
*     Station identifiers C may be up to 10 characters long,
*     and station names NAME may be up to 40 characters long.
*
*     C and N are alternative ways of specifying the observing
*     station.  The C option, which is the most generally useful,
*     may be selected by specifying an N value of zero or less.
*     If N is 1 or more, the parameters of the Nth station
*     in the currently supported list are interrogated, and
*     the station identifier C is returned as well as NAME, W,
*     P and H.
*
*     If the station parameters are not available, either because
*     the station identifier C is not recognized, or because an
*     N value greater than the number of stations supported is
*     given, a name of '?' is returned and C, W, P and H are left
*     in their current states.
*
*     Programs can obtain a list of all currently supported
*     stations by calling the routine repeatedly, with N=1,2,3...
*     When NAME='?' is seen, the list of stations has been
*     exhausted.
*
*     Station numbers, identifiers, names and other details are
*     subject to change and should not be hardwired into
*     application programs.
*
*     All station identifiers C are uppercase only;  lowercase
*     characters must be converted to uppercase by the calling
*     program.  The station names returned may contain both upper-
*     and lowercase.  All characters up to the first space are
*     checked;  thus an abbreviated ID will return the parameters
*     for the first station in the list which matches the
*     abbreviation supplied, and no station in the list will ever
*     contain embedded spaces.  C must not have leading spaces.
*
*     IMPORTANT -- BEWARE OF THE LONGITUDE SIGN CONVENTION.  The
*     longitude returned by sla_OBS is west-positive in accordance
*     with astronomical usage.  However, this sign convention is
*     left-handed and is the opposite of the one used by geographers;
*     elsewhere in SLALIB the preferable east-positive convention is
*     used.  In particular, note that for use in sla_AOP, sla_AOPPA
*     and sla_OAP the sign of the longitude must be reversed.
*
*     Users are urged to inform the author of any improvements
*     they would like to see made.  For example:
*
*         typographical corrections
*         more accurate parameters
*         better station identifiers or names
*         additional stations
*
*  P.T.Wallace   Starlink   15 March 2002
*
*  Copyright (C) 2002 Rutherford Appleton Laboratory
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

      INTEGER N
      CHARACTER C*(*),NAME*(*)
      DOUBLE PRECISION W,P,H

      INTEGER NMAX,M,NS,I
      CHARACTER*10 CC

      DOUBLE PRECISION AS2R,WEST,NORTH,EAST,SOUTH
      INTEGER ID,IAM
      REAL AS
      PARAMETER (AS2R=0.484813681109535994D-5)

*  Table of station identifiers
      PARAMETER (NMAX=84)
      CHARACTER*10 CTAB(NMAX)
      DATA CTAB  (1) /'AAT       '/
      DATA CTAB  (2) /'LPO4.2    '/
      DATA CTAB  (3) /'LPO2.5    '/
      DATA CTAB  (4) /'LPO1      '/
      DATA CTAB  (5) /'LICK120   '/
      DATA CTAB  (6) /'MMT       '/
      DATA CTAB  (7) /'DAO72     '/
      DATA CTAB  (8) /'DUPONT    '/
      DATA CTAB  (9) /'MTHOP1.5  '/
      DATA CTAB (10) /'STROMLO74 '/
      DATA CTAB (11) /'ANU2.3    '/
      DATA CTAB (12) /'GBVA140   '/
      DATA CTAB (13) /'TOLOLO4M  '/
      DATA CTAB (14) /'TOLOLO1.5M'/
      DATA CTAB (15) /'TIDBINBLA '/
      DATA CTAB (16) /'BLOEMF    '/
      DATA CTAB (17) /'BOSQALEGRE'/
      DATA CTAB (18) /'FLAGSTF61 '/
      DATA CTAB (19) /'LOWELL72  '/
      DATA CTAB (20) /'HARVARD   '/
      DATA CTAB (21) /'OKAYAMA   '/
      DATA CTAB (22) /'KPNO158   '/
      DATA CTAB (23) /'KPNO90    '/
      DATA CTAB (24) /'KPNO84    '/
      DATA CTAB (25) /'KPNO36FT  '/
      DATA CTAB (26) /'KOTTAMIA  '/
      DATA CTAB (27) /'ESO3.6    '/
      DATA CTAB (28) /'MAUNAK88  '/
      DATA CTAB (29) /'UKIRT     '/
      DATA CTAB (30) /'QUEBEC1.6 '/
      DATA CTAB (31) /'MTEKAR    '/
      DATA CTAB (32) /'MTLEMMON60'/
      DATA CTAB (33) /'MCDONLD2.7'/
      DATA CTAB (34) /'MCDONLD2.1'/
      DATA CTAB (35) /'PALOMAR200'/
      DATA CTAB (36) /'PALOMAR60 '/
      DATA CTAB (37) /'DUNLAP74  '/
      DATA CTAB (38) /'HPROV1.93 '/
      DATA CTAB (39) /'HPROV1.52 '/
      DATA CTAB (40) /'SANPM83   '/
      DATA CTAB (41) /'SAAO74    '/
      DATA CTAB (42) /'TAUTNBG   '/
      DATA CTAB (43) /'CATALINA61'/
      DATA CTAB (44) /'STEWARD90 '/
      DATA CTAB (45) /'USSR6     '/
      DATA CTAB (46) /'ARECIBO   '/
      DATA CTAB (47) /'CAMB5KM   '/
      DATA CTAB (48) /'CAMB1MILE '/
      DATA CTAB (49) /'EFFELSBERG'/
      DATA CTAB (50) /'GBT       '/
      DATA CTAB (51) /'JODRELL1  '/
      DATA CTAB (52) /'PARKES    '/
      DATA CTAB (53) /'VLA       '/
      DATA CTAB (54) /'SUGARGROVE'/
      DATA CTAB (55) /'USSR600   '/
      DATA CTAB (56) /'NOBEYAMA  '/
      DATA CTAB (57) /'JCMT      '/
      DATA CTAB (58) /'ESONTT    '/
      DATA CTAB (59) /'ST.ANDREWS'/
      DATA CTAB (60) /'APO3.5    '/
      DATA CTAB (61) /'KECK1     '/
      DATA CTAB (62) /'TAUTSCHM  '/
      DATA CTAB (63) /'PALOMAR48 '/
      DATA CTAB (64) /'UKST      '/
      DATA CTAB (65) /'KISO      '/
      DATA CTAB (66) /'ESOSCHM   '/
      DATA CTAB (67) /'ATCA      '/
      DATA CTAB (68) /'MOPRA     '/
      DATA CTAB (69) /'SUBARU    '/
      DATA CTAB (70) /'CFHT      '/
      DATA CTAB (71) /'KECK2     '/
      DATA CTAB (72) /'GEMININ   '/
      DATA CTAB (73) /'FCRAO     '/
      DATA CTAB (74) /'IRTF      '/
      DATA CTAB (75) /'CSO       '/
      DATA CTAB (76) /'VLT1      '/
      DATA CTAB (77) /'VLT2      '/
      DATA CTAB (78) /'VLT3      '/
      DATA CTAB (79) /'VLT4      '/
      DATA CTAB (80) /'GEMINIS   '/
      DATA CTAB (81) /'KOSMA3M   '/
      DATA CTAB (82) /'MAGELLAN1 '/
      DATA CTAB (83) /'MAGELLAN2 '/
      DATA CTAB (84) /'NRT       '/

*  Degrees, arcminutes, arcseconds to radians
      WEST(ID,IAM,AS)=AS2R*(DBLE(60*(60*ID+IAM))+DBLE(AS))
      NORTH(ID,IAM,AS)=WEST(ID,IAM,AS)
      EAST(ID,IAM,AS)=-WEST(ID,IAM,AS)
      SOUTH(ID,IAM,AS)=-WEST(ID,IAM,AS)




*  Station specified by number or identifier?
      IF (N.GT.0) THEN

*     Station specified by number
         M=N
         IF (M.LE.NMAX) C=CTAB(M)

      ELSE

*     Station specified by identifier:  determine corresponding number
         CC=C
         DO NS=1,NMAX
            DO I=1,10
               IF (CC(I:I).EQ.' ') GO TO 5
               IF (CC(I:I).NE.CTAB(NS)(I:I)) GO TO 1
            END DO
            GO TO 5
 1          CONTINUE
         END DO
         NS=NMAX+1
 5       CONTINUE
         IF (C(1:1).NE.' ') THEN
            M=NS
         ELSE
            M=NMAX+1
         END IF

      END IF

*
*  Return parameters of Mth station
*  --------------------------------

      GO TO (10,20,30,40,50,60,70,80,90,100,
     :       110,120,130,140,150,160,170,180,190,200,
     :       210,220,230,240,250,260,270,280,290,300,
     :       310,320,330,340,350,360,370,380,390,400,
     :       410,420,430,440,450,460,470,480,490,500,
     :       510,520,530,540,550,560,570,580,590,600,
     :       610,620,630,640,650,660,670,680,690,700,
     :       710,720,730,740,750,760,770,780,790,800,
     :       810,820,830,840) M
      GO TO 9000

*  AAT (Observer's Guide)                                            AAT
 10   CONTINUE
      NAME='Anglo-Australian 3.9m Telescope'
      W=EAST(149,03,57.91)
      P=SOUTH(31,16,37.34)
      H=1164D0
      GO TO 9999

*  WHT (Gemini, April 1987)                                       LPO4.2
 20   CONTINUE
      NAME='William Herschel 4.2m Telescope'
      W=WEST(17,52,53.9)
      P=NORTH(28,45,38.1)
      H=2332D0
      GO TO 9999

*  INT (Gemini, April 1987)                                       LPO2.5
 30   CONTINUE
      NAME='Isaac Newton 2.5m Telescope'
      W=WEST(17,52,39.5)
      P=NORTH(28,45,43.2)
      H=2336D0
      GO TO 9999

*  JKT (Gemini, April 1987)                                         LPO1
 40   CONTINUE
      NAME='Jacobus Kapteyn 1m Telescope'
      W=WEST(17,52,41.2)
      P=NORTH(28,45,39.9)
      H=2364D0
      GO TO 9999

*  Lick 120" (S.L.Allen, private communication, 2002)            LICK120
 50   CONTINUE
      NAME='Lick 120 inch'
      W=WEST(121,38,13.689)
      P=NORTH(37,20,34.931)
      H=1286D0
      GO TO 9999

*  MMT 6.5m conversion (MMT Observatory website)                     MMT
 60   CONTINUE
      NAME='MMT 6.5m, Mt Hopkins'
      W=WEST(110,53,04.4)
      P=NORTH(31,41,19.6)
      H=2608D0
      GO TO 9999

*  Victoria B.C. 1.85m (1984 Almanac)                              DAO72
 70   CONTINUE
      NAME='DAO Victoria BC 1.85 metre'
      W=WEST(123,25,01.18)
      P=NORTH(48,31,11.9)
      H=238D0
      GO TO 9999

*  Las Campanas (1983 Almanac)                                    DUPONT
 80   CONTINUE
      NAME='Du Pont 2.5m Telescope, Las Campanas'
      W=WEST(70,42,9.)
      P=SOUTH(29,00,11.)
      H=2280D0
      GO TO 9999

*  Mt Hopkins 1.5m (1983 Almanac)                               MTHOP1.5
 90   CONTINUE
      NAME='Mt Hopkins 1.5 metre'
      W=WEST(110,52,39.00)
      P=NORTH(31,40,51.4)
      H=2344D0
      GO TO 9999

*  Mt Stromlo 74" (1983 Almanac)                               STROMLO74
 100  CONTINUE
      NAME='Mount Stromlo 74 inch'
      W=EAST(149,00,27.59)
      P=SOUTH(35,19,14.3)
      H=767D0
      GO TO 9999

*  ANU 2.3m, SSO (Gary Hovey)                                     ANU2.3
 110  CONTINUE
      NAME='Siding Spring 2.3 metre'
      W=EAST(149,03,40.3)
      P=SOUTH(31,16,24.1)
      H=1149D0
      GO TO 9999

*  Greenbank 140' (1983 Almanac)                                 GBVA140
 120  CONTINUE
      NAME='Greenbank 140 foot'
      W=WEST(79,50,09.61)
      P=NORTH(38,26,15.4)
      H=881D0
      GO TO 9999

*  Cerro Tololo 4m (1982 Almanac)                               TOLOLO4M
 130  CONTINUE
      NAME='Cerro Tololo 4 metre'
      W=WEST(70,48,53.6)
      P=SOUTH(30,09,57.8)
      H=2235D0
      GO TO 9999

*  Cerro Tololo 1.5m (1982 Almanac)                           TOLOLO1.5M
 140  CONTINUE
      NAME='Cerro Tololo 1.5 metre'
      W=WEST(70,48,54.5)
      P=SOUTH(30,09,56.3)
      H=2225D0
      GO TO 9999

*  Tidbinbilla 64m (1982 Almanac)                              TIDBINBLA
 150  CONTINUE
      NAME='Tidbinbilla 64 metre'
      W=EAST(148,58,48.20)
      P=SOUTH(35,24,14.3)
      H=670D0
      GO TO 9999

*  Bloemfontein 1.52m (1981 Almanac)                              BLOEMF
 160  CONTINUE
      NAME='Bloemfontein 1.52 metre'
      W=EAST(26,24,18.)
      P=SOUTH(29,02,18.)
      H=1387D0
      GO TO 9999

*  Bosque Alegre 1.54m (1981 Almanac)                         BOSQALEGRE
 170  CONTINUE
      NAME='Bosque Alegre 1.54 metre'
      W=WEST(64,32,48.0)
      P=SOUTH(31,35,53.)
      H=1250D0
      GO TO 9999

*  USNO 61" astrographic reflector, Flagstaff (1981 Almanac)   FLAGSTF61
 180  CONTINUE
      NAME='USNO 61 inch astrograph, Flagstaff'
      W=WEST(111,44,23.6)
      P=NORTH(35,11,02.5)
      H=2316D0
      GO TO 9999

*  Lowell 72" (1981 Almanac)                                    LOWELL72
 190  CONTINUE
      NAME='Perkins 72 inch, Lowell'
      W=WEST(111,32,09.3)
      P=NORTH(35,05,48.6)
      H=2198D0
      GO TO 9999

*  Harvard 1.55m (1981 Almanac)                                  HARVARD
 200  CONTINUE
      NAME='Harvard College Observatory 1.55m'
      W=WEST(71,33,29.32)
      P=NORTH(42,30,19.0)
      H=185D0
      GO TO 9999

*  Okayama 1.88m (1981 Almanac)                                  OKAYAMA
 210  CONTINUE
      NAME='Okayama 1.88 metre'
      W=EAST(133,35,47.29)
      P=NORTH(34,34,26.1)
      H=372D0
      GO TO 9999

*  Kitt Peak Mayall 4m (1981 Almanac)                            KPNO158
 220  CONTINUE
      NAME='Kitt Peak 158 inch'
      W=WEST(111,35,57.61)
      P=NORTH(31,57,50.3)
      H=2120D0
      GO TO 9999

*  Kitt Peak 90 inch (1981 Almanac)                               KPNO90
 230  CONTINUE
      NAME='Kitt Peak 90 inch'
      W=WEST(111,35,58.24)
      P=NORTH(31,57,46.9)
      H=2071D0
      GO TO 9999

*  Kitt Peak 84 inch (1981 Almanac)                               KPNO84
 240  CONTINUE
      NAME='Kitt Peak 84 inch'
      W=WEST(111,35,51.56)
      P=NORTH(31,57,29.2)
      H=2096D0
      GO TO 9999

*  Kitt Peak 36 foot (1981 Almanac)                             KPNO36FT
 250  CONTINUE
      NAME='Kitt Peak 36 foot'
      W=WEST(111,36,51.12)
      P=NORTH(31,57,12.1)
      H=1939D0
      GO TO 9999

*  Kottamia 74" (1981 Almanac)                                  KOTTAMIA
 260  CONTINUE
      NAME='Kottamia 74 inch'
      W=EAST(31,49,30.)
      P=NORTH(29,55,54.)
      H=476D0
      GO TO 9999

*  La Silla 3.6m (1981 Almanac)                                   ESO3.6
 270  CONTINUE
      NAME='ESO 3.6 metre'
      W=WEST(70,43,36.)
      P=SOUTH(29,15,36.)
      H=2428D0
      GO TO 9999

*  Mauna Kea 88 inch                                            MAUNAK88
*  (IfA website, Richard Wainscoat)
 280  CONTINUE
      NAME='Mauna Kea 88 inch'
      W=WEST(155,28,09.96)
      P=NORTH(19,49,22.77)
      H=4213.6D0
      GO TO 9999

*  UKIRT (IfA website, Richard Wainscoat)                          UKIRT
 290  CONTINUE
      NAME='UK Infra Red Telescope'
      W=WEST(155,28,13.18)
      P=NORTH(19,49,20.75)
      H=4198.5D0
      GO TO 9999

*  Quebec 1.6m (1981 Almanac)                                  QUEBEC1.6
 300  CONTINUE
      NAME='Quebec 1.6 metre'
      W=WEST(71,09,09.7)
      P=NORTH(45,27,20.6)
      H=1114D0
      GO TO 9999

*  Mt Ekar 1.82m (1981 Almanac)                                   MTEKAR
 310  CONTINUE
      NAME='Mt Ekar 1.82 metre'
      W=EAST(11,34,15.)
      P=NORTH(45,50,48.)
      H=1365D0
      GO TO 9999

*  Mt Lemmon 60" (1981 Almanac)                               MTLEMMON60
 320  CONTINUE
      NAME='Mt Lemmon 60 inch'
      W=WEST(110,42,16.9)
      P=NORTH(32,26,33.9)
      H=2790D0
      GO TO 9999

*  Mt Locke 2.7m (1981 Almanac)                               MCDONLD2.7
 330  CONTINUE
      NAME='McDonald 2.7 metre'
      W=WEST(104,01,17.60)
      P=NORTH(30,40,17.7)
      H=2075D0
      GO TO 9999

*  Mt Locke 2.1m (1981 Almanac)                               MCDONLD2.1
 340  CONTINUE
      NAME='McDonald 2.1 metre'
      W=WEST(104,01,20.10)
      P=NORTH(30,40,17.7)
      H=2075D0
      GO TO 9999

*  Palomar 200" (1981 Almanac)                                PALOMAR200
 350  CONTINUE
      NAME='Palomar 200 inch'
      W=WEST(116,51,50.)
      P=NORTH(33,21,22.)
      H=1706D0
      GO TO 9999

*  Palomar 60" (1981 Almanac)                                  PALOMAR60
 360  CONTINUE
      NAME='Palomar 60 inch'
      W=WEST(116,51,31.)
      P=NORTH(33,20,56.)
      H=1706D0
      GO TO 9999

*  David Dunlap 74" (1981 Almanac)                              DUNLAP74
 370  CONTINUE
      NAME='David Dunlap 74 inch'
      W=WEST(79,25,20.)
      P=NORTH(43,51,46.)
      H=244D0
      GO TO 9999

*  Haute Provence 1.93m (1981 Almanac)                         HPROV1.93
 380  CONTINUE
      NAME='Haute Provence 1.93 metre'
      W=EAST(5,42,46.75)
      P=NORTH(43,55,53.3)
      H=665D0
      GO TO 9999

*  Haute Provence 1.52m (1981 Almanac)                         HPROV1.52
 390  CONTINUE
      NAME='Haute Provence 1.52 metre'
      W=EAST(5,42,43.82)
      P=NORTH(43,56,00.2)
      H=667D0
      GO TO 9999

*  San Pedro Martir 83" (1981 Almanac)                           SANPM83
 400  CONTINUE
      NAME='San Pedro Martir 83 inch'
      W=WEST(115,27,47.)
      P=NORTH(31,02,38.)
      H=2830D0
      GO TO 9999

*  Sutherland 74" (1981 Almanac)                                  SAAO74
 410  CONTINUE
      NAME='Sutherland 74 inch'
      W=EAST(20,48,44.3)
      P=SOUTH(32,22,43.4)
      H=1771D0
      GO TO 9999

*  Tautenburg 2m (1981 Almanac)                                  TAUTNBG
 420  CONTINUE
      NAME='Tautenburg 2 metre'
      W=EAST(11,42,45.)
      P=NORTH(50,58,51.)
      H=331D0
      GO TO 9999

*  Catalina 61" (1981 Almanac)                                CATALINA61
 430  CONTINUE
      NAME='Catalina 61 inch'
      W=WEST(110,43,55.1)
      P=NORTH(32,25,00.7)
      H=2510D0
      GO TO 9999

*  Steward 90" (1981 Almanac)                                  STEWARD90
 440  CONTINUE
      NAME='Steward 90 inch'
      W=WEST(111,35,58.24)
      P=NORTH(31,57,46.9)
      H=2071D0
      GO TO 9999

*  Russian 6m (1981 Almanac)                                       USSR6
 450  CONTINUE
      NAME='USSR 6 metre'
      W=EAST(41,26,30.0)
      P=NORTH(43,39,12.)
      H=2100D0
      GO TO 9999

*  Arecibo 1000' (1981 Almanac)                                  ARECIBO
 460  CONTINUE
      NAME='Arecibo 1000 foot'
      W=WEST(66,45,11.1)
      P=NORTH(18,20,36.6)
      H=496D0
      GO TO 9999

*  Cambridge 5km (1981 Almanac)                                  CAMB5KM
 470  CONTINUE
      NAME='Cambridge 5km'
      W=EAST(0,02,37.23)
      P=NORTH(52,10,12.2)
      H=17D0
      GO TO 9999

*  Cambridge 1 mile (1981 Almanac)                             CAMB1MILE
 480  CONTINUE
      NAME='Cambridge 1 mile'
      W=EAST(0,02,21.64)
      P=NORTH(52,09,47.3)
      H=17D0
      GO TO 9999

*  Bonn 100m (1981 Almanac)                                   EFFELSBERG
 490  CONTINUE
      NAME='Effelsberg 100 metre'
      W=EAST(6,53,01.5)
      P=NORTH(50,31,28.6)
      H=366D0
      GO TO 9999

*  Green Bank Telescop 100m
 500  CONTINUE
      NAME='Green Bank Telescope'
      W=WEST(79,50,23.406)
      P=NORTH(38,25,59.236)
      H=880D0
      GO TO 9999

*  Jodrell Bank Mk 1 (1981 Almanac)                             JODRELL1
 510  CONTINUE
      NAME='Jodrell Bank 250 foot'
      W=WEST(2,18,25.)
      P=NORTH(53,14,10.5)
      H=78D0
      GO TO 9999

*  Australia Telescope Parkes Observatory                         PARKES
*  (Peter te Lintel Hekkert)
 520  CONTINUE
      NAME='Parkes 64 metre'
      W=EAST(148,15,44.3591)
      P=SOUTH(32,59,59.8657)
      H=391.79D0
      GO TO 9999

*  VLA (1981 Almanac)                                                VLA
 530  CONTINUE
      NAME='Very Large Array'
      W=WEST(107,37,03.82)
      P=NORTH(34,04,43.5)
      H=2124D0
      GO TO 9999

*  Sugar Grove 150' (1981 Almanac)                            SUGARGROVE
 540  CONTINUE
      NAME='Sugar Grove 150 foot'
      W=WEST(79,16,23.)
      P=NORTH(38,31,14.)
      H=705D0
      GO TO 9999

*  Russian 600' (1981 Almanac)                                   USSR600
 550  CONTINUE
      NAME='USSR 600 foot'
      W=EAST(41,35,25.5)
      P=NORTH(43,49,32.)
      H=973D0
      GO TO 9999

*  Nobeyama 45 metre mm dish (based on 1981 Almanac entry)      NOBEYAMA
 560  CONTINUE
      NAME='Nobeyama 45 metre'
      W=EAST(138,29,12.)
      P=NORTH(35,56,19.)
      H=1350D0
      GO TO 9999

*  James Clerk Maxwell 15 metre mm telescope, Mauna Kea             JCMT
*  (IfA website, Richard Wainscoat, height from I.Coulson)
 570  CONTINUE
      NAME='JCMT 15 metre'
      W=WEST(155,28,37.20)
      P=NORTH(19,49,22.11)
      H=4111D0
      GO TO 9999

*  ESO 3.5 metre NTT, La Silla (K.Wirenstrand)                    ESONTT
 580  CONTINUE
      NAME='ESO 3.5 metre NTT'
      W=WEST(70,43,07.)
      P=SOUTH(29,15,30.)
      H=2377D0
      GO TO 9999

*  St Andrews University Observatory (1982 Almanac)           ST.ANDREWS
 590  CONTINUE
      NAME='St Andrews'
      W=WEST(2,48,52.5)
      P=NORTH(56,20,12.)
      H=30D0
      GO TO 9999

*  Apache Point 3.5 metre (R.Owen)                                APO3.5
 600  CONTINUE
      NAME='Apache Point 3.5m'
      W=WEST(105,49,11.56)
      P=NORTH(32,46,48.96)
      H=2809D0
      GO TO 9999

*  W.M.Keck Observatory, Telescope 1                               KECK1
*  (William Lupton)
 610  CONTINUE
      NAME='Keck 10m Telescope #1'
      W=WEST(155,28,28.99)
      P=NORTH(19,49,33.41)
      H=4160D0
      GO TO 9999

*  Tautenberg Schmidt (1983 Almanac)                            TAUTSCHM
 620  CONTINUE
      NAME='Tautenberg 1.34 metre Schmidt'
      W=EAST(11,42,45.0)
      P=NORTH(50,58,51.0)
      H=331D0
      GO TO 9999

*  Palomar Schmidt (1981 Almanac)                              PALOMAR48
 630  CONTINUE
      NAME='Palomar 48-inch Schmidt'
      W=WEST(116,51,32.0)
      P=NORTH(33,21,26.0)
      H=1706D0
      GO TO 9999

*  UK Schmidt, Siding Spring (1983 Almanac)                         UKST
 640  CONTINUE
      NAME='UK 1.2 metre Schmidt, Siding Spring'
      W=EAST(149,04,12.8)
      P=SOUTH(31,16,27.8)
      H=1145D0
      GO TO 9999

*  Kiso Schmidt, Japan (1981 Almanac)                               KISO
 650  CONTINUE
      NAME='Kiso 1.05 metre Schmidt, Japan'
      W=EAST(137,37,42.2)
      P=NORTH(35,47,38.7)
      H=1130D0
      GO TO 9999

*  ESO Schmidt, La Silla (1981 Almanac)                          ESOSCHM
 660  CONTINUE
      NAME='ESO 1 metre Schmidt, La Silla'
      W=WEST(70,43,46.5)
      P=SOUTH(29,15,25.8)
      H=2347D0
      GO TO 9999

*  Australia Telescope Compact Array                                ATCA
*  (WGS84 coordinates of Station 35, Mark Calabretta)
 670  CONTINUE
      NAME='Australia Telescope Compact Array'
      W=EAST(149,33,00.500)
      P=SOUTH(30,18,46.385)
      H=236.9D0
      GO TO 9999

*  Australia Telescope Mopra Observatory                           MOPRA
*  (Peter te Lintel Hekkert)
 680  CONTINUE
      NAME='ATNF Mopra Observatory'
      W=EAST(149,05,58.732)
      P=SOUTH(31,16,04.451)
      H=850D0
      GO TO 9999

*  Subaru telescope, Mauna Kea                                     SUBARU
*  (IfA website, Richard Wainscoat)
 690  CONTINUE
      NAME='Subaru 8m telescope'
      W=WEST(155,28,33.67)
      P=NORTH(19,49,31.81)
      H=4163D0
      GO TO 9999

*  Canada-France-Hawaii Telescope, Mauna Kea                         CFHT
*  (IfA website, Richard Wainscoat)
 700  CONTINUE
      NAME='Canada-France-Hawaii 3.6m Telescope'
      W=WEST(155,28,07.95)
      P=NORTH(19,49,30.91)
      H=4204.1D0
      GO TO 9999

*  W.M.Keck Observatory, Telescope 2                                KECK2
*  (William Lupton)
 710  CONTINUE
      NAME='Keck 10m Telescope #2'
      W=WEST(155,28,27.24)
      P=NORTH(19,49,35.62)
      H=4159.6D0
      GO TO 9999

*  Gemini North, Mauna Kea                                        GEMININ
*  (IfA website, Richard Wainscoat)
 720  CONTINUE
      NAME='Gemini North 8-m telescope'
      W=WEST(155,28,08.57)
      P=NORTH(19,49,25.69)
      H=4213.4D0
      GO TO 9999

*  Five College Radio Astronomy Observatory                        FCRAO
*  (Tim Jenness)
 730  CONTINUE
      NAME='Five College Radio Astronomy Obs'
      W=WEST(72,20,42.0)
      P=NORTH(42,23,30.0)
      H=314D0
      GO TO 9999

*  NASA Infra Red Telescope Facility                                IRTF
*  (IfA website, Richard Wainscoat)
 740  CONTINUE
      NAME='NASA IR Telescope Facility, Mauna Kea'
      W=WEST(155,28,19.20)
      P=NORTH(19,49,34.39)
      H=4168.1D0
      GO TO 9999

*  Caltech Submillimeter Observatory                                 CSO
*  (IfA website, Richard Wainscoat; height estimated)
 750  CONTINUE
      NAME='Caltech Sub-mm Observatory, Mauna Kea'
      W=WEST(155,28,31.79)
      P=NORTH(19,49,20.78)
      H=4080D0
      GO TO 9999

* ESO VLT, UT1                                                       VLT1
* (ESO website, VLT Whitebook Chapter 2)
 760  CONTINUE
      NAME='ESO VLT, Paranal, Chile: UT1'
      W=WEST(70,24,11.642)
      P=SOUTH(24,37,33.117)
      H=2635.43
      GO TO 9999

* ESO VLT, UT2                                                       VLT2
* (ESO website, VLT Whitebook Chapter 2)
 770  CONTINUE
      NAME='ESO VLT, Paranal, Chile: UT2'
      W=WEST(70,24,10.855)
      P=SOUTH(24,37,31.465)
      H=2635.43
      GO TO 9999

* ESO VLT, UT3                                                       VLT3
* (ESO website, VLT Whitebook Chapter 2)
 780  CONTINUE
      NAME='ESO VLT, Paranal, Chile: UT3'
      W=WEST(70,24,09.896)
      P=SOUTH(24,37,30.300)
      H=2635.43
      GO TO 9999

* ESO VLT, UT4                                                       VLT4
* (ESO website, VLT Whitebook Chapter 2)
 790  CONTINUE
      NAME='ESO VLT, Paranal, Chile: UT4'
      W=WEST(70,24,08.000)
      P=SOUTH(24,37,31.000)
      H=2635.43
      GO TO 9999

*  Gemini South, Cerro Pachon                                     GEMINIS
*  (GPS readings by Patrick Wallace)
 800  CONTINUE
      NAME='Gemini South 8-m telescope'
      W=WEST(70,44,11.5)
      P=SOUTH(30,14,26.7)
      H=2738D0
      GO TO 9999

*  Cologne Observatory for Submillimeter Astronomy (KOSMA)        KOSMA3M
*  (Holger Jakob)
 810  CONTINUE
      NAME='KOSMA 3m telescope, Gornergrat'
      W=EAST(7,47,3.48)
      P=NORTH(45,58,59.772)
      H=3141D0
      GO TO 9999

*  Magellan 1, 6.5m telescope at Las Campanas, Chile            MAGELLAN1
*  (Skip Schaller)
 820  CONTINUE
      NAME='Magellan 1, 6.5m, Las Campanas'
      W=WEST(70,41,31.9)
      P=SOUTH(29,00,51.7)
      H=2408D0
      GO TO 9999

*  Magellan 2, 6.5m telescope at Las Campanas, Chile            MAGELLAN2
*  (Skip Schaller)
 830  CONTINUE
      NAME='Magellan 2, 6.5m, Las Campanas'
      W=WEST(70,41,33.5)
      P=SOUTH(29,00,50.3)
      H=2408D0
      GO TO 9999

*  Nancay Radio Telescope, France, 94m equivalent               NRT
 840  CONTINUE
      NAME='Nancay Radio Telescope, France, 94m equivalent'
      W=EAST(2,11,50.92)
      P=NORTH(47,22,24.97)
      H=191.0
      GO TO 9999

*  Unrecognized station
 9000 CONTINUE
      NAME='?'

*  Exit
 9999 CONTINUE

      END
