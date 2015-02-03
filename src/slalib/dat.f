      DOUBLE PRECISION FUNCTION sla_DAT (UTC)
*+
*     - - - -
*      D A T
*     - - - -
*
*  Increment to be applied to Coordinated Universal Time UTC to give
*  International Atomic Time TAI (double precision)
*
*  Given:
*     UTC      d      UTC date as a modified JD (JD-2400000.5)
*
*  Result:  TAI-UTC in seconds
*
*  Notes:
*
*  1  The UTC is specified to be a date rather than a time to indicate
*     that care needs to be taken not to specify an instant which lies
*     within a leap second.  Though in most cases UTC can include the
*     fractional part, correct behaviour on the day of a leap second
*     can only be guaranteed up to the end of the second 23:59:59.
*
*  2  For epochs from 1961 January 1 onwards, the expressions from the
*     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
*
*  3  The 5ms time step at 1961 January 1 is taken from 2.58.1 (p87) of
*     the 1992 Explanatory Supplement.
*
*  4  UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
*     to call the routine with an earlier epoch.  However, if this
*     is attempted, the TAI-UTC expression for the year 1960 is used.
*
*
*     :-----------------------------------------:
*     :                                         :
*     :                IMPORTANT                :
*     :                                         :
*     :  This routine must be updated on each   :
*     :     occasion that a leap second is      :
*     :                announced                :
*     :                                         :
*     :  Latest leap second:  2015 July 1       :
*     :                                         :
*     :-----------------------------------------:
*
*  Last revision:   31 January 2015
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

      DOUBLE PRECISION UTC

      DOUBLE PRECISION DT



      IF (.FALSE.) THEN

* - - - - - - - - - - - - - - - - - - - - - - *
*  Add new code here on each occasion that a  *
*  leap second is announced, and update the   *
*  preamble comments appropriately.           *
* - - - - - - - - - - - - - - - - - - - - - - *

*     2015 July 1
      ELSE IF (UTC.GE.57204D0) THEN
         DT=36D0

*     2012 July 1
      ELSE IF (UTC.GE.56109D0) THEN
         DT=35D0

*     2009 January 1
      ELSE IF (UTC.GE.54832D0) THEN
         DT=34D0

*     2006 January 1
      ELSE IF (UTC.GE.53736D0) THEN
         DT=33D0

*     1999 January 1
      ELSE IF (UTC.GE.51179D0) THEN
         DT=32D0

*     1997 July 1
      ELSE IF (UTC.GE.50630D0) THEN
         DT=31D0

*     1996 January 1
      ELSE IF (UTC.GE.50083D0) THEN
         DT=30D0

*     1994 July 1
      ELSE IF (UTC.GE.49534D0) THEN
         DT=29D0

*     1993 July 1
      ELSE IF (UTC.GE.49169D0) THEN
         DT=28D0

*     1992 July 1
      ELSE IF (UTC.GE.48804D0) THEN
         DT=27D0

*     1991 January 1
      ELSE IF (UTC.GE.48257D0) THEN
         DT=26D0

*     1990 January 1
      ELSE IF (UTC.GE.47892D0) THEN
         DT=25D0

*     1988 January 1
      ELSE IF (UTC.GE.47161D0) THEN
         DT=24D0

*     1985 July 1
      ELSE IF (UTC.GE.46247D0) THEN
         DT=23D0

*     1983 July 1
      ELSE IF (UTC.GE.45516D0) THEN
         DT=22D0

*     1982 July 1
      ELSE IF (UTC.GE.45151D0) THEN
         DT=21D0

*     1981 July 1
      ELSE IF (UTC.GE.44786D0) THEN
         DT=20D0

*     1980 January 1
      ELSE IF (UTC.GE.44239D0) THEN
         DT=19D0

*     1979 January 1
      ELSE IF (UTC.GE.43874D0) THEN
         DT=18D0

*     1978 January 1
      ELSE IF (UTC.GE.43509D0) THEN
         DT=17D0

*     1977 January 1
      ELSE IF (UTC.GE.43144D0) THEN
         DT=16D0

*     1976 January 1
      ELSE IF (UTC.GE.42778D0) THEN
         DT=15D0

*     1975 January 1
      ELSE IF (UTC.GE.42413D0) THEN
         DT=14D0

*     1974 January 1
      ELSE IF (UTC.GE.42048D0) THEN
         DT=13D0

*     1973 January 1
      ELSE IF (UTC.GE.41683D0) THEN
         DT=12D0

*     1972 July 1
      ELSE IF (UTC.GE.41499D0) THEN
         DT=11D0

*     1972 January 1
      ELSE IF (UTC.GE.41317D0) THEN
         DT=10D0

*     1968 February 1
      ELSE IF (UTC.GE.39887D0) THEN
         DT=4.2131700D0+(UTC-39126D0)*0.002592D0

*     1966 January 1
      ELSE IF (UTC.GE.39126D0) THEN
         DT=4.3131700D0+(UTC-39126D0)*0.002592D0

*     1965 September 1
      ELSE IF (UTC.GE.39004D0) THEN
         DT=3.8401300D0+(UTC-38761D0)*0.001296D0

*     1965 July 1
      ELSE IF (UTC.GE.38942D0) THEN
         DT=3.7401300D0+(UTC-38761D0)*0.001296D0

*     1965 March 1
      ELSE IF (UTC.GE.38820D0) THEN
         DT=3.6401300D0+(UTC-38761D0)*0.001296D0

*     1965 January 1
      ELSE IF (UTC.GE.38761D0) THEN
         DT=3.5401300D0+(UTC-38761D0)*0.001296D0

*     1964 September 1
      ELSE IF (UTC.GE.38639D0) THEN
         DT=3.4401300D0+(UTC-38761D0)*0.001296D0

*     1964 April 1
      ELSE IF (UTC.GE.38486D0) THEN
         DT=3.3401300D0+(UTC-38761D0)*0.001296D0

*     1964 January 1
      ELSE IF (UTC.GE.38395D0) THEN
         DT=3.2401300D0+(UTC-38761D0)*0.001296D0

*     1963 November 1
      ELSE IF (UTC.GE.38334D0) THEN
         DT=1.9458580D0+(UTC-37665D0)*0.0011232D0

*     1962 January 1
      ELSE IF (UTC.GE.37665D0) THEN
         DT=1.8458580D0+(UTC-37665D0)*0.0011232D0

*     1961 August 1
      ELSE IF (UTC.GE.37512D0) THEN
         DT=1.3728180D0+(UTC-37300D0)*0.001296D0

*     1961 January 1
      ELSE IF (UTC.GE.37300D0) THEN
         DT=1.4228180D0+(UTC-37300D0)*0.001296D0

*     Before that
      ELSE
         DT=1.4178180D0+(UTC-37300D0)*0.001296D0

      END IF

      sla_DAT=DT

      END
