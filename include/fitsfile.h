/*** File fitsfile.h  FITS and IRAF file access subroutines
 *** March 5, 2003
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 1996-2003
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Correspondence concerning WCSTools should be addressed as follows:
           Internet email: dmink@cfa.harvard.edu
           Postal address: Doug Mink
                           Smithsonian Astrophysical Observatory
                           60 Garden St.
                           Cambridge, MA 02138 USA
 */

#ifndef fitsfile_h_
#define fitsfile_h_
#include "fitshead.h"

/* Declarations for subroutines in fitsfile.c, imhfile.c, imio.c,
 * fileutil.c, and dateutil.c */

/* FITS table keyword structure */
struct Keyword {
    char kname[10];	/* Keyword for table entry */
    int lname;		/* Length of keyword name */
    int kn;		/* Index of entry on line */
    int kf;		/* Index in line of first character of entry */
    int kl;		/* Length of entry value */
};

#define FITSBLOCK 2880

/* FITS file access subroutines in fitsfile.c */
extern int fitsropen();
extern char *fitsrhead();
extern char *fitsrimage();
extern char *fitsrsect();
extern int fitswhead();
extern int fitswext();
extern int fitswhdu();
extern int fitswimage();
extern int fitscimage();
extern int isfits();
extern void fitserr();          /* Print FITS error message to stderr */

/* FITS table file access subroutines in fitsfile.c */
extern int fitsrtopen();
extern int fitsrthead();
extern void fitsrtlset();
extern int fitsrtline();
extern short ftgeti2();
extern int ftgeti4();
extern float ftgetr4();
extern double ftgetr8();
extern int ftgetc();

/* IRAF file access subroutines in imhfile.c */
extern char *irafrhead();
extern char *irafrimage();
extern int irafwhead();
extern int irafwimage();
extern int isiraf();
extern char *iraf2fits();
extern char *fits2iraf();

/* Image pixel access subroutines in imio.c */
extern double getpix();	/* Read one pixel from any data type 2-D array (0,0)*/
extern double getpix1(); /* Read one pixel from any data type 2-D array (1,1)*/
extern void putpix();	/* Write one pixel to any data type 2-D array (0,0)*/
extern void putpix1();	/* Write one pixel to any data type 2-D array (1,1) */
extern void addpix();	/* Add to one pixel in any data type 2-D array (0,0)*/
extern void addpix1();	/* Add to one pixel in any data type 2-D array (1,1)*/
extern void movepix();	/* Move one pixel value between two 2-D arrays (0,0) */
extern void movepix1();	/* Move one pixel value between two 2-D arrays (1,1) */
extern void getvec();	/* Read vector from 2-D array */
extern void putvec();	/* Write vector into 2-D array */
extern void imswap();	/* Swap alternating bytes in a vector */
extern void imswap2();	/* Swap bytes in a vector of 2-byte (short) integers */
extern void imswap4();	/* Reverse bytes in a vector of 4-byte numbers */
extern void imswap8();	/* Reverse bytes in a vector of 8-byte numbers */
extern int imswapped();	/* Return 1 if machine byte order is not FITS order */

/* File utilities from fileutil.c */
extern int getfilelines();
extern char *getfilebuff();
extern int getfilesize();
extern int isimlist();
extern int isimlistd();
extern int isfilelist();
extern int isfile();
extern int first_token();

/* Subroutines for translating dates and times in dateutil.c */
void doy2dt();	/* year and day of year to yyyy.mmdd hh.mmss */
double doy2ep(); /* year and day of year to fractional year (epoch) */
double doy2epb(); /* year and day of year to Besselian epoch */
double doy2epj(); /* year and day of year to Julian epoch */
char *doy2fd();	/* year and day of year to FITS date */
double doy2jd(); /* year and day of year to Julian date */
double doy2mjd(); /* year and day of year to modified Julian date */
double doy2ts(); /* year and day of year to seconds since 1950.0 */ 
int doy2tsi();	/* year and day of year to IRAF seconds since 1980-01-01 */
time_t doy2tsu();	/* year and day of year to Unix seconds since 1970-01-01 */
void dt2doy();	/* yyyy.mmdd hh.mmss to year and day of year */
double dt2ep();	/* yyyy.ddmm and hh.mmsss to fractional year (epoch) */
double dt2epb(); /* yyyy.ddmm and hh.mmsss to Besselian epoch */
double dt2epj(); /* yyyy.ddmm and hh.mmsss to Julian epoch */
char *dt2fd();	/* yyyy.ddmm and hh.mmsss to FITS date string */
void dt2i();	/* yyyy.ddmm and hh.mmsss to year, month, day, hrs, min, sec */
double dt2jd();	/* yyyy.ddmm and hh.mmsss to Julian date */
double dt2mjd(); /* yyyy.ddmm and hh.mmsss to modified Julian date */
double dt2ts();	/* yyyy.ddmm and hh.mmsss to seconds since 1950.0 */ 
int dt2tsi();	/* yyyy.ddmm and hh.mmsss to IRAF seconds since 1980-01-01 */
time_t dt2tsu();	/* yyyy.ddmm and hh.mmsss to Unix seconds since 1970-01-01 */
void ep2dt();	/* Fractional year to yyyy.mmdd hh.mmssss */
void epb2dt();	/* Besselian epoch to yyyy.mmdd hh.mmssss */
void epj2dt();	/* Julian epoch to yyyy.mmdd hh.mmssss */
char *ep2fd();	/* Fractional year to FITS date string yyyy-mm-ddThh:mm:ss.ss */
char *epb2fd();	/* Besselian epoch to FITS date string yyyy-mm-ddThh:mm:ss.ss */
char *epj2fd();	/* Julian epoch to FITS date string yyyy-mm-ddThh:mm:ss.ss */
void ep2i();	/* Fractional year to year, month, day, hours, min., sec. */
void epb2i();	/* Besselian epoch to year, month, day, hours, min., sec. */
void epj2i();	/* Julian epoch to year, month, day, hours, min., sec. */
double ep2jd();	/* Fractional year to Julian Date */
double epb2jd(); /* Besselian epoch to Julian Date */
double epj2jd(); /* Julian epoch to Julian Date */
double ep2mjd(); /* Fractional year to modified Julian Date */
double epb2mjd(); /* Besselian epoch to modified Julian Date */
double epj2mjd(); /* Julian epoch to modified Julian Date */
double ep2epb(); /* Fractional year to Besselian epoch */
double ep2epj(); /* Fractional year to Julian epoch */
double epb2epj(); /* Besselian epoch to Julian epoch */
double epj2epb(); /* Julian epoch to Besselian epoch */
double ep2ts();	/* Fractional year to seconds since 1950.0 */
double epb2ts(); /* Besselian epoch to seconds since 1950.0 */
double epj2ts(); /* Julian epoch to seconds since 1950.0 */
void fd2dt();	/* FITS standard date string to Julian date */
void fd2doy();	/* FITS standard date string to year, day of year */
double fd2ep();	/* FITS standard date string to fractional year (epoch) */
double fd2epb(); /* FITS standard date string to Besselian epoch */
double fd2epj(); /* FITS standard date string to Julian epoch */
char *fd2fd();	/* Any FITS standard date string to ISO FITS date string */
char *fd2of();	/* Any FITS standard date string to old FITS date and time */
char *fd2ofd();	/* Any FITS standard date string to old FITS date string */
char *fd2oft(); /* Any FITS standard date string to old FITS time string */
void fd2i();	/* FITS standard date string to year, mon, day, hrs, min, sec */
double fd2jd();	/* FITS standard date string to Julian date */
double fd2mjd(); /* FITS standard date string to modified Julian date */
double fd2ts();	/* FITS standard date to seconds since 1950-01-01 */
int fd2tsi();	/* FITS standard date to IRAF seconds since 1980-01-01 */
time_t fd2tsu();	/* FITS standard date to Unix seconds since 1970-01-01 */
void jd2doy();	/* Julian date to year and day of year */
void jd2dt();	/* Julian date to yyyy.mmdd hh.mmssss */
double jd2ep();	/* Julian date to fractional year */
double jd2epb(); /* Julian date to Besselian epoch */
double jd2epj(); /* Julian date to Julian epoch */
char *jd2fd();	/* Julian date to FITS date string yyyy-mm-ddThh:mm:ss.ss */
void jd2i();	/* Julian date to year, month, day, hours, min., sec. */
double jd2mjd(); /* Julian date to modified Julian date */
double jd2ts();	/* Julian date to seconds since 1950.0 */
void lt2dt();	/* Current local time to date (yyyy.mmdd), time (hh.mmsss) */
char *lt2fd();	/* Current local time to FITS ISO date string */
int lt2tsi();	/* Current local time to IRAF seconds since 1980-01-01T00:00 */
time_t lt2tsu();	/* Current local time to Unix seconds since 1970-01-01T00:00 */
double lt2ts(); /* Current local time to IRAF seconds since 1950-01-01T00:00 */
void mjd2dt();	/* Modified Julian date to yyyy.mmdd hh.mmssss */
double mjd2ep(); /* Modified Julian date to fractional year */
double mjd2epb(); /* Modified Julian date to Besselian epoch */
double mjd2epj(); /* Modified Julian date to Julian epoch */
char *mjd2fd();	/* Modified Julian date to FITS date yyyy-mm-ddThh:mm:ss.ss */
void mjd2i();	/* Modified Julian date to year, month, day, hours, min, sec */
double mjd2jd(); /* Modified Julian date to Julian date */
double mjd2ts(); /* Modified Julian date to seconds since 1950.0 */
void ts2dt();	/* Seconds since 1950.0 to yyyy.mmdd hh.mmssss */
double ts2ep();	/* Seconds since 1950.0 to fractional year */
double ts2epb(); /* Seconds since 1950.0 to Besselian epoch */
double ts2epj(); /* Seconds since 1950.0 to Julian epoch */
char *ts2fd();	/* Seconds since 1950.0 to FITS date, yyyy-mm-ddT00:00:00.000 */
void ts2i();	/* Seconds since 1950.0 to year, month, day, hours, min, sec */
double ts2jd();	/* Seconds since 1950.0 to Julian date */
double ts2mjd(); /* Seconds since 1950.0 to modified Julian date */
char *tsi2fd();	/* Seconds since 1980-01-01 to FITS standard date string */
double tsi2ts(); /* Seconds since 1980-01-01 to seconds since 1950-01-01 */
double tsi2ts(); /* Seconds since 1980-01-01 to seconds since 1950-01-01 */
void tsu2dt();	/* Seconds since 1970-01-01 to date yyyy.ddmm, time hh.mmsss */
char *tsu2fd();	/* Seconds since 1970-01-01 to FITS standard date string */
double tsu2ts(); /* Seconds since 1970-01-01 to seconds since 1950-01-01 */
int tsu2tsi();	/* Seconds since 1970-01-01 to local seconds since 1980-01-01 */
int isdate();	/* Return 1 if string is FITS old or ISO date */
void ut2dt();	/* Current Universal Time to date (yyyy.mmdd), time (hh.mmsss) */
void ut2doy(); /* Current Universal Time to year, day of year */
double ut2ep(); /* Current Universal Time to fractional year */
double ut2epb(); /* Current Universal Time to Besselian Epoch */
double ut2epj(); /* Current Universal Time to Julian Epoch */
char *ut2fd();	/* Current Universal Time to FITS ISO date string */
double ut2jd();	/* Current Universal Time to Julian Date */
double ut2mjd(); /* Current Universal Time to Modified Julian Date */
int ut2tsi();	/* Current UT to IRAF seconds since 1980-01-01T00:00 */
time_t ut2tsu();	/* Current UT to Unix seconds since 1970-01-01T00:00 */
double ut2ts(); /* Current UT to IRAF seconds since 1950-01-01T00:00 */
int sts2c();	/* Replaces spaces in a string with a specified character */
int stc2s();	/* Replaces a specified character in a string with spaces */
char *et2fd();	/* ET (or TDT or TT) in FITS format to UT in FITS format */
char *fd2et();	/* UT in FITS format to ET (or TDT or TT) in FITS format */
double jd2jed(); /* Convert from Julian Date to Julian Ephemeris Date */
double jed2jd(); /* Convert from Julian Ephemeris Date to Julian Date */
double ets2ts(); /* ET in seconds since 1950-01-01 to UT in same format */
double ts2ets(); /* UT in seconds since 1950-01-01 to ET in same format */
void dt2et();	/* yyyy.ddmm and hh.mmsss to Ephemeris Time */ 
void edt2dt(); /* yyyy.ddmm and hh.mmsss Ephemeris Time to UT */ 
double utdt();	/* Compute difference between UT and dynamical time (ET-UT) */
char *fd2gst();	/* Convert from FITS UT date to Greenwich Sidereal Time */
void dt2gst();	/* Convert from UT as yyyy.mmdd hh.mmssss to Greenwich Sidereal Time */
double ts2gst(); /* Calculate Greenwich Sidereal Time given Universal Time */
char *fd2mst();	/* Convert from FITS UT date to Mean Sidereal Time */
void dt2mst();	/* Convert from UT as yyyy.mmdd hh.mmssss to Mean Sidereal Time */
double ts2mst(); /* Calculate Mean Sidereal Time given Universal Time */
void compnut();	/* Compute nutation in longitude and obliquity and mean obliquity*/

#endif /* fitsfile_h_ */

/* May 31 1996	Use stream I/O for reading as well as writing
 * Jun 12 1996	Add byte-swapping subroutines
 * Jul 10 1996	FITS header now allocated in subroutines
 * Jul 17 1996	Add FITS table column extraction subroutines
 * Aug  6 1996	Add MOVEPIX, HDEL and HCHANGE declarations
 *
 * Oct 10 1997	FITS file opening subroutines now return int instead of FILE *
 *
 * May 27 1998	Split off fitsio and imhio subroutines to fitsio.h
 * Jun  4 1998	Change fits2iraf from int to int *
 * Jul 24 1998	Make IRAF header char instead of int
 * Aug 18 1998	Change name to fitsfile.h from fitsio.h
 * Oct  5 1998	Add isiraf() and isfits()
 * Oct  7 1998	Note separation of imhfile.c into two files
 *
 * Jul 15 1999	Add fileutil.c subroutines
 * Sep 28 1999	Add (1,1)-based image access subroutines
 * Oct 21 1999	Add fitswhead()
 * Nov  2 1999	Add date utilities from wcscat.h
 * Nov 23 1999	Add fitscimage()
 * Dec 15 1999	Fix misdeclaration of *2fd() subroutines, add fd2i(), dt2i()
 * Dec 20 1999	Add isdate()
 *
 * Jan 20 2000	Add conversions to and from Besselian and Julian epochs
 * Jan 21 2000	Add conversions to old FITS date and time
 * Jan 26 2000	Add conversion to modified Julian date (JD - 2400000.5
 * Mar 22 2000  Add lt2* and ut2* to get current time as local and UT
 * Mar 24 2000	Add tsi2* and tsu2* to convert IRAF and Unix seconds
 * Sep  8 2000	Improve comments
 *
 * Apr 24 2001	Add length of column name to column data structure
 * May 22 2001	Add day of year date conversion subroutines
 * Sep 25 2001	Add isfilelist() and isfile()
 *
 * Jan  8 2002	Add sts2c() and stc2s()
 * Apr  8 2002	Change all long declarations to time_t for compatibility
 * Jun 18 2002	Add fitserr() to print error messages
 * Aug 30 2002	Add Ephemeris Time date conversions
 * Sep 10 2002	Add Sidereal Time conversions
 * Oct 21 2002	Add fitsrsect() to read sections of FITS images
 *
 * Mar  5 2003	Add isimlistd() to check image lists with root directory
 */
