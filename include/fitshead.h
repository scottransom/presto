/*** File fitshead.h  FITS header access subroutines
 *** August 30, 2002
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 1996-2002
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

/* Declarations for subroutine sin hget.c, hput.c, and iget.c */

#ifndef _fitshead_h_
#define _fitshead_h_

#include <sys/types.h>

#ifdef __cplusplus /* C++ prototypes */
extern "C" {

/* Subroutines in hget.c */
    int hgeti2(			/* Extract short value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	short* val);		/* short integer value (returned) */
    int hgeti4c(		/* Extract int value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	const char mchar,	/* WCS to use */
	int* val);		/* integer value (returned) */
    int hgeti4(			/* Extract int value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	int* val);		/* integer value (returned) */
    int hgetr4(			/* Extract float value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	float* val);		/* float value (returned) */
    int hgetr8c(			/* Extract double value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	const char mchar,	/* WCS to use */
	double* val);		/* double value (returned) */
    int hgetr8(			/* Extract double value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	double* val);		/* double value (returned) */
    int hgetra(			/* Extract right ascension from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	double* ra);		/* RA in degrees (returned) */
    int hgetdec(		/* Extract declination from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	double* dec);		/* Dec in degrees (returned) */
    int hgetdate(		/* Extract date from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	double* date);		/* Date in fractional years (returned) */
    int hgetl(			/* Extract boolean value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	int* lval);		/* 1 if T, 0 if F (returned) */
    int hgetsc(			/* Extract string value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	const char mchar,	/* WCS to use */
	const int lstr,		/* maximum length of returned string */
	char* string);		/* null-terminated string value (returned) */
    int hgets(			/* Extract string value from FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	const int lstr,		/* maximum length of returned string */
	char* string);		/* null-terminated string value (returned) */
    int hgetm (			/* Extract string from multiple keywords */
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	const int lstr,		/* maximum length of returned string */
	char* string);		/* null-terminated string value (returned) */
    int hgetndec(		/* Find number of decimal places in FITS value*/
	const char* hstring,	/* FITS header string */
	const char* keyword,	/* FITS keyword */
	int* ndec);		/* number of decimal places (returned) */

    char* hgetc(		/* Return pointer to value for FITS keyword */
	const char* hstring,	/* FITS header string */
	const char* keyword);	/* FITS keyword */

    char* ksearch(		/* Return pointer to keyword in FITS header */
	const char* hstring,	/* FITS header string */
	const char* keyword);	/* FITS keyword */
    char *blsearch (
	const char* hstring,	/* FITS header string */
	const char* keyword);	/* FITS keyword */

    char *strsrch (		/* Find string s2 within string s1 */
	const char* s1,		/* String to search */
	const char* s2);	/* String to look for */
    char *strnsrch (		/* Find string s2 within string s1 */
	const char* s1,		/* String to search */
	const char* s2,		/* String to look for */
	const int ls1);		/* Length of string being searched */

    char *strcsrch (		/* Find string s2 within string s1 (no case) */
	const char* s1,		/* String to search */
	const char* s2);	/* String to look for */
    char *strncsrch (		/* Find string s2 within string s1 (no case) */
	const char* s1,		/* String to search */
	const char* s2,		/* String to look for */
	const int ls1);		/* Length of string being searched */

    int hlength(		/* Set length of unterminated FITS header */
        char    *header,        /* FITS header */
        const int lhead);	/* Allocated length of FITS header */
    int gethlength(		/* Get length of current FITS header */
        char    *header);	/* FITS header */

    double str2ra(		/* Return RA in degrees from string */
	const char* in);	/* Character string (hh:mm:ss.sss or dd.dddd) */
    double str2dec(		/* Return Dec in degrees from string */
	const char* in);	/* Character string (dd:mm:ss.sss or dd.dddd) */

    int isnum(			/* Return 1 if number, else 0 */
	const char* string);	/* Character string which may be a number */
    int notnum(			/* Return 0 if number, else 1 */
	const char* string);	/* Character string which may be a number */

    char *getltime();		/* Return current local time in ISO format */
    char *getutime();		/* Return current UT as an ISO-format string */

/* Subroutines in iget.c */
    int mgets(			/* Extract string from multiline FITS keyword */
	const char* hstring,	/* FITS header string */
	const char* mkey,	/* FITS keyword root _n added for extra lines */
	const char* keyword,	/* IRAF keyword */
	const int lstr,		/* maximum length of returned string */
	char* string);		/* null-terminated string value (returned) */
    int mgeti4(			/* Extract int from multiline FITS keyword */
	const char* hstring,	/* FITS header string */
	const char* mkey,	/* FITS keyword root _n added for extra lines */
	const char* keyword,	/* IRAF keyword */
	int* ival);		/* int keyword value (returned) */
    int mgetr8(			/* Extract double from multiline FITS keyword */
	const char* hstring,	/* FITS header string */
	const char* mkey,	/* FITS keyword root _n added for extra lines */
	const char* keyword,	/* IRAF keyword */
	double dval);		/* double keyword value (returned) */
    int igeti4(			/* Extract int from IRAF keyword string */
	const char* hstring,	/* Multiline IRAF keyword string value */
	const char* keyword,	/* IRAF keyword */
	int* val);		/* int value (returned) */
    int igetr4(			/* Extract float from IRAF keyword string */
	const char* hstring,	/* Multiline IRAF keyword string value */
	const char* keyword,	/* IRAF keyword */
	float* val);		/* float value (returned) */
    int igetr8(			/* Extract double from IRAF keyword string */
	const char* hstring,	/* Multiline IRAF keyword string value */
	const char* keyword,	/* IRAF keyword */
	double* val);		/* double value (returned) */
    int igets(			/* Extract string from IRAF keyword string */
	const char* hstring,	/* Multiline IRAF keyword string value */
	const char* keyword,	/* IRAF keyword */
	const int lstr,		/* maximum length of returned string */
	char* string);		/* null-terminated string value (returned) */

/* Subroutines in hput.c */
/* All hput* routines return 0 if successful, else -1 */
    int hputi2(		/* Implant short value into FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const short ival);	/* short value */
    int hputi4(		/* Implant int value into FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const int ival);	/* int value */
    int hputr4(		/* Implant float value into FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const float rval);	/* float value */
    int hputr8(		/* Implant short into FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const double dval);	/* double value */
    int hputnr8(	/* double with specified number of decimal places */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const int ndec,		/* Number of decimal places in keyword value */
	const double dval);	/* double value */
    int hputs(			/* Quoted character string into FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const char* cval);	/* Character string value */
    int hputm(			/* Quoted character string, mutiple keywords */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const char* cval);	/* Character string value */
    int hputcom(		/* Add comment to keyword line in FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const char* comment);	/* Comment string */
    int hputra(	/* Right ascension in degrees into hh:mm:ss.sss */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const double ra);	/* Right ascension in degrees */
    int hputdec(		/* Declination in degrees into dd:mm:ss.ss */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const double dec);	/* Declination in degrees */
    int hputl(			/* Implant boolean value into FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const int lval);	/* 0->F, else ->T */
    int hputc(			/* Implant character string without quotes */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword,	/* FITS keyword */
	const char* cval);	/* Character string value */

    int hdel(			/* Delete a keyword line from a FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword);	/* FITS keyword to delete */
    int hadd(			/* Add a keyword line from a FITS header */
	char* hplace,		/* Location in FITS header string (modified) */
	const char* keyword);	/* FITS keyword to add */
    int hchange(		/* Change a keyword name in a FITS header */
	char* hstring,		/* FITS header string (modified) */
	const char* keyword1,	/* Current FITS keyword name */
	const char* keyword2);	/* New FITS keyword name */

    void ra2str(		/* Convert degrees to hh:mm:ss.ss */
        char *string,		/* Character string (returned) */
	int lstr,		/* Length of string */
        const double ra,	/* Right ascension in degrees */
        const int ndec);	/* Number of decimal places in seconds */
    void dec2str(		/* Convert degrees to dd:mm:ss.ss */
        char *string,		/* Character string (returned) */
	int lstr,		/* Length of string */
        const double dec,	/* Declination in degrees */
        const int ndec);	/* Number of decimal places in arcseconds */
    void deg2str(		/* Format angle into decimal degrees string */
        char *string,		/* Character string (returned) */
	int lstr,		/* Length of string */
        const double deg,	/* Angle in degrees */
        const int ndec);	/* Number of decimal places in degrees */
    void num2str(		/* Format number into string */
        char *string,		/* Character string (returned) */
        const double  num,	/* Number */
	const int field,	/* Total field size in characters */
        const int ndec);	/* Number of decimal places */
};
#else /* __cplusplus */

/* Subroutines in hget.c */

/* Extract a value from a FITS header for given keyword */
extern int hgeti4();	/* int (Multiple WCS) */
extern int hgeti4c();	/* int */
extern int hgeti2();	/* short */
extern int hgetr4();	/* float */
extern int hgetr8();	/* double */
extern int hgetr8c();	/* double (Multiple WCS) */
extern int hgetra();	/* Right ascension in degrees from string */
extern int hgetdec();	/* Declination in degrees from string */
extern int hgetdate();	/* Date in years from FITS date string */
extern int hgetl();	/* T->1, F->0 from FITS logical entry */
extern int hgets();	/* Previously allocated string */
extern int hgetsc();	/* Previously allocated string (Multiple WCS) */
extern int hgetm();	/* Previously allocated string from multiple keywords */
extern char *hgetc();	/* Return pointer to string */
extern int hgetndec();	/* Number of decimal places in keyword value */

/* Subroutines to convert strings to RA and Dec in degrees */
extern double str2ra();
extern double str2dec();

/* Check to see whether a string is a number or not */
extern int isnum();
extern int notnum();

/* Find given keyword entry in FITS header */
extern char *ksearch();

/* Find beginning of fillable blank line before FITS header keyword */
extern char *blsearch();

/* Search for substring s2 within string s1 */
extern char *strsrch ();	/* s1 null-terminated */
extern char *strnsrch ();	/* s1 ls1 characters long */
extern char *strcsrch ();	/* s1 null-terminated (case-insensitive) */
extern char *strncsrch ();	/* s1 ls1 characters long (case-insensitive) */

/* Set length of header which is not null-terminated */
extern int hlength();

/* Get length of current FITS header */
extern int gethlength();

/* Subroutines in iget.c */
extern int mgets();	/* Previously allocated string from multiline keyword */
extern int mgetr8();	/* double from multiline keyword */
extern int mgeti4();	/* int from multiline keyword */
extern int igeti4();	/* long integer from IRAF compound keyword value */
extern int igetr4();	/* real from IRAF compound keyword value */
extern int igetr8();	/* double from IRAF compound keyword value */
extern int igets();	/* character string from IRAF compound keyword value */

/* Subroutines in hput.c */

/* Implant a value into a FITS header for given keyword */
extern int hputi4();	/* int */
extern int hputi2();	/* short */
extern int hputr4();	/* float */
extern int hputr8();	/* double */
extern int hputnr8();	/* double with specified number of decimal places */
extern int hputra();	/* Right ascension in degrees into hh:mm:ss.sss */
extern int hputdec();	/* Declination in degrees into dd:mm:ss.ss */
extern int hputl();	/* 0 -> F, else T FITS logical entry */
extern int hputs();	/* Quoted character string */
extern int hputm();	/* Quoted character string into mutiple keywords */
extern int hputc();	/* Character string without quotes (returns 0 if OK) */
extern int hputcom();	/* Comment after keyword=value (returns 0 if OK) */

extern int hdel();	/* Delete a keyword line from a FITS header */
extern int hadd();	/* Add a keyword line to a FITS header */
extern int hchange();	/* Change a keyword name in a FITS header */

/* Subroutines to convert RA and Dec in degrees to strings */
extern void ra2str();
extern void dec2str();

extern void deg2str();
extern void num2str();

extern char *getltime(); /* Return current local time in ISO format */
extern char *getutime(); /* Return current UT as an ISO-format string */

#endif	/* __cplusplus */
#endif	/* fitshead_h_ */

/* Apr 26 1996	Add HGETDATE to get year from date string
 * May 22 1996	Return double from STR2RA and STR2DEC
 * May 31 1996	Use stream I/O for reading as well as writing
 * Jun 12 1996	Add byte-swapping subroutines
 * Jul 10 1996	FITS header now allocated in subroutines
 * Jul 17 1996	Add FITS table column extraction subroutines
 * Jul 19 1996	Add declarations for header implanting subroutines
 * Aug  5 1996	Add HLENGTH for FITS headers which are not null-terminated
 * Aug  5 1996	Add STRNSRCH for FITS headers which are not null-terminated
 * Aug  6 1996	Add HPUTNR8 to save a specified number of decimal places
 * Aug  6 1996	Add MOVEPIX, HDEL and HCHANGE declarations
 * Nov  1 1996	Add DEG2STR
 * Dec 12 1996	Add ISNUM
 *
 * Oct 10 1997	FITS file opening subroutines now return int instead of FILE *
 *
 * Mar 12 1998	Add NOTNUM
 * Apr 30 1998	Clean up declarations and add more comments
 * May 12 1998	Add MGETS, MGETR8, MGETI4 for IRAF multi-line keywords
 * May 26 1998	Add HGETNDEC for number of decimal places in keyword value
 * May 27 1998	Add BLSEARCH to find usable blank lines in header
 * May 27 1998	Split off fitsio and imhio subroutines to fitsio.h
 * May 27 1998	Add all subroutines in hget.c, hput.c, and iget.c to C++ dec.
 * Jun 24 1998	Add string lengths to ra2str(), dec2str, and deg2str() calls
 * Jun 25 1998	Fix other C++ declarations with added string lengths
 * Aug 31 1998	Add current date subroutines getltime() and getutime()
 * Oct 28 1998	Add missing hgetc() to non c++ declarations
 *
 * Oct  6 1999	Add gethlength() to return current size of header
 * Oct 14 1999	All HPUT subroutines now return an error code, 0 if OK, else -1
 * Oct 15 1999	Add hputcom() declaration
 * Oct 21 1999	Add hgetm() declaration
 *
 * Mar 22 2000	Add int to iget*() declarations
 * Mar 27 2000	Add hputm() declaration
 *
 * Apr  3 2002	Add hgeti4c(), hgetr8c(), and hgetsc()
 * Apr  8 2002	Include sys/types.h
 * Aug 30 2002	Add strcsrch() and strncsrch()
 */
