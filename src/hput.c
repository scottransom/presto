/*** File libwcs/hput.c
 *** January 4, 2002
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 1995-2002
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

 * Module:	hput.c (Put FITS Header parameter values)
 * Purpose:	Implant values for parameters into FITS header string
 * Subroutine:	hputi2 (hstring,keyword,ival) sets integer*2 ival
 * Subroutine:	hputi4 (hstring,keyword,ival) sets int ival
 * Subroutine:	hputr4 (hstring,keyword,rval) sets real*4 rval
 * Subroutine:	hputr8 (hstring,keyword,dval) sets real*8 dval
 * Subroutine:	hputnr8 (hstring,keyword,ndec,dval) sets real*8 dval
 * Subroutine:	hputra (hstring,keyword,lval) sets right ascension as string
 * Subroutine:	hputdec (hstring,keyword,lval) sets declination as string
 * Subroutine:	hputl  (hstring,keyword,lval) sets logical lval
 * Subroutine:	hputs  (hstring,keyword,cval) sets character string adding ''
 * Subroutine:	hputm  (hstring,keyword,cval) sets multi-line character string
 * Subroutine:	hputc  (hstring,keyword,cval) sets character string cval
 * Subroutine:	hdel   (hstring,keyword) deletes entry for keyword keyword
 * Subroutine:	hadd   (hplace,keyword) adds entry for keyword at hplace
 * Subroutine:	hchange (hstring,keyword1,keyword2) changes keyword for entry
 * Subroutine:	hputcom (hstring,keyword,comment) sets comment for parameter keyword
 * Subroutine:	ra2str (out, lstr, ra, ndec) converts RA from degrees to string
 * Subroutine:	dec2str (out, lstr, dec, ndec) converts Dec from degrees to string
 * Subroutine:	deg2str (out, lstr, deg, ndec) converts degrees to string
 * Subroutine:	num2str (out, num, field, ndec) converts number to string
 * Subroutine:  getltime () returns current local time as ISO-style string
 * Subroutine:  getutime () returns current UT as ISO-style string
 */
#include <sys/time.h>
#include <string.h>             /* NULL, strlen, strstr, strcpy */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitshead.h"

static int verbose = 0;         /* Set to 1 to print error messages and other info */


/*  HPUTI4 - Set int keyword = ival in FITS header string */

int hputi4(hstring, keyword, ival)

char *hstring;                  /* character string containing FITS-style header
                                   information in the format
                                   <keyword>= <value> {/ <comment>}
                                   each entry is padded with spaces to 80 characters */

char *keyword;                  /* character string containing the name of the variable
                                   to be returned.  hput searches for a line beginning
                                   with this string, and if there isn't one, creates one.
                                   The first 8 characters of keyword must be unique. */
int ival;                       /* int number */
{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf(value, "%d", ival);

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTI2 - Set short keyword = ival in FITS header string */

int hputi2(hstring, keyword, ival)

char *hstring;                  /* FITS header string */
char *keyword;                  /* Keyword name */
short ival;                     /* short number */

{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf(value, "%d", ival);

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTR4 - Set float keyword = rval in FITS header string */

int hputr4(hstring, keyword, rval)

char *hstring;                  /* FITS header string */
char *keyword;                  /* Keyword name */
float rval;                     /* float number */
{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf(value, "%f", rval);

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTR8 - Set double keyword = dval in FITS header string */

int hputr8(hstring, keyword, dval)

char *hstring;                  /* FITS header string */
char *keyword;                  /* Keyword name */
double dval;                    /* double number */
{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf(value, "%g", dval);

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTNR8 - Set double keyword = dval in FITS header string */

int hputnr8(hstring, keyword, ndec, dval)

char *hstring;                  /* FITS header string */
char *keyword;                  /* Keyword name */
int ndec;                       /* Number of decimal places to print */
double dval;                    /* double number */
{
    char value[30];
    char format[8];
    unsigned int i;

    /* Translate value from binary to ASCII */
    if (ndec < 0) {
        sprintf(format, "%%.%dg", -ndec);
        sprintf(value, format, dval);
        for (i = 0; i < strlen(value); i++)
            if (value[i] == 'e')
                value[i] = 'E';
    } else {
        sprintf(format, "%%.%df", ndec);
        sprintf(value, format, dval);
    }

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTRA - Set double keyword = hh:mm:ss.sss in FITS header string */

int hputra(hstring, keyword, ra)

char *hstring;                  /* FITS header string */
char *keyword;                  /* Keyword name */
double ra;                      /* Right ascension in degrees */
{
    char value[30];

    /* Translate value from binary to ASCII */
    ra2str(value, 30, ra, 3);

    /* Put value into header string */
    return (hputs(hstring, keyword, value));
}


/*  HPUTDEC - Set double keyword = dd:mm:ss.sss in FITS header string */

int hputdec(hstring, keyword, dec)

char *hstring;                  /* FITS header string */
char *keyword;                  /* Keyword name */
double dec;                     /* Declination in degrees */
{
    char value[30];

    /* Translate value from binary to ASCII */
    dec2str(value, 30, dec, 2);

    /* Put value into header string */
    return (hputs(hstring, keyword, value));
}



/*  HPUTL - Set keyword = F if lval=0, else T, in FITS header string */

int hputl(hstring, keyword, lval)

char *hstring;                  /* FITS header */
char *keyword;                  /* Keyword name */
int lval;                       /* logical variable (0=false, else true) */
{
    char value[8];

    /* Translate value from binary to ASCII */
    if (lval)
        strcpy(value, "T");
    else
        strcpy(value, "F");

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTM - Set multi-line character string in FITS header string */
/*          return number of keywords written */

int hputm(hstring, keyword, cval)

char *hstring;                  /* FITS header */
char *keyword;                  /* Keyword name root (6 characters or less) */
char *cval;                     /* character string containing the value for variable
                                   keyword.  trailing and leading blanks are removed.  */
{
    int lroot = 0, lcv, i, ii, nkw, lkw, lval;
    int comment = 0;
    char keyroot[8], newkey[12], *v, value[80];
    char squot = 39;

    /*  If COMMENT or HISTORY, use the same keyword on every line */
    lkw = strlen(keyword);
    if (lkw == 7 && (strncmp(keyword, "COMMENT", 7) == 0 ||
                     strncmp(keyword, "HISTORY", 7) == 0))
        comment = 1;

    /* Set up keyword root, shortening it to 6 characters, if necessary */
    else {
        comment = 0;
        strcpy(keyroot, keyword);
        lroot = strlen(keyroot);
        if (lroot > 6) {
            keyroot[6] = (char) 0;
            lroot = 6;
        }
    }

    /* Write keyword value one line of up to 67 characters at a time */
    ii = '1';
    nkw = 0;
    lcv = strlen(cval);
    strcpy(newkey, keyroot);
    strcat(newkey, "_");
    newkey[lroot + 2] = (char) 0;
    v = cval;
    while (lcv > 0) {
        if (lcv > 67)
            lval = 67;
        else
            lval = lcv;
        value[0] = squot;
        for (i = 1; i <= lval; i++)
            value[i] = *v++;

        /* Pad short strings to 8 characters */
        if (lval < 8) {
            for (i = lval + 1; i < 9; i++)
                value[i] = ' ';
            lval = 8;
        }
        value[lval + 1] = squot;
        value[lval + 2] = (char) 0;

        /* Add this line to the header */
        if (comment)
            i = hputc(hstring, keyroot, value);
        else {
            newkey[lroot + 1] = ii;
            ii++;
            i = hputc(hstring, newkey, value);
        }
        if (i != 0)
            return (i);
        nkw++;
        if (lcv > 67)
            lcv = lcv - 67;
        else
            break;
    }
    return (nkw);
}


/*  HPUTS - Set character string keyword = 'cval' in FITS header string */

int hputs(hstring, keyword, cval)

char *hstring;                  /* FITS header */
char *keyword;                  /* Keyword name */
char *cval;                     /* character string containing the value for variable
                                   keyword.  trailing and leading blanks are removed.  */
{
    char squot = 39;
    char value[80];
    int lcval, i, lkeyword;

    /*  If COMMENT or HISTORY, just add it as is */
    lkeyword = strlen(keyword);
    if (lkeyword == 7 && (strncmp(keyword, "COMMENT", 7) == 0 ||
                          strncmp(keyword, "HISTORY", 7) == 0))
        return (hputc(hstring, keyword, cval));

    /*  find length of variable string */
    lcval = strlen(cval);
    if (lcval > 67)
        lcval = 67;

    /* Put single quote at start of string */
    value[0] = squot;
    strncpy(&value[1], cval, lcval);

    /* If string is less than eight characters, pad it with spaces */
    if (lcval < 8) {
        for (i = lcval; i < 8; i++) {
            value[i + 1] = ' ';
        }
        lcval = 8;
    }

    /* Add single quote and null to end of string */
    value[lcval + 1] = squot;
    value[lcval + 2] = (char) 0;

    /* Put value into header string */
    return (hputc(hstring, keyword, value));
}


/*  HPUTC - Set character string keyword = value in FITS header string */
/*          Return -1 if error, 0 if OK */

int hputc(hstring, keyword, value)

char *hstring;
char *keyword;
char *value;                    /* character string containing the value for variable
                                   keyword.  trailing and leading blanks are removed.  */
{
    char squot = 39;
    char line[100];
    char newcom[50];
    char *v, *vp, *v1, *v2, *q1, *q2, *c1, *ve;
    int lkeyword, lcom, lval, lc, lv1, lhead;

    /* Find length of keyword, value, and header */
    lkeyword = strlen(keyword);
    lval = strlen(value);
    lhead = gethlength(hstring);

    /*  If COMMENT or HISTORY, always add it just before the END */
    if (lkeyword == 7 && (strncmp(keyword, "COMMENT", 7) == 0 ||
                          strncmp(keyword, "HISTORY", 7) == 0)) {

        /* First look for blank lines before END */
        v1 = blsearch(hstring, "END");

        /*  Otherwise, create a space for it at the end of the header */
        if (v1 == NULL) {

            /* Find end of header */
            v1 = ksearch(hstring, "END");
            v2 = v1 + 80;

            /* If header length is exceeded, return error code */
            if (v2 - hstring > lhead) {
                return (-1);
            }

            /* Move END down 80 characters */
            strncpy(v2, v1, 80);
        } else
            v2 = v1 + 80;

        /* Insert keyword */
        strncpy(v1, keyword, 7);

        /* Pad with spaces */
        for (vp = v1 + lkeyword; vp < v2; vp++)
            *vp = ' ';

        if (lval > 71)
            lv1 = 71;
        else
            lv1 = lval;

        /* Insert comment */
        strncpy(v1 + 9, value, lv1);
        return (0);
    }

    /* Otherwise search for keyword */
    else
        v1 = ksearch(hstring, keyword);

    /*  If parameter is not found, find a place to put it */
    if (v1 == NULL) {

        /* First look for blank lines before END */
        v1 = blsearch(hstring, "END");

        /*  Otherwise, create a space for it at the end of the header */
        if (v1 == NULL) {
            ve = ksearch(hstring, "END");
            v1 = ve;
            v2 = v1 + 80;

            /* If header length is exceeded, return error code */
            if (v2 - hstring > lhead) {
                return (-1);
            }

            strncpy(v2, ve, 80);
        } else
            v2 = v1 + 80;
        lcom = 0;
        newcom[0] = 0;
    }

    /*  Otherwise, extract the entry for this keyword from the header */
    else {
        strncpy(line, v1, 80);
        line[80] = 0;
        v2 = v1 + 80;

        /*  check for quoted value */
        q1 = strchr(line, squot);
        if (q1 != NULL)
            q2 = strchr(q1 + 1, squot);
        else
            q2 = line;

        /*  extract comment and remove trailing spaces */

        c1 = strchr(q2, '/');
        if (c1 != NULL) {
            lcom = 80 - (c1 - line);
            strncpy(newcom, c1 + 1, lcom);
            vp = newcom + lcom - 1;
            while (vp-- > newcom && *vp == ' ')
                *vp = 0;
            lcom = strlen(newcom);
        } else {
            newcom[0] = 0;
            lcom = 0;
        }
    }

    /* Fill new entry with spaces */
    for (vp = v1; vp < v2; vp++)
        *vp = ' ';

    /*  Copy keyword to new entry */
    strncpy(v1, keyword, lkeyword);

    /*  Add parameter value in the appropriate place */
    vp = v1 + 8;
    *vp = '=';
    vp = v1 + 9;
    *vp = ' ';
    vp = vp + 1;
    if (*value == squot) {
        strncpy(vp, value, lval);
        if (lval + 12 > 31)
            lc = lval + 12;
        else
            lc = 30;
    } else {
        vp = v1 + 30 - lval;
        strncpy(vp, value, lval);
        lc = 30;
    }

    /* Add comment in the appropriate place */
    if (lcom > 0) {
        if (lc + 2 + lcom > 80)
            lcom = 78 - lc;
        vp = v1 + lc + 2;       /* Jul 16 1997: was vp = v1 + lc * 2 */
        *vp = '/';
        vp = vp + 1;
        strncpy(vp, newcom, lcom);
        for (v = vp + lcom; v < v2; v++)
            *v = ' ';
    }

    if (verbose) {
        if (lcom > 0)
            fprintf(stderr, "HPUT: %s  = %s  / %s\n", keyword, value, newcom);
        else
            fprintf(stderr, "HPUT: %s  = %s\n", keyword, value);
    }

    return (0);
}


/*  HPUTCOM - Set comment for keyword or on line in FITS header string */

int hputcom(hstring, keyword, comment)

char *hstring;
char *keyword;
char *comment;
{
    char squot;
    char line[100];
    int lkeyword, lcom, lhead, i, lblank;
    char *vp, *v1, *v2, *c0 = NULL, *c1, *q1, *q2;
    char *ksearch();

    squot = 39;

    /*  Find length of variable name */
    lkeyword = strlen(keyword);
    lhead = gethlength(hstring);
    lcom = strlen(comment);

    /*  If COMMENT or HISTORY, always add it just before the END */
    if (lkeyword == 7 && (strncmp(keyword, "COMMENT", 7) == 0 ||
                          strncmp(keyword, "HISTORY", 7) == 0)) {

        /* Find end of header */
        v1 = ksearch(hstring, "END");
        v2 = v1 + 80;

        /* If header length is exceeded, return error code */
        if (v2 - hstring > lhead) {
            return (-1);
        }

        /* Move END down 80 characters */
        strncpy(v2, v1, 80);

        /*  blank out new line and insert keyword */
        for (vp = v1; vp < v2; vp++)
            *vp = ' ';
        strncpy(v1, keyword, lkeyword);
    }

    /* Search header string for variable name */
    else {
        v1 = ksearch(hstring, keyword);
        v2 = v1 + 80;

        /* If parameter is not found, return without doing anything */
        if (v1 == NULL) {
            if (verbose)
                fprintf(stderr, "HPUTCOM: %s not found\n", keyword);
            return (-1);
        }

        /* Otherwise, extract entry for this variable from the header */
        strncpy(line, v1, 80);
        line[80] = '\0';        /* Null-terminate line before strchr call */

        /* check for quoted value */
        q1 = strchr(line, squot);
        if (q1 != NULL)
            q2 = strchr(q1 + 1, squot);
        else
            q2 = NULL;

        if (q2 == NULL || q2 - line < 31)
            c0 = v1 + 31;
        else
            c0 = v1 + (q2 - line) + 2;  /* allan: 1997-09-30, was c0=q2+2 */

        /* If comment will not fit, do not add it */
        if (c0 - v1 > 77)
            return (-1);
        strncpy(c0, "/ ", 2);
    }

    /* Create new entry */
    if (lcom > 0) {
        c1 = c0 + 2;
        if (c1 + lcom > v2)
            lcom = v2 - c1 - 2;
        strncpy(c1, comment, lcom);
        lblank = v2 - c1 - lcom;
        c1 = c1 + lcom;
        for (i = 0; i < lblank; i++)
            *c1++ = ' ';
    }

    if (verbose) {
        fprintf(stderr, "HPUTCOM: %s / %s\n", keyword, comment);
    }
    return (0);
}


/*  HDEL - Set character string keyword = value in FITS header string
 *	    returns 1 if entry deleted, else 0
 */

int hdel(hstring, keyword)

char *hstring;                  /* FITS header */
char *keyword;                  /* Keyword of entry to be deleted */
{
    char *v, *v1, *v2, *ve;
    char *ksearch();

    /* Search for keyword */
    v1 = ksearch(hstring, keyword);

    /*  If keyword is not found, return header unchanged */
    if (v1 == NULL) {
        return (0);
    }

    /*  Find end of header */
    ve = ksearch(hstring, "END");

    /* Shift rest of header up one line */
    for (v = v1; v < ve; v = v + 80) {
        v2 = v + 80;
        strncpy(v, v2, 80);
    }

    /* Cover former last line with spaces */
    v2 = ve + 80;
    for (v = ve; v < v2; v++)
        *v = ' ';

    return (1);
}


/*  HADD - Add character string keyword = value to FITS header string
 *	    returns 1 if entry added, else 0
 *	    Call hputx() to put value into entry
 */

int hadd(hplace, keyword)

char *hplace;                   /* FITS header position for new keyword */
char *keyword;                  /* Keyword of entry to be deleted */
{
    char *v, *v1, *v2, *ve;
    char *ksearch();
    int i, lkey;

    /*  Find end of header */
    ve = ksearch(hplace, "END");

    /*  If END is not found, return header unchanged */
    if (ve == NULL) {
        return (0);
    }

    v1 = hplace;

    /* Shift rest of header down one line */
    /* limit bug found by Paolo Montegriffo fixed 2000-04-19 */
    for (v = ve; v >= v1; v = v - 80) {
        v2 = v + 80;
        strncpy(v2, v, 80);
    }

    /* Cover former first line with new keyword */
    lkey = strlen(keyword);
    strncpy(hplace, keyword, lkey);
    if (lkey < 8) {
        for (i = lkey; i < 8; i++)
            hplace[i] = ' ';
        hplace[8] = '=';
    }
    for (i = 9; i < 80; i++)
        hplace[i] = ' ';

    return (1);
}


/*  HCHANGE - Changes keyword for entry from keyword1 to keyword2 in FITS
              header string
 *	      returns 1 if entry changed, else 0
 */

int hchange(hstring, keyword1, keyword2)

char *hstring;                  /* FITS header */
char *keyword1;                 /* Keyword to be changed */
char *keyword2;                 /* New keyword name */
{
    char *v, *v1, *v2;
    int lv2, i;
    char *ksearch();

    /* Search for keyword */
    v1 = ksearch(hstring, keyword1);

    /*  If keyword is not found, return header unchanged */
    if (!v1)
        return (0);

    else {
        lv2 = strlen(keyword2);
        v = v1;
        v2 = keyword2;
        for (i = 0; i < 8; i++) {
            if (i < lv2)
                v[i] = v2[i];
            else
                v[i] = ' ';
        }
    }

    return (1);
}


/* Write the right ascension ra in sexagesimal format into string*/

void ra2str(string, lstr, ra, ndec)

char *string;                   /* Character string (returned) */
int lstr;                       /* Maximum number of characters in string */
double ra;                      /* Right ascension in degrees */
int ndec;                       /* Number of decimal places in seconds */

{
    double a, b;
    double seconds;
    char tstring[64];
    int hours;
    int minutes;
    int isec;
    double dsgn;

    /* Keep RA between 0 and 360 */
    if (ra < 0.0) {
        ra = -ra;
        dsgn = -1.0;
    } else
        dsgn = 1.0;
    ra = fmod(ra, 360.0);
    ra *= dsgn;
    if (ra < 0.0)
        ra = ra + 360.0;

    a = ra / 15.0;

    /* Convert to hours */
    hours = (int) a;

    /* Compute minutes */
    b = (a - (double) hours) * 60.0;
    minutes = (int) b;

    /* Compute seconds */
    seconds = (b - (double) minutes) * 60.0;
    isec = (int) (seconds + 0.5);

    if (ndec > 5) {
        if (seconds > 59.999999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%09.6f", hours, minutes, seconds);
    } else if (ndec > 4) {
        if (seconds > 59.99999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%08.5f", hours, minutes, seconds);
    } else if (ndec > 3) {
        if (seconds > 59.9999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%07.4f", hours, minutes, seconds);
    } else if (ndec > 2) {
        if (seconds > 59.999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%06.3f", hours, minutes, seconds);
    } else if (ndec > 1) {
        if (seconds > 59.99) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%05.2f", hours, minutes, seconds);
    } else if (ndec > 0) {
        if (seconds > 59.9) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%04.1f", hours, minutes, seconds);
    } else if (ndec > -1) {
        if (isec > 59) {
            isec = 0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            hours = hours + 1;
        }
        hours = hours % 24;
        (void) sprintf(tstring, "%02d:%02d:%04.1f", hours, minutes, seconds);
    }

    /* Move formatted string to returned string */
    if (strlen(tstring) < (unsigned int) (lstr - 1))
        strcpy(string, tstring);
    else {
        strncpy(string, tstring, lstr - 1);
        string[lstr - 1] = 0;
    }
    return;
}


/* Write the variable a in sexagesimal format into string */

void dec2str(string, lstr, dec, ndec)

char *string;                   /* Character string (returned) */
int lstr;                       /* Maximum number of characters in string */
double dec;                     /* Declination in degrees */
int ndec;                       /* Number of decimal places in arcseconds */

{
    double a, b, dsgn, deg1;
    double seconds;
    char sign;
    int degrees;
    int minutes;
    int isec;
    char tstring[64];

    /* Keep angle between -180 and 360 degrees */
    deg1 = dec;
    if (deg1 < 0.0) {
        deg1 = -deg1;
        dsgn = -1.0;
    } else
        dsgn = 1.0;
    deg1 = fmod(deg1, 360.0);
    deg1 *= dsgn;
    if (deg1 <= -180.0)
        deg1 = deg1 + 360.0;

    a = deg1;

    /* Set sign and do all the rest with a positive */
    if (a < 0) {
        sign = '-';
        a = -a;
    } else
        sign = '+';

    /* Convert to degrees */
    degrees = (int) a;

    /* Compute minutes */
    b = (a - (double) degrees) * 60.0;
    minutes = (int) b;

    /* Compute seconds */
    seconds = (b - (double) minutes) * 60.0;
    isec = (int) (seconds + 0.5);

    if (ndec > 5) {
        if (seconds > 59.999999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%09.6f", sign, degrees, minutes,
                       seconds);
    } else if (ndec > 4) {
        if (seconds > 59.99999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%08.5f", sign, degrees, minutes,
                       seconds);
    } else if (ndec > 3) {
        if (seconds > 59.9999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%07.4f", sign, degrees, minutes,
                       seconds);
    } else if (ndec > 2) {
        if (seconds > 59.999) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%06.3f", sign, degrees, minutes,
                       seconds);
    } else if (ndec > 1) {
        if (seconds > 59.99) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%05.2f", sign, degrees, minutes,
                       seconds);
    } else if (ndec > 0) {
        if (seconds > 59.9) {
            seconds = 0.0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%04.1f", sign, degrees, minutes,
                       seconds);
    } else if (ndec > -1) {
        if (isec > 59) {
            isec = 0;
            minutes = minutes + 1;
        }
        if (minutes > 59) {
            minutes = 0;
            degrees = degrees + 1;
        }
        (void) sprintf(tstring, "%c%02d:%02d:%04.1f", sign, degrees, minutes,
                       seconds);
    }

    /* Move formatted string to returned string */
    if (strlen(tstring) < (unsigned int) (lstr - 1))
        strcpy(string, tstring);
    else {
        strncpy(string, tstring, lstr - 1);
        string[lstr - 1] = 0;
    }
    return;
}


/* Write the angle a in decimal format into string */

void deg2str(string, lstr, deg, ndec)

char *string;                   /* Character string (returned) */
int lstr;                       /* Maximum number of characters in string */
double deg;                     /* Angle in degrees */
int ndec;                       /* Number of decimal places in degree string */

{
    char degform[8];
    int field;
    char tstring[64];
    double deg1;
    double dsgn;

    /* Keep angle between -180 and 360 degrees */
    deg1 = deg;
    if (deg1 < 0.0) {
        deg1 = -deg1;
        dsgn = -1.0;
    } else
        dsgn = 1.0;
    deg1 = fmod(deg1, 360.0);
    deg1 *= dsgn;
    if (deg1 <= -180.0)
        deg1 = deg1 + 360.0;

    /* Write angle to string, adding 4 digits to number of decimal places */
    field = ndec + 4;
    if (ndec > 0) {
        sprintf(degform, "%%%d.%df", field, ndec);
        sprintf(tstring, degform, deg1);
    } else {
        sprintf(degform, "%%%4d", field);
        sprintf(tstring, degform, (int) deg1);
    }

    /* Move formatted string to returned string */
    if (strlen(tstring) < (unsigned int) (lstr - 1))
        strcpy(string, tstring);
    else {
        strncpy(string, tstring, (unsigned int) (lstr - 1));
        string[lstr - 1] = 0;
    }
    return;
}


/* Write the variable a in decimal format into field-character string  */

void num2str(string, num, field, ndec)

char *string;                   /* Character string (returned) */
double num;                     /* Number */
int field;                      /* Number of characters in output field (0=any) */
int ndec;                       /* Number of decimal places in degree string */

{
    char numform[8];

    if (field > 0) {
        if (ndec > 0) {
            sprintf(numform, "%%%d.%df", field, ndec);
            sprintf(string, numform, num);
        } else {
            sprintf(numform, "%%%dd", field);
            sprintf(string, numform, (int) num);
        }
    } else {
        if (ndec > 0) {
            sprintf(numform, "%%.%df", ndec);
            sprintf(string, numform, num);
        } else {
            sprintf(string, "%d", (int) num);
        }
    }
    return;
}

/* Dec 14 1995	Original subroutines

 * Feb  5 1996	Added HDEL to delete keyword entry from FITS header
 * Feb  7 1996	Add EOS to LINE in HPUTC
 * Feb 21 1996	Add RA2STR and DEC2STR string routines
 * Jul 19 1996	Add HPUTRA and HPUTDEC
 * Jul 22 1996	Add HCHANGE to change keywords
 * Aug  5 1996	Add HPUTNR8 to save specific number of decimal places
 * Oct 15 1996	Fix spelling
 * Nov  1 1996	Add DEG2STR to set specific number of decimal places
 * Nov  1 1996	Allow DEC2STR to handle upt to 6 decimal places
 *
 * Mar 20 1997	Fix format error in DEG2STR
 * Jul  7 1997	Fix 2 errors in HPUTCOM found by Allan Brighton
 * Jul 16 1997	Fix error in HPUTC found by Allan Brighton
 * Jul 17 1997	Fix error in HPUTC found by Allan Brighton
 * Sep 30 1997	Fix error in HPUTCOM found by Allan Brighton
 * Dec 15 1997	Fix minor bugs after lint
 * Dec 31 1997	Always put two hour digits in RA2STR
 *
 * Feb 25 1998	Add HADD to insert keywords at specific locations
 * Mar 27 1998	If n is negative, write g format in HPUTNR8()
 * Apr 24 1998	Add NUM2STR() for easy output formatting
 * Apr 30 1998	Use BLSEARCH() to overwrite blank lines before END
 * May 27 1998	Keep Dec between -90 and +90 in DEC2STR()
 * May 28 1998	Keep RA between 0 and 360 in RA2STR()
 * Jun  2 1998	Fix bug when filling in blank lines before END
 * Jun 24 1998	Add string length to ra2str(), dec2str(), and deg2str()
 * Jun 25 1998	Make string converstion subroutines more robust
 * Aug 31 1998	Add getltime() and getutime()
 * Sep 28 1998	Null-terminate comment in HPUTCOM (Allan Brighton)
 * Oct  1 1998	Change clock declaration in getltime() from int (Allan Brighton)
 *
 * Jan 28 1999	Fix bug to avoid writing HISTORY or COMMENT past 80 characters
 * Jul 14 1999	Pad string in hputs() to minimum of 8 characters
 * Aug 16 1999	Keep angle between -180 and +360 in dec2str()
 * Oct  6 1999	Reallocate header buffer if it is too small in hputc()
 * Oct 14 1999	Do not reallocate header; return error if not successful
 *
 * Mar  2 2000	Do not add quotes if adding HISTORY or COMMENT with hputs()
 * Mar 22 2000	Move getutime() and getltime() to dateutil.c
 * Mar 27 2000	Add hputm() for muti-line keywords
 * Mar 27 2000	Fix bug testing for space to fit comment in hputcom()
 * Apr 19 2000	Fix bug in hadd() which overwrote line
 * Jun  2 2000	Dropped unused variable lv in hputm() after lint
 * Jul 20 2000	Drop unused variables blank and i in hputc()
 *
 * Jan 11 2001	Print all messages to stderr
 * Jan 18 2001	Drop declaration of blsearch(); it is in fitshead.h
 *
 * Jan  4 2002	Fix placement of comments
 */
