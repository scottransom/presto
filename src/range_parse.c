/*
  ranges_to_ivect() is Copyright 2004 by Scott Ransom
  Also released under the GPL.
*/

/*
 * Copyright (c) 2000 Silicon Graphics, Inc.  All Rights Reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of version 2 of the GNU General Public License as
 * published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it would be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * Further, this software is distributed without any warranty that it is
 * free of the rightful claim of any third person regarding infringement
 * or the like.  Any license provided herein, whether implied or
 * otherwise, applies only to this software file.  Patent licenses, if
 * any, provided herein do not apply to combinations of this program with
 * other software, or any other product whatsoever.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write the Free Software Foundation, Inc., 59
 * Temple Place - Suite 330, Boston MA 02111-1307, USA.
 * 
 * Contact information: Silicon Graphics, Inc., 1600 Amphitheatre Pkwy,
 * Mountain View, CA  94043, or:
 * 
 * http://www.sgi.com 
 * 
 * For further information regarding this notice, see: 
 * 
 * http://oss.sgi.com/projects/GenInfo/NoticeExplan/
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <malloc.h>

/* #define DEBUGPRINT */

/*
 * Internal format of the range array set up by parse_range()
 */

struct range {
    int min;
    int max;
    int mult;
};

/*
 * parse_ranges() is a function to parse a comma-separated list of range
 * tokens each having the following form:
 *
 *		num
 *	or
 *		min:max[:mult]
 *
 * any of the values may be blank (ie. min::mult, :max, etc.) and default
 * values for missing arguments may be supplied by the caller.
 *
 * The special first form is short hand for 'num:num'.
 *
 * After parsing the string, the ranges are put into an array of integers,
 * which is malloc'd by the routine.  The min, max, and mult entries of each
 * range can be extracted from the array using the range_min(), range_max(),
 * and range_mult() functions.
 *
 * It is the responsibility of the caller to free the space allocated by
 * parse_ranges() - a single call to free() will free the space.
 *
 *	str		The string to parse - assumed to be a comma-separated
 *			list of tokens having the above format.
 *	defmin		default value to plug in for min, if it is missing
 *	defmax		default value to plug in for max, if it is missing     
 *	defmult		default value to plug in for mult, if missing
 *	parse_func	A user-supplied function pointer, which parse_ranges()
 *			can call to parse the min, max, and mult strings.  This
 *			allows for customized number formats.  The function
 *			MUST have the following prototype:
 *				parse_func(char *str, int *val)
 *			The function should return -1 if str cannot be parsed
 *			into an integer, or >= 0 if it was successfully 
 *			parsed.  The resulting integer will be stored in
 *			*val.  If parse_func is NULL, parse_ranges will parse
 *			the tokens in a manner consistent with the the sscanf
 *			%i format.
 *	range_ptr	A user-supplied char **, which will be set to point
 *			at malloc'd space which holds the parsed range
 *			values.   If range_ptr is NULL, parse_ranges() just
 *			parses the string.  The data returned in range_ptr
 *			should not be processed directly - use the functions
 *			range_min(), range_max(), and range_mult() to access
 *			data for a given range.
 *	errptr		user-supplied char ** which can be set to point to a
 *			static error string.  If errptr is NULL, it is ignored.
 *
 * parse_range() returns -1 on error, or the number of ranges parsed.
 */

static int str_to_int();

int parse_ranges(str, defmin, defmax, defmult, parse_func, rangeptr, errptr)
char *str;
int defmin;
int defmax;
int defmult;
int (*parse_func) ();
char **rangeptr;
char **errptr;
{
    int ncommas;
    char *tmpstr, *cp, *tok, *n1str, *n2str, *multstr;
    struct range *rp, *ranges;
    static char errmsg[256];

    if (errptr != NULL) {
        *errptr = errmsg;
    }

    for (ncommas = 0, cp = str; *cp != '\0'; cp++) {
        if (*cp == ',') {
            ncommas++;
        }
    }

    if (parse_func == NULL) {
        parse_func = str_to_int;
    }

    tmpstr = strdup(str);
    ranges = (struct range *) malloc((ncommas + 1) * sizeof(struct range));
    rp = ranges;

    tok = strtok(tmpstr, ",");
    while (tok != NULL) {
        n1str = tok;
        n2str = NULL;
        multstr = NULL;

        rp->min = defmin;
        rp->max = defmax;
        rp->mult = defmult;

        if ((cp = strchr(n1str, ':')) != NULL) {
            *cp = '\0';
            n2str = cp + 1;

            if ((cp = strchr(n2str, ':')) != NULL) {
                *cp = '\0';
                multstr = cp + 1;
            }
        }

        /*
         * Parse the 'min' field - if it is zero length (:n2[:mult]
         * format), retain the default value, otherwise, pass the
         * string to the parse function.
         */

        if ((int) strlen(n1str) > 0) {
            if ((*parse_func) (n1str, &rp->min) < 0) {
                sprintf(errmsg, "error parsing string %s into an integer", n1str);
                free(tmpstr);
                free(ranges);
                return -1;
            }
        }

        /*
         * Process the 'max' field - if one was not present (n1 format)
         * set max equal to min.  If the field was present, but 
         * zero length (n1: format), retain the default.  Otherwise
         * pass the string to the parse function.
         */

        if (n2str == NULL) {
            rp->max = rp->min;
        } else if ((int) strlen(n2str) > 0) {
            if ((*parse_func) (n2str, &rp->max) < 0) {
                sprintf(errmsg, "error parsing string %s into an integer", n2str);
                free(tmpstr);
                free(ranges);
                return -1;
            }
        }

        /*
         * Process the 'mult' field - if one was not present 
         * (n1:n2 format), or the field was zero length (n1:n2: format)
         * then set the mult field to defmult - otherwise pass then
         * mult field to the parse function.
         */

        if (multstr != NULL && (int) strlen(multstr) > 0) {
            if ((*parse_func) (multstr, &rp->mult) < 0) {
                sprintf(errmsg, "error parsing string %s into an integer", multstr);
                free(tmpstr);
                free(ranges);
                return -1;
            }
        }

        rp++;
        tok = strtok(NULL, ",");
    }

    free(tmpstr);

    if (rangeptr != NULL) {
        *rangeptr = (char *) ranges;
    } else {
        free(ranges);           /* just running in parse mode */
    }

    return (rp - ranges);
}

/*
 * The default integer-parsing function
 */

static int str_to_int(str, ip)
char *str;
int *ip;
{
    char c;

    if (sscanf(str, "%i%c", ip, &c) != 1) {
        return -1;
    } else {
        return 0;
    }
}

/*
 * Three simple functions to return the min, max, and mult values for a given
 * range.  It is assumed that rbuf is a range buffer set up by parse_ranges(),
 * and that r is a valid range within that buffer.
 */

int range_min(rbuf, r)
char *rbuf;
int r;
{
    return ((struct range *) rbuf)[r].min;
}

int range_max(rbuf, r)
char *rbuf;
int r;
{
    return ((struct range *) rbuf)[r].max;
}

int range_mult(rbuf, r)
char *rbuf;
int r;
{
    return ((struct range *) rbuf)[r].mult;
}

void range_vals(rbuf, r, mn, mx, mult)
char *rbuf;
int r;
int *mn;
int *mx;
int *mult;
{
    *mn = ((struct range *) rbuf)[r].min;
    *mx = ((struct range *) rbuf)[r].max;
    *mult = ((struct range *) rbuf)[r].mult;
}


int *ranges_to_ivect(char *str, int minval, int maxval, int *numvals)
{
    int numranges = 0, numvalues = 0, ii, jj, kk = 0;
    int *values = NULL;
    char *ranges = NULL;

    numranges = parse_ranges(str, minval, maxval, 1, NULL, &ranges, NULL);
#ifdef DEBUGPRINT
    printf("\nThere are %d ranges...\n\n", numranges);
#endif
    if (numranges == 0) {
        *numvals = 0;
        return NULL;
    } else if (numranges > 0) {
        int mn, mx, mult;

        /* Count the number of values */
        for (ii = 0; ii < numranges; ii++) {
            int numinrange = 0;

            range_vals(ranges, ii, &mn, &mx, &mult);
            if (mn < minval)
                mn = minval;
            if (mn > maxval)
                continue;
            if (mx > maxval)
                mx = maxval;
            if (mx < minval)
                continue;
#ifdef DEBUGPRINT
            printf("%d-%d by %d's:  ", mn, mx, mult);
#endif
            for (jj = mn; jj <= mx; jj += mult) {
                numinrange++;
#ifdef DEBUGPRINT
                printf("%d ", jj);
#endif
            }
#ifdef DEBUGPRINT
            printf("\n");
#endif
            numvalues += numinrange;
        }
#ifdef DEBUGPRINT
        printf("\n %d values total.\n\n", numvalues);
#endif

        values = (int *) malloc((size_t) (sizeof(int) * numvalues));
        if (!values) {
            perror("\nAllocation error in ranges_to_ivect()");
            printf("\n");
            exit(-1);
        }

        /* Set the values in the array */
        for (ii = 0; ii < numranges; ii++) {
            range_vals(ranges, ii, &mn, &mx, &mult);
            if (mn < minval)
                mn = minval;
            if (mn > maxval)
                continue;
            if (mx > maxval)
                mx = maxval;
            if (mx < minval)
                continue;
            for (jj = mn; jj <= mx; jj += mult, kk++)
                values[kk] = jj;
        }
    }
    free(ranges);
    *numvals = numvalues;
    return values;
}


/*
int main(int argc, char *argv[])
{
  if (argc > 1){
    int ii, *values, numvalues;

    values = ranges_to_ivect(argv[1], 0, 2048, &numvalues);
    printf("\nThere are %d values...\n", numvalues);
    for (ii=0; ii<numvalues; ii++){
      printf(" %d", values[ii]);
    }
    printf("\n\n");
    free(values);
  } else {
    printf("\nusage:  range_parse  range_string\n");
  }
  exit(0);
}
*/
