/*
*  Name:
*     sla.c

*  Purpose:
*     Implement a C interface to the Fortran SLALIB library.

*  Description:
*     This file implements a C interface to the Fortran version of the
*     SLALIB library. 

*  Notes:
*     This interface only supports a subset of the functions provided by
*     SLALIB. It should be extended as and when necessary.

*  Copyright:
*     Copyright (C) 2006 Council for the Central Laboratory of the
*     Research Councils

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public Licence as
*     published by the Free Software Foundation; either version 2 of
*     the Licence, or (at your option) any later version.
*     
*     This program is distributed in the hope that it will be
*     useful,but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public Licence for more details.
*     
*     You should have received a copy of the GNU General Public Licence
*     along with this program; if not, write to the Free Software
*     Foundation, Inc., 59 Temple Place,Suite 330, Boston, MA
*     02111-1307, USA

*  Authors:
*     RFWS: R.F. Warren-Smith (STARLINK)
*     DSB: David S. Berry (STARLINK)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     SMR: Scott M. RAnsom (NRAO)

*  History:
*     12-NOV-1996 (RFWS):
*        Original version.
*     28-APR-1997 (RFWS):
*        Added SLA_DJCAL.
*     26-SEP-1997 (DSB):
*        Added SLA_DD2TF, SLA_DJCL.
*     21-JUN-2001 (DSB):
*        Added SLA_DBEAR, SLA_DVDV.
*     23-AUG-2001 (DSB):
*        Added SLA_SVD and SLA_SVDSOL
*     11-NOV-2002 (DSB):
*        Added SLA_RVEROT, SLA_GMST, SLA_EQEQX, SLA_RVLSRK, SLA_RVLSRD, 
*        SLA_RVLG, SLA_RVGALC.
*     11-JUN-2003 (DSB):
*        Added SLA_GEOC, SLA_HFK5Z and SLA_FK5HZ.
*     2-DEC-2004 (DSB):
*        Added SLA_DEULER.
*     29-SEP-2005 (DSB):
*        Added SLA_DE2H and SLA_DH2E
*     12-JUN-2006 (DSB):
*        Moved from AST to SLALIB.
*     25-JUN-2006 (TIMJ):
*        Add SLA_AIRMAS.
*     07-AUG-2006 (TIMJ):
*        Import cnfImprt from CNF.
*        Add SLA_OBS
*     08-AUG-2006 (TIMJ):
*        Add SLA_PA
*     12-DEC-2006 (TIMJ):
*        Add SLA_DTT and SLA_DAT
*     21-DEC-2006 (TIMJ):
*        Add SLA_RDPLAN
*     04-JUN-2007 (SMR):
*        Add SLA_OAP and SLA_AMP

*-
*/

/* Header files. */
/* ============= */
#include "f77.h"                /* FORTRAN <-> C interface macros (SUN/209) */
#include "slalib.h"             /* Prototypes for C SLALIB functions */
#include <stdlib.h>             /* Malloc etc */
#include <string.h>             /* string manipulation */


/* Functions needed to avoid a dependence on CNF. */
/* ============================================== */

static void slaStringExport(const char *, char *, int);
static void slaStringImport(const char *source_f, int source_len, char *dest_c);

static void slaStringExport(const char *source_c, char *dest_f, int dest_len)
{
/*
*+
*  Name:
*     slaStringExport

*  Purpose:
*     Export a C string to a FORTRAN string.

*  Type:
*     Protected function.

*  Synopsis:
*     void slaStringExport( const char *source_c, char *dest_f, int dest_len )

*  Description:
*     This function creates a FORTRAN string from a C string, storing
*     it in the supplied memory. If the C string is shorter than the
*     space allocated for the FORTRAN string, then it is padded with
*     blanks. If the C string is longer than the space allocated for
*     the FORTRAN string, then the string is truncated.

*  Parameters:
*     source_c 
*        A pointer to the input C string.
*     dest_f 
*        A pointer to the output FORTRAN string.
*     dest_len 
*        The length of the output FORTRAN string.

*  Notes:
*     - This function is potentially platform-specific. For example,
*     if FORTRAN strings were passed by descriptor, then the
*     descriptor address would be passed as "dest_f" and this must
*     then be used to locate the actual FORTRAN character data.
*     - This function is equivalent to cnfExprt but is included here to
*     avoid SLALIB becoming dependent on CNF.
*-
*/

/* Local Variables:*/
    int i;                      /* Loop counter for characters */

/* Check the supplied pointers. */
    if (!source_c || !dest_f)
        return;

/* Copy the characters of the input C string to the output FORTRAN
   string, taking care not to go beyond the end of the FORTRAN
   string.*/
    for (i = 0; source_c[i] && (i < dest_len); i++) {
        dest_f[i] = source_c[i];
    }

/* Fill the rest of the output FORTRAN string with blanks. */
    for (; i < dest_len; i++)
        dest_f[i] = ' ';
}

void slaStringImport(const char *source_f, int source_len, char *dest_c)

/*
*+
*  Name:
*     slaStringImportt

*  Purpose:
*     Import a FORTRAN string into a C string

*  Type:
*     Protected function.

*  Language:
*     ANSI C

*  Invocation:
*     slaStringImport( source_f, source_len, dest_c )

*  Description:
*     Import a FORTRAN string into a C string, discarding trailing
*     blanks. The NUL character is appended to the C string after
*     the last non-blank character. The input string and output string
*     pointers can point to the same location if the string is to be
*     modified in place (but care must be taken to allow for the additional
*     C terminator when allocating memory).

*  Arguments:
*     const char *source_f (Given)
*        A pointer to the input FORTRAN string
*     int source_len (Given)
*        The length of the input FORTRAN string
*     char *dest_c (Returned via pointer)
*        A pointer to the output C string. Can be same as source.

*  Notes:
*     -  No check is made that there is sufficient space allocated to
*        the C string to hold the FORTRAN string and a terminating null.
*        It is responsibility of the programmer to check this.
*     -  This function is equivalent to cnfImprt but is included here to
*        avoid SLALIB becoming dependent on CNF.

*  Authors:
*     PMA: Peter Allan (Starlink, RAL)
*     AJC: Alan Chipperfield (Starlink, RAL)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     27-MAR-1991 (PMA):
*        Original version.
*     22-MAY-1996 (AJC):
*        Correct description re trailing blanks
*     24-SEP-1998 (AJC):
*        Specify const char * for input strings
*     25-NOV-2005 (TIMJ):
*        Allow the strings to be identical
*     {enter_changes_here}

*-

*...........................................................................*/
{
/* Local Variables:							    */

    int i;                      /* Loop counter                            */


/* Find the last non blank character in the input FORTRAN string.	    */

    for (i = source_len - 1; (i >= 0) && (source_f[i] == ' '); i--);

/* Put a null character at the end of the output C string.		    */

    dest_c[i + 1] = '\0';

/* Copy the characters from the input FORTRAN string to the output C	    */
/* string if the strings are different.				      	    */

    if (dest_c != source_f) {
        memmove(dest_c, source_f, (size_t) i + 1);
    }
}




/* SLALIB wrapper implementations. */
/* =============================== */
/* Fortran routine prototype. */
F77_SUBROUTINE(sla_addet) (DOUBLE(RM),
                           DOUBLE(DM), DOUBLE(EQ), DOUBLE(RC), DOUBLE(DC));

/* C interface implementation. */
void slaAddet(double rm, double dm, double eq, double *rc, double *dc)
{
    DECLARE_DOUBLE(RM);
    DECLARE_DOUBLE(DM);
    DECLARE_DOUBLE(EQ);
    DECLARE_DOUBLE(RC);
    DECLARE_DOUBLE(DC);
    RM = rm;
    DM = dm;
    EQ = eq;
    F77_CALL(sla_addet) (DOUBLE_ARG(&RM),
                         DOUBLE_ARG(&DM),
                         DOUBLE_ARG(&EQ), DOUBLE_ARG(&RC), DOUBLE_ARG(&DC));
    *rc = RC;
    *dc = DC;
}

/* etc... */
F77_SUBROUTINE(sla_amp) (DOUBLE(RA),
                         DOUBLE(DA),
                         DOUBLE(DATE), DOUBLE(EQ), DOUBLE(RM), DOUBLE(DM));

void slaAmp(double ra, double da, double date, double eq, double *rm, double *dm)
{
    DECLARE_DOUBLE(RA);
    DECLARE_DOUBLE(DA);
    DECLARE_DOUBLE(DATE);
    DECLARE_DOUBLE(EQ);
    DECLARE_DOUBLE(RM);
    DECLARE_DOUBLE(DM);
    RA = ra;
    DA = da;
    DATE = date;
    EQ = eq;
    F77_CALL(sla_amp) (DOUBLE_ARG(&RA),
                       DOUBLE_ARG(&DA),
                       DOUBLE_ARG(&DATE),
                       DOUBLE_ARG(&EQ), DOUBLE_ARG(&RM), DOUBLE_ARG(&DM));
    *rm = RM;
    *dm = DM;
}

F77_SUBROUTINE(sla_ampqk) (DOUBLE(RA),
                           DOUBLE(DA), DOUBLE_ARRAY(AMPRMS), DOUBLE(RM), DOUBLE(DM));

void slaAmpqk(double ra, double da, double amprms[21], double *rm, double *dm)
{
    DECLARE_DOUBLE(RA);
    DECLARE_DOUBLE(DA);
    DECLARE_DOUBLE_ARRAY(AMPRMS, 21);
    DECLARE_DOUBLE(RM);
    DECLARE_DOUBLE(DM);
    int i;
    RA = ra;
    DA = da;
    for (i = 0; i < 21; i++)
        AMPRMS[i] = amprms[i];
    F77_CALL(sla_ampqk) (DOUBLE_ARG(&RA),
                         DOUBLE_ARG(&DA),
                         DOUBLE_ARRAY_ARG(AMPRMS), DOUBLE_ARG(&RM), DOUBLE_ARG(&DM));
    *rm = RM;
    *dm = DM;
}

F77_DOUBLE_FUNCTION(sla_airmas) (DOUBLE(ZD));

double slaAirmas(double zd)
{
    DECLARE_DOUBLE(ZD);
    double result;
    ZD = zd;
    result = F77_CALL(sla_airmas) (DOUBLE_ARG(&ZD));
    return result;
}

F77_SUBROUTINE(sla_caldj) (INTEGER(IY),
                           INTEGER(IM), INTEGER(ID), DOUBLE(DJM), INTEGER(J));

void slaCaldj(int iy, int im, int id, double *djm, int *j)
{
    DECLARE_INTEGER(IY);
    DECLARE_INTEGER(IM);
    DECLARE_INTEGER(ID);
    DECLARE_DOUBLE(DJM);
    DECLARE_INTEGER(J);
    IY = iy;
    IM = im;
    ID = id;
    F77_CALL(sla_caldj) (INTEGER_ARG(&IY),
                         INTEGER_ARG(&IM),
                         INTEGER_ARG(&ID), DOUBLE_ARG(&DJM), INTEGER_ARG(&J));
    *djm = DJM;
    *j = J;
}

F77_SUBROUTINE(sla_daf2r) (INTEGER(IDEG),
                           INTEGER(IAMIN), DOUBLE(ASEC), DOUBLE(RAD), INTEGER(J));

void slaDaf2r(int ideg, int iamin, double asec, double *rad, int *j)
{
    DECLARE_INTEGER(IDEG);
    DECLARE_INTEGER(IAMIN);
    DECLARE_DOUBLE(ASEC);
    DECLARE_DOUBLE(RAD);
    DECLARE_INTEGER(J);
    IDEG = ideg;
    IAMIN = iamin;
    ASEC = asec;
    F77_CALL(sla_daf2r) (INTEGER_ARG(&IDEG),
                         INTEGER_ARG(&IAMIN),
                         DOUBLE_ARG(&ASEC), DOUBLE_ARG(&RAD), INTEGER_ARG(&J));
    *rad = RAD;
    *j = J;
}

F77_SUBROUTINE(sla_dav2m) (DOUBLE_ARRAY(AXVEC), DOUBLE_ARRAY(RMAT));

void slaDav2m(double axvec[3], double rmat[3][3])
{
    DECLARE_DOUBLE_ARRAY(AXVEC, 3);
    DECLARE_DOUBLE_ARRAY(RMAT, 9);
    int i;
    int j;
    for (i = 0; i < 3; i++)
        AXVEC[i] = axvec[i];
    F77_CALL(sla_dav2m) (DOUBLE_ARRAY_ARG(AXVEC), DOUBLE_ARRAY_ARG(RMAT));
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            rmat[i][j] = RMAT[i + 3 * j];
    }
}

F77_SUBROUTINE(sla_dcc2s) (DOUBLE_ARRAY(V), DOUBLE(A), DOUBLE(B));

void slaDcc2s(double v[3], double *a, double *b)
{
    DECLARE_DOUBLE_ARRAY(V, 3);
    DECLARE_DOUBLE(A);
    DECLARE_DOUBLE(B);
    int i;
    for (i = 0; i < 3; i++)
        V[i] = v[i];
    F77_CALL(sla_dcc2s) (DOUBLE_ARRAY_ARG(V), DOUBLE_ARG(&A), DOUBLE_ARG(&B));
    *a = A;
    *b = B;
}

F77_SUBROUTINE(sla_dcs2c) (DOUBLE(A), DOUBLE(B), DOUBLE_ARRAY(V));

void slaDcs2c(double a, double b, double v[3])
{
    DECLARE_DOUBLE(A);
    DECLARE_DOUBLE(B);
    DECLARE_DOUBLE_ARRAY(V, 3);
    int i;
    A = a;
    B = b;
    F77_CALL(sla_dcs2c) (DOUBLE_ARG(&A), DOUBLE_ARG(&B), DOUBLE_ARRAY_ARG(V));
    for (i = 0; i < 3; i++)
        v[i] = V[i];
}

F77_SUBROUTINE(sla_dd2tf) (INTEGER(NDP),
                           DOUBLE(DAYS), CHARACTER(SIGN), INTEGER_ARRAY(IHMSF)
                           TRAIL(SIGN));

void slaDd2tf(int ndp, double days, char *sign, int ihmsf[4])
{
    DECLARE_INTEGER(NDP);
    DECLARE_DOUBLE(DAYS);
    DECLARE_CHARACTER(SIGN, 2);
    DECLARE_INTEGER_ARRAY(IHMSF, 4);
    int i;

    NDP = ndp;
    DAYS = days;
    F77_CALL(sla_dd2tf) (INTEGER_ARG(&NDP),
                         DOUBLE_ARG(&DAYS),
                         CHARACTER_ARG(SIGN), INTEGER_ARRAY_ARG(IHMSF)
                         TRAIL_ARG(SIGN));
    sign[0] = SIGN[0];
    sign[1] = 0;
    for (i = 0; i < 4; i++)
        ihmsf[i] = IHMSF[i];
}

F77_SUBROUTINE(sla_dimxv) (DOUBLE_ARRAY(DM), DOUBLE_ARRAY(VA), DOUBLE_ARRAY(VB));

void slaDimxv(double dm[3][3], double va[3], double vb[3])
{
    DECLARE_DOUBLE_ARRAY(DM, 9);
    DECLARE_DOUBLE_ARRAY(VA, 3);
    DECLARE_DOUBLE_ARRAY(VB, 3);
    int i;
    int j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            DM[i + j * 3] = dm[i][j];
        VA[i] = va[i];
    }
    F77_CALL(sla_dimxv) (DOUBLE_ARRAY_ARG(DM),
                         DOUBLE_ARRAY_ARG(VA), DOUBLE_ARRAY_ARG(VB));
    for (i = 0; i < 3; i++)
        vb[i] = VB[i];
}

F77_SUBROUTINE(sla_djcal) (INTEGER(NDP),
                           DOUBLE(DJM), INTEGER_ARRAY(IYMDF), INTEGER(J));

void slaDjcal(int ndp, double djm, int iymdf[4], int *j)
{
    DECLARE_INTEGER(NDP);
    DECLARE_DOUBLE(DJM);
    DECLARE_INTEGER_ARRAY(IYMDF, 4);
    DECLARE_INTEGER(J);
    int i;

    NDP = ndp;
    DJM = djm;
    F77_CALL(sla_djcal) (INTEGER_ARG(&NDP),
                         DOUBLE_ARG(&DJM),
                         INTEGER_ARRAY_ARG(IYMDF), INTEGER_ARG(&J));
    for (i = 0; i < 4; i++)
        iymdf[i] = IYMDF[i];
    *j = J;
}

F77_SUBROUTINE(sla_djcl) (DOUBLE(DJM),
                          INTEGER(IY),
                          INTEGER(IM), INTEGER(ID), DOUBLE(FD), INTEGER(J));

void slaDjcl(double djm, int *iy, int *im, int *id, double *fd, int *j)
{
    DECLARE_DOUBLE(DJM);
    DECLARE_INTEGER(IY);
    DECLARE_INTEGER(IM);
    DECLARE_INTEGER(ID);
    DECLARE_DOUBLE(FD);
    DECLARE_INTEGER(J);

    DJM = djm;
    F77_CALL(sla_djcl) (DOUBLE_ARG(&DJM),
                        INTEGER_ARG(&IY),
                        INTEGER_ARG(&IM),
                        INTEGER_ARG(&ID), DOUBLE_ARG(&FD), INTEGER_ARG(&J));
    *iy = IY;
    *im = IM;
    *id = ID;
    *fd = FD;
    *j = J;
}

F77_SUBROUTINE(sla_dmat) (INTEGER(N),
                          DOUBLE_ARRAY(A),
                          DOUBLE_ARRAY(Y),
                          DOUBLE(D), INTEGER(JF), INTEGER_ARRAY(IW));

void slaDmat(int n, double *a, double *y, double *d, int *jf, int *iw)
{
    DECLARE_INTEGER(N);
    F77_DOUBLE_TYPE *A;
    F77_DOUBLE_TYPE *Y;
    DECLARE_DOUBLE(D);
    DECLARE_INTEGER(JF);
    F77_INTEGER_TYPE *IW;
    int i;
    int j;
    A = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) (n * n));
    Y = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) n);
    if (sizeof(F77_INTEGER_TYPE) > sizeof(int)) {
        IW = malloc(sizeof(F77_INTEGER_TYPE) * (size_t) n);
    } else {
        IW = (F77_INTEGER_TYPE *) iw;
    }
    if (IW) {
        N = n;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                A[i + n * j] = a[n * i + j];
            Y[i] = y[i];
        }
        F77_CALL(sla_dmat) (INTEGER_ARG(&N), DOUBLE_ARRAY_ARG(A),
                            DOUBLE_ARRAY_ARG(Y), DOUBLE_ARG(&D),
                            INTEGER_ARG(&JF), INTEGER_ARG(IW));
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                a[n * i + j] = A[i + n * j];
            y[i] = Y[i];
        }
        *d = D;
        *jf = JF;
    }
    free(A);
    free(Y);
    if (sizeof(F77_INTEGER_TYPE) > sizeof(int))
        free(IW);
}

F77_SUBROUTINE(sla_dmxm) (DOUBLE_ARRAY(A), DOUBLE_ARRAY(B), DOUBLE_ARRAY(C));

void slaDmxm(double a[3][3], double b[3][3], double c[3][3])
{
    DECLARE_DOUBLE_ARRAY(A, 9);
    DECLARE_DOUBLE_ARRAY(B, 9);
    DECLARE_DOUBLE_ARRAY(C, 9);
    int i;
    int j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            A[i + 3 * j] = a[i][j];
            B[i + 3 * j] = b[i][j];
        }
    }
    F77_CALL(sla_dmxm) (DOUBLE_ARRAY_ARG(A),
                        DOUBLE_ARRAY_ARG(B), DOUBLE_ARRAY_ARG(C));
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            c[i][j] = C[i + 3 * j];
    }
}

F77_SUBROUTINE(sla_dmxv) (DOUBLE_ARRAY(DM), DOUBLE_ARRAY(VA), DOUBLE_ARRAY(VB));

void slaDmxv(double dm[3][3], double va[3], double vb[3])
{
    DECLARE_DOUBLE_ARRAY(DM, 9);
    DECLARE_DOUBLE_ARRAY(VA, 3);
    DECLARE_DOUBLE_ARRAY(VB, 3);
    int i;
    int j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            DM[i + 3 * j] = dm[i][j];
        VA[i] = va[i];
    }
    F77_CALL(sla_dmxv) (DOUBLE_ARRAY_ARG(DM),
                        DOUBLE_ARRAY_ARG(VA), DOUBLE_ARRAY_ARG(VB));
    for (i = 0; i < 3; i++)
        vb[i] = VB[i];
}

F77_DOUBLE_FUNCTION(sla_dbear) (DOUBLE(A1), DOUBLE(B1), DOUBLE(A2), DOUBLE(B2));

double slaDbear(double a1, double b1, double a2, double b2)
{
    DECLARE_DOUBLE(A1);
    DECLARE_DOUBLE(B1);
    DECLARE_DOUBLE(A2);
    DECLARE_DOUBLE(B2);
    double result;
    A1 = a1;
    B1 = b1;
    A2 = a2;
    B2 = b2;
    result = F77_CALL(sla_dbear) (DOUBLE_ARG(&A1), DOUBLE_ARG(&B1),
                                  DOUBLE_ARG(&A2), DOUBLE_ARG(&B2));
    return result;
}

F77_DOUBLE_FUNCTION(sla_drange) (DOUBLE(ANGLE));

double slaDrange(double angle)
{
    DECLARE_DOUBLE(ANGLE);
    double result;
    ANGLE = angle;
    result = F77_CALL(sla_drange) (DOUBLE_ARG(&ANGLE));
    return result;
}

F77_DOUBLE_FUNCTION(sla_dranrm) (DOUBLE(ANGLE));

double slaDranrm(double angle)
{
    DECLARE_DOUBLE(ANGLE);
    double result;
    ANGLE = angle;
    result = F77_CALL(sla_dranrm) (DOUBLE_ARG(&ANGLE));
    return result;
}

F77_DOUBLE_FUNCTION(sla_dsep) (DOUBLE(A1), DOUBLE(B1), DOUBLE(A2), DOUBLE(B2));

double slaDsep(double a1, double b1, double a2, double b2)
{
    DECLARE_DOUBLE(A1);
    DECLARE_DOUBLE(B1);
    DECLARE_DOUBLE(A2);
    DECLARE_DOUBLE(B2);
    double result;
    A1 = a1;
    B1 = b1;
    A2 = a2;
    B2 = b2;
    result = F77_CALL(sla_dsep) (DOUBLE_ARG(&A1),
                                 DOUBLE_ARG(&B1), DOUBLE_ARG(&A2), DOUBLE_ARG(&B2));
    return result;
}

F77_DOUBLE_FUNCTION(sla_dvdv) (DOUBLE_ARRAY(VA), DOUBLE_ARRAY(VB));

double slaDvdv(double va[3], double vb[3])
{
    DECLARE_DOUBLE_ARRAY(VA, 3);
    DECLARE_DOUBLE_ARRAY(VB, 3);
    double result;
    int i;
    for (i = 0; i < 3; i++) {
        VA[i] = va[i];
        VB[i] = vb[i];
    }
    result = F77_CALL(sla_dvdv) (DOUBLE_ARRAY_ARG(VA), DOUBLE_ARRAY_ARG(VB));
    return result;
}

F77_SUBROUTINE(sla_dtf2d) (INTEGER(IHOUR),
                           INTEGER(IMIN), DOUBLE(SEC), DOUBLE(DAYS), INTEGER(J));

void slaDtf2d(int ihour, int imin, double sec, double *days, int *j)
{
    DECLARE_INTEGER(IHOUR);
    DECLARE_INTEGER(IMIN);
    DECLARE_DOUBLE(SEC);
    DECLARE_DOUBLE(DAYS);
    DECLARE_INTEGER(J);
    IHOUR = ihour;
    IMIN = imin;
    SEC = sec;
    F77_CALL(sla_dtf2d) (INTEGER_ARG(&IHOUR),
                         INTEGER_ARG(&IMIN),
                         DOUBLE_ARG(&SEC), DOUBLE_ARG(&DAYS), INTEGER_ARG(&J));
    *days = DAYS;
    *j = J;
}

F77_SUBROUTINE(sla_dtf2r) (INTEGER(IHOUR),
                           INTEGER(IMIN), DOUBLE(SEC), DOUBLE(RAD), INTEGER(J));

void slaDtf2r(int ihour, int imin, double sec, double *rad, int *j)
{
    DECLARE_INTEGER(IHOUR);
    DECLARE_INTEGER(IMIN);
    DECLARE_DOUBLE(SEC);
    DECLARE_DOUBLE(RAD);
    DECLARE_INTEGER(J);
    IHOUR = ihour;
    IMIN = imin;
    SEC = sec;
    F77_CALL(sla_dtf2r) (INTEGER_ARG(&IHOUR),
                         INTEGER_ARG(&IMIN),
                         DOUBLE_ARG(&SEC), DOUBLE_ARG(&RAD), INTEGER_ARG(&J));
    *rad = RAD;
    *j = J;
}

F77_DOUBLE_FUNCTION(sla_dt) (DOUBLE(EPOCH));

double slaDt(double epoch)
{
    DECLARE_DOUBLE(EPOCH);
    double result;
    EPOCH = epoch;
    result = F77_CALL(sla_dt) (DOUBLE_ARG(&EPOCH));
    return result;
}

F77_SUBROUTINE(sla_dvn) (DOUBLE_ARRAY(V), DOUBLE_ARRAY(UV), DOUBLE(VM));

void slaDvn(double v[3], double uv[3], double *vm)
{
    DECLARE_DOUBLE_ARRAY(V, 3);
    DECLARE_DOUBLE_ARRAY(UV, 3);
    DECLARE_DOUBLE(VM);
    int i;
    for (i = 0; i < 3; i++)
        V[i] = v[i];
    F77_CALL(sla_dvn) (DOUBLE_ARRAY_ARG(V), DOUBLE_ARRAY_ARG(UV), DOUBLE_ARG(&VM));
    for (i = 0; i < 3; i++)
        uv[i] = UV[i];
    *vm = VM;
}

F77_SUBROUTINE(sla_dvxv) (DOUBLE_ARRAY(VA), DOUBLE_ARRAY(VB), DOUBLE_ARRAY(VC));

void slaDvxv(double va[3], double vb[3], double vc[3])
{
    DECLARE_DOUBLE_ARRAY(VA, 3);
    DECLARE_DOUBLE_ARRAY(VB, 3);
    DECLARE_DOUBLE_ARRAY(VC, 3);
    int i;
    for (i = 0; i < 3; i++) {
        VA[i] = va[i];
        VB[i] = vb[i];
    }
    F77_CALL(sla_dvxv) (DOUBLE_ARRAY_ARG(VA),
                        DOUBLE_ARRAY_ARG(VB), DOUBLE_ARRAY_ARG(VC));
    for (i = 0; i < 3; i++)
        vc[i] = VC[i];
}

F77_SUBROUTINE(sla_ecmat) (DOUBLE(DATE), DOUBLE_ARRAY(RMAT));

void slaEcmat(double date, double rmat[3][3])
{
    DECLARE_DOUBLE(DATE);
    DECLARE_DOUBLE_ARRAY(RMAT, 9);
    int i;
    int j;
    DATE = date;
    F77_CALL(sla_ecmat) (DOUBLE_ARG(&DATE), DOUBLE_ARRAY_ARG(RMAT));
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            rmat[i][j] = RMAT[i + 3 * j];
    }
}

F77_DOUBLE_FUNCTION(sla_epb) (DOUBLE(DATE));

double slaEpb(double date)
{
    DECLARE_DOUBLE(DATE);
    double result;
    DATE = date;
    result = F77_CALL(sla_epb) (DOUBLE_ARG(&DATE));
    return result;
}

F77_DOUBLE_FUNCTION(sla_epb2d) (DOUBLE(EPB));

double slaEpb2d(double epb)
{
    DECLARE_DOUBLE(EPB);
    double result;
    EPB = epb;
    result = F77_CALL(sla_epb2d) (DOUBLE_ARG(&EPB));
    return result;
}

F77_DOUBLE_FUNCTION(sla_epj) (DOUBLE(DATE));

double slaEpj(double date)
{
    DECLARE_DOUBLE(DATE);
    double result;
    DATE = date;
    result = F77_CALL(sla_epj) (DOUBLE_ARG(&DATE));
    return result;
}

F77_DOUBLE_FUNCTION(sla_epj2d) (DOUBLE(EPJ));

double slaEpj2d(double epj)
{
    DECLARE_DOUBLE(EPJ);
    double result;
    EPJ = epj;
    result = F77_CALL(sla_epj2d) (DOUBLE_ARG(&EPJ));
    return result;
}

F77_DOUBLE_FUNCTION(sla_eqeqx) (DOUBLE(DATE));

double slaEqeqx(double date)
{
    DECLARE_DOUBLE(DATE);
    double result;
    DATE = date;
    result = F77_CALL(sla_eqeqx) (DOUBLE_ARG(&DATE));
    return result;
}

F77_SUBROUTINE(sla_eqgal) (DOUBLE(DR), DOUBLE(DD), DOUBLE(DL), DOUBLE(DB));

void slaEqgal(double dr, double dd, double *dl, double *db)
{
    DECLARE_DOUBLE(DR);
    DECLARE_DOUBLE(DD);
    DECLARE_DOUBLE(DL);
    DECLARE_DOUBLE(DB);
    DR = dr;
    DD = dd;
    F77_CALL(sla_eqgal) (DOUBLE_ARG(&DR),
                         DOUBLE_ARG(&DD), DOUBLE_ARG(&DL), DOUBLE_ARG(&DB));
    *dl = DL;
    *db = DB;
}

F77_SUBROUTINE(sla_fk45z) (DOUBLE(R1950),
                           DOUBLE(D1950),
                           DOUBLE(BEPOCH), DOUBLE(R2000), DOUBLE(D2000));

void slaFk45z(double r1950, double d1950, double bepoch,
              double *r2000, double *d2000)
{
    DECLARE_DOUBLE(R1950);
    DECLARE_DOUBLE(D1950);
    DECLARE_DOUBLE(BEPOCH);
    DECLARE_DOUBLE(R2000);
    DECLARE_DOUBLE(D2000);
    R1950 = r1950;
    D1950 = d1950;
    BEPOCH = bepoch;
    F77_CALL(sla_fk45z) (DOUBLE_ARG(&R1950),
                         DOUBLE_ARG(&D1950),
                         DOUBLE_ARG(&BEPOCH),
                         DOUBLE_ARG(&R2000), DOUBLE_ARG(&D2000));
    *r2000 = R2000;
    *d2000 = D2000;
}

F77_SUBROUTINE(sla_fk54z) (DOUBLE(R2000),
                           DOUBLE(D2000),
                           DOUBLE(BEPOCH),
                           DOUBLE(R1950),
                           DOUBLE(D1950), DOUBLE(DR1950), DOUBLE(DD1950));

void slaFk54z(double r2000, double d2000, double bepoch,
              double *r1950, double *d1950, double *dr1950, double *dd1950)
{
    DECLARE_DOUBLE(R2000);
    DECLARE_DOUBLE(D2000);
    DECLARE_DOUBLE(BEPOCH);
    DECLARE_DOUBLE(R1950);
    DECLARE_DOUBLE(D1950);
    DECLARE_DOUBLE(DR1950);
    DECLARE_DOUBLE(DD1950);
    R2000 = r2000;
    D2000 = d2000;
    BEPOCH = bepoch;
    F77_CALL(sla_fk54z) (DOUBLE_ARG(&R2000),
                         DOUBLE_ARG(&D2000),
                         DOUBLE_ARG(&BEPOCH),
                         DOUBLE_ARG(&R1950),
                         DOUBLE_ARG(&D1950),
                         DOUBLE_ARG(&DR1950), DOUBLE_ARG(&DD1950));
    *r1950 = R1950;
    *d1950 = D1950;
    *dr1950 = DR1950;
    *dd1950 = DD1950;
}

F77_SUBROUTINE(sla_galeq) (DOUBLE(DL), DOUBLE(DB), DOUBLE(DR), DOUBLE(DD));

void slaGaleq(double dl, double db, double *dr, double *dd)
{
    DECLARE_DOUBLE(DL);
    DECLARE_DOUBLE(DB);
    DECLARE_DOUBLE(DR);
    DECLARE_DOUBLE(DD);
    DL = dl;
    DB = db;
    F77_CALL(sla_galeq) (DOUBLE_ARG(&DL),
                         DOUBLE_ARG(&DB), DOUBLE_ARG(&DR), DOUBLE_ARG(&DD));
    *dr = DR;
    *dd = DD;
}

F77_SUBROUTINE(sla_galsup) (DOUBLE(DL), DOUBLE(DB), DOUBLE(DSL), DOUBLE(DSB));

void slaGalsup(double dl, double db, double *dsl, double *dsb)
{
    DECLARE_DOUBLE(DL);
    DECLARE_DOUBLE(DB);
    DECLARE_DOUBLE(DSL);
    DECLARE_DOUBLE(DSB);
    DL = dl;
    DB = db;
    F77_CALL(sla_galsup) (DOUBLE_ARG(&DL),
                          DOUBLE_ARG(&DB), DOUBLE_ARG(&DSL), DOUBLE_ARG(&DSB));
    *dsl = DSL;
    *dsb = DSB;
}

F77_DOUBLE_FUNCTION(sla_gmst) (DOUBLE(UT1));

double slaGmst(double ut1)
{
    DECLARE_DOUBLE(UT1);
    double result;
    UT1 = ut1;
    result = F77_CALL(sla_gmst) (DOUBLE_ARG(&UT1));
    return result;
}

F77_SUBROUTINE(sla_mappa) (DOUBLE(EQ), DOUBLE(DATE), DOUBLE_ARRAY(AMPRMS));

void slaMappa(double eq, double date, double amprms[21])
{
    DECLARE_DOUBLE(EQ);
    DECLARE_DOUBLE(DATE);
    DECLARE_DOUBLE_ARRAY(AMPRMS, 21);
    int i;
    EQ = eq;
    DATE = date;
    F77_CALL(sla_mappa) (DOUBLE_ARG(&EQ),
                         DOUBLE_ARG(&DATE), DOUBLE_ARRAY_ARG(AMPRMS));
    for (i = 0; i < 21; i++)
        amprms[i] = AMPRMS[i];
}

F77_SUBROUTINE(sla_mapqkz) (DOUBLE(RM),
                            DOUBLE(DM),
                            DOUBLE_ARRAY(AMPRMS), DOUBLE(RA), DOUBLE(DA));

void slaMapqkz(double rm, double dm, double amprms[21], double *ra, double *da)
{
    DECLARE_DOUBLE(RM);
    DECLARE_DOUBLE(DM);
    DECLARE_DOUBLE_ARRAY(AMPRMS, 21);
    DECLARE_DOUBLE(RA);
    DECLARE_DOUBLE(DA);
    int i;
    RM = rm;
    DM = dm;
    for (i = 0; i < 21; i++)
        AMPRMS[i] = amprms[i];
    F77_CALL(sla_mapqkz) (DOUBLE_ARG(&RM),
                          DOUBLE_ARG(&DM),
                          DOUBLE_ARRAY_ARG(AMPRMS),
                          DOUBLE_ARG(&RA), DOUBLE_ARG(&DA));
    *ra = RA;
    *da = DA;
}

F77_SUBROUTINE(sla_prebn) (DOUBLE(BEP0), DOUBLE(BEP1), DOUBLE_ARRAY(RMATP));

void slaPrebn(double bep0, double bep1, double rmatp[3][3])
{
    DECLARE_DOUBLE(BEP0);
    DECLARE_DOUBLE(BEP1);
    DECLARE_DOUBLE_ARRAY(RMATP, 9);
    int i;
    int j;
    BEP0 = bep0;
    BEP1 = bep1;
    F77_CALL(sla_prebn) (DOUBLE_ARG(&BEP0),
                         DOUBLE_ARG(&BEP1), DOUBLE_ARRAY_ARG(RMATP));
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            rmatp[i][j] = RMATP[i + 3 * j];
    }
}

F77_SUBROUTINE(sla_prec) (DOUBLE(EP0), DOUBLE(EP1), DOUBLE_ARRAY(RMATP));

void slaPrec(double ep0, double ep1, double rmatp[3][3])
{
    DECLARE_DOUBLE(EP0);
    DECLARE_DOUBLE(EP1);
    DECLARE_DOUBLE_ARRAY(RMATP, 9);
    int i;
    int j;
    EP0 = ep0;
    EP1 = ep1;
    F77_CALL(sla_prec) (DOUBLE_ARG(&EP0), DOUBLE_ARG(&EP1), DOUBLE_ARRAY_ARG(RMATP));
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            rmatp[i][j] = RMATP[i + 3 * j];
    }
}

F77_REAL_FUNCTION(sla_rverot) (REAL(PHI), REAL(RA), REAL(DEC), REAL(ST));

float slaRverot(float phi, float ra, float dec, float st)
{
    DECLARE_REAL(PHI);
    DECLARE_REAL(RA);
    DECLARE_REAL(DEC);
    DECLARE_REAL(ST);
    float result;
    PHI = phi;
    RA = ra;
    DEC = dec;
    ST = st;
    result = F77_CALL(sla_rverot) (REAL_ARG(&PHI),
                                   REAL_ARG(&RA), REAL_ARG(&DEC), REAL_ARG(&ST));
    return result;
}

F77_REAL_FUNCTION(sla_rvgalc) (REAL(RA), REAL(DEC));

float slaRvgalc(float ra, float dec)
{
    DECLARE_REAL(RA);
    DECLARE_REAL(DEC);
    float result;
    RA = ra;
    DEC = dec;
    result = F77_CALL(sla_rvgalc) (REAL_ARG(&RA), REAL_ARG(&DEC));
    return result;
}

F77_REAL_FUNCTION(sla_rvlg) (REAL(RA), REAL(DEC));

float slaRvlg(float ra, float dec)
{
    DECLARE_REAL(RA);
    DECLARE_REAL(DEC);
    float result;
    RA = ra;
    DEC = dec;
    result = F77_CALL(sla_rvlg) (REAL_ARG(&RA), REAL_ARG(&DEC));
    return result;
}

F77_REAL_FUNCTION(sla_rvlsrd) (REAL(RA), REAL(DEC));

float slaRvlsrd(float ra, float dec)
{
    DECLARE_REAL(RA);
    DECLARE_REAL(DEC);
    float result;
    RA = ra;
    DEC = dec;
    result = F77_CALL(sla_rvlsrd) (REAL_ARG(&RA), REAL_ARG(&DEC));
    return result;
}

F77_REAL_FUNCTION(sla_rvlsrk) (REAL(RA), REAL(DEC));

float slaRvlsrk(float ra, float dec)
{
    DECLARE_REAL(RA);
    DECLARE_REAL(DEC);
    float result;
    RA = ra;
    DEC = dec;
    result = F77_CALL(sla_rvlsrk) (REAL_ARG(&RA), REAL_ARG(&DEC));
    return result;
}


F77_SUBROUTINE(sla_subet) (DOUBLE(RC),
                           DOUBLE(DC), DOUBLE(EQ), DOUBLE(RM), DOUBLE(DM));

void slaSubet(double rc, double dc, double eq, double *rm, double *dm)
{
    DECLARE_DOUBLE(RC);
    DECLARE_DOUBLE(DC);
    DECLARE_DOUBLE(EQ);
    DECLARE_DOUBLE(RM);
    DECLARE_DOUBLE(DM);
    RC = rc;
    DC = dc;
    EQ = eq;
    F77_CALL(sla_subet) (DOUBLE_ARG(&RC),
                         DOUBLE_ARG(&DC),
                         DOUBLE_ARG(&EQ), DOUBLE_ARG(&RM), DOUBLE_ARG(&DM));
    *rm = RM;
    *dm = DM;
}

F77_SUBROUTINE(sla_supgal) (DOUBLE(DSL), DOUBLE(DSB), DOUBLE(DL), DOUBLE(DB));

void slaSupgal(double dsl, double dsb, double *dl, double *db)
{
    DECLARE_DOUBLE(DSL);
    DECLARE_DOUBLE(DSB);
    DECLARE_DOUBLE(DL);
    DECLARE_DOUBLE(DB);
    DSL = dsl;
    DSB = dsb;
    F77_CALL(sla_supgal) (DOUBLE_ARG(&DSL),
                          DOUBLE_ARG(&DSB), DOUBLE_ARG(&DL), DOUBLE_ARG(&DB));
    *dl = DL;
    *db = DB;
}



F77_SUBROUTINE(sla_svd) (INTEGER(M),
                         INTEGER(N),
                         INTEGER(MP),
                         INTEGER(NP),
                         DOUBLE_ARRAY(A),
                         DOUBLE_ARRAY(W),
                         DOUBLE_ARRAY(V), DOUBLE_ARRAY(WORK), INTEGER(JSTAT));

void slaSvd(int m, int n, int mp, int np,
            double *a, double *w, double *v, double *work, int *jstat)
{
    DECLARE_INTEGER(M);
    DECLARE_INTEGER(N);
    DECLARE_INTEGER(MP);
    DECLARE_INTEGER(NP);
    F77_DOUBLE_TYPE *A;
    F77_DOUBLE_TYPE *W;
    F77_DOUBLE_TYPE *V;
    F77_DOUBLE_TYPE *WORK;
    DECLARE_INTEGER(JSTAT);


    int i;
    int j;

    A = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) (mp * np));
    W = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) n);
    V = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) (np * np));
    WORK = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) n);

    if (WORK) {
        M = m;
        N = n;
        MP = mp;
        NP = np;

        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++)
                A[i + mp * j] = a[np * i + j];
        }

        F77_CALL(sla_svd) (INTEGER_ARG(&M),
                           INTEGER_ARG(&N),
                           INTEGER_ARG(&MP),
                           INTEGER_ARG(&NP),
                           DOUBLE_ARRAY_ARG(A),
                           DOUBLE_ARRAY_ARG(W),
                           DOUBLE_ARRAY_ARG(V),
                           DOUBLE_ARRAY_ARG(WORK), INTEGER_ARG(&JSTAT));


        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++)
                a[np * i + j] = A[i + mp * j];
        }

        for (i = 0; i < n; i++) {
            w[i] = W[i];
            work[i] = WORK[i];
            for (j = 0; j < n; j++)
                v[np * i + j] = V[i + np * j];
        }

        *jstat = JSTAT;
    }

    free(A);
    free(W);
    free(V);
    free(WORK);
}

F77_SUBROUTINE(sla_svdsol) (INTEGER(M),
                            INTEGER(N),
                            INTEGER(MP),
                            INTEGER(NP),
                            DOUBLE_ARRAY(B),
                            DOUBLE_ARRAY(U),
                            DOUBLE_ARRAY(W),
                            DOUBLE_ARRAY(V), DOUBLE_ARRAY(WORK), DOUBLE_ARRAY(X));

void slaSvdsol(int m, int n, int mp, int np,
               double *b, double *u, double *w, double *v, double *work, double *x)
{

    DECLARE_INTEGER(M);
    DECLARE_INTEGER(N);
    DECLARE_INTEGER(MP);
    DECLARE_INTEGER(NP);
    F77_DOUBLE_TYPE *B;
    F77_DOUBLE_TYPE *U;
    F77_DOUBLE_TYPE *W;
    F77_DOUBLE_TYPE *V;
    F77_DOUBLE_TYPE *WORK;
    F77_DOUBLE_TYPE *X;

    int i;
    int j;

    B = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) (m));
    U = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) (mp * np));
    W = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) n);
    V = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) (np * np));
    WORK = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) n);
    X = malloc(sizeof(F77_DOUBLE_TYPE) * (size_t) n);

    if (X) {
        M = m;
        N = n;
        MP = mp;
        NP = np;

        for (i = 0; i < m; i++) {
            B[i] = b[i];
            for (j = 0; j < n; j++)
                U[i + mp * j] = u[np * i + j];
        }
        for (i = 0; i < n; i++) {
            W[i] = w[i];
            for (j = 0; j < n; j++)
                V[i + np * j] = v[np * i + j];
        }

        F77_CALL(sla_svdsol) (INTEGER_ARG(&M),
                              INTEGER_ARG(&N),
                              INTEGER_ARG(&MP),
                              INTEGER_ARG(&NP),
                              DOUBLE_ARRAY_ARG(B),
                              DOUBLE_ARRAY_ARG(U),
                              DOUBLE_ARRAY_ARG(W),
                              DOUBLE_ARRAY_ARG(V),
                              DOUBLE_ARRAY_ARG(WORK), DOUBLE_ARRAY_ARG(X));

        for (i = 0; i < n; i++) {
            x[i] = X[i];
            work[i] = WORK[i];
        }
    }

    free(B);
    free(U);
    free(W);
    free(V);
    free(WORK);
    free(X);
}



F77_SUBROUTINE(sla_evp) (DOUBLE(DATE),
                         DOUBLE(DEQX),
                         DOUBLE_ARRAY(DVB),
                         DOUBLE_ARRAY(DPB), DOUBLE_ARRAY(DVH), DOUBLE_ARRAY(DPH));

void slaEvp(double date, double deqx, double dvb[3], double dpb[3],
            double dvh[3], double dph[3])
{
    DECLARE_DOUBLE(DATE);
    DECLARE_DOUBLE(DEQX);
    DECLARE_DOUBLE_ARRAY(DVB, 3);
    DECLARE_DOUBLE_ARRAY(DPB, 3);
    DECLARE_DOUBLE_ARRAY(DVH, 3);
    DECLARE_DOUBLE_ARRAY(DPH, 3);

    int i;
    DATE = date;
    DEQX = deqx;
    F77_CALL(sla_evp) (DOUBLE_ARG(&DATE),
                       DOUBLE_ARG(&DEQX),
                       DOUBLE_ARRAY_ARG(DVB),
                       DOUBLE_ARRAY_ARG(DPB),
                       DOUBLE_ARRAY_ARG(DVH), DOUBLE_ARRAY_ARG(DPH));
    for (i = 0; i < 3; i++) {
        dvb[i] = DVB[i];
        dpb[i] = DPB[i];
        dvh[i] = DVH[i];
        dph[i] = DPH[i];
    }

}

F77_SUBROUTINE(sla_fk5hz) (DOUBLE(R5),
                           DOUBLE(D5), DOUBLE(JEPOCH), DOUBLE(RH), DOUBLE(DH));

void slaFk5hz(double r5, double d5, double jepoch, double *rh, double *dh)
{
    DECLARE_DOUBLE(R5);
    DECLARE_DOUBLE(D5);
    DECLARE_DOUBLE(JEPOCH);
    DECLARE_DOUBLE(RH);
    DECLARE_DOUBLE(DH);
    R5 = r5;
    D5 = d5;
    JEPOCH = jepoch;
    F77_CALL(sla_fk5hz) (DOUBLE_ARG(&R5),
                         DOUBLE_ARG(&D5),
                         DOUBLE_ARG(&JEPOCH), DOUBLE_ARG(&RH), DOUBLE_ARG(&DH));
    *rh = RH;
    *dh = DH;
}

F77_SUBROUTINE(sla_hfk5z) (DOUBLE(RH),
                           DOUBLE(DH),
                           DOUBLE(JEPOCH),
                           DOUBLE(R5), DOUBLE(D5), DOUBLE(DR5), DOUBLE(DD5));

void slaHfk5z(double rh, double dh, double jepoch,
              double *r5, double *d5, double *dr5, double *dd5)
{
    DECLARE_DOUBLE(RH);
    DECLARE_DOUBLE(DH);
    DECLARE_DOUBLE(JEPOCH);
    DECLARE_DOUBLE(R5);
    DECLARE_DOUBLE(D5);
    DECLARE_DOUBLE(DR5);
    DECLARE_DOUBLE(DD5);
    RH = rh;
    DH = dh;
    JEPOCH = jepoch;
    F77_CALL(sla_hfk5z) (DOUBLE_ARG(&RH),
                         DOUBLE_ARG(&DH),
                         DOUBLE_ARG(&JEPOCH),
                         DOUBLE_ARG(&R5),
                         DOUBLE_ARG(&D5), DOUBLE_ARG(&DR5), DOUBLE_ARG(&DD5));
    *r5 = R5;
    *d5 = D5;
    *dr5 = DR5;
    *dd5 = DD5;
}

F77_SUBROUTINE(sla_geoc) (DOUBLE(P), DOUBLE(H), DOUBLE(R), DOUBLE(Z));

void slaGeoc(double p, double h, double *r, double *z)
{
    DECLARE_DOUBLE(P);
    DECLARE_DOUBLE(H);
    DECLARE_DOUBLE(R);
    DECLARE_DOUBLE(Z);
    P = p;
    H = h;
    F77_CALL(sla_geoc) (DOUBLE_ARG(&P),
                        DOUBLE_ARG(&H), DOUBLE_ARG(&R), DOUBLE_ARG(&Z));
    *r = R;
    *z = Z;
}


F77_SUBROUTINE(sla_deuler) (CHARACTER(ORDER),
                            DOUBLE(PHI),
                            DOUBLE(THETA), DOUBLE(PSI), DOUBLE_ARRAY(RMAT)
                            TRAIL(ORDER));

void slaDeuler(char *order, double phi, double theta, double psi, double rmat[3][3])
{

    DECLARE_CHARACTER(ORDER, 4);
    DECLARE_DOUBLE(PHI);
    DECLARE_DOUBLE(THETA);
    DECLARE_DOUBLE(PSI);
    DECLARE_DOUBLE_ARRAY(RMAT, 9);
    int i, j;

    PHI = phi;
    THETA = theta;
    PSI = psi;

    slaStringExport(order, ORDER, 4);

    F77_CALL(sla_deuler) (CHARACTER_ARG(ORDER),
                          DOUBLE_ARG(&PHI),
                          DOUBLE_ARG(&THETA),
                          DOUBLE_ARG(&PSI), DOUBLE_ARRAY_ARG(RMAT)
                          TRAIL_ARG(ORDER));

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            rmat[i][j] = RMAT[i + 3 * j];
    }

}


F77_SUBROUTINE(sla_de2h) (DOUBLE(HA),
                          DOUBLE(DEC), DOUBLE(PHI), DOUBLE(AZ), DOUBLE(EL));

void slaDe2h(double ha, double dec, double phi, double *az, double *el)
{
    DECLARE_DOUBLE(HA);
    DECLARE_DOUBLE(DEC);
    DECLARE_DOUBLE(PHI);
    DECLARE_DOUBLE(AZ);
    DECLARE_DOUBLE(EL);
    HA = ha;
    DEC = dec;
    PHI = phi;
    F77_CALL(sla_de2h) (DOUBLE_ARG(&HA),
                        DOUBLE_ARG(&DEC),
                        DOUBLE_ARG(&PHI), DOUBLE_ARG(&AZ), DOUBLE_ARG(&EL));
    *az = AZ;
    *el = EL;
}

F77_SUBROUTINE(sla_dh2e) (DOUBLE(AZ),
                          DOUBLE(EL), DOUBLE(PHI), DOUBLE(HA), DOUBLE(DEC));

void slaDh2e(double az, double el, double phi, double *ha, double *dec)
{
    DECLARE_DOUBLE(AZ);
    DECLARE_DOUBLE(EL);
    DECLARE_DOUBLE(PHI);
    DECLARE_DOUBLE(HA);
    DECLARE_DOUBLE(DEC);
    AZ = az;
    EL = el;
    PHI = phi;
    F77_CALL(sla_dh2e) (DOUBLE_ARG(&AZ),
                        DOUBLE_ARG(&EL),
                        DOUBLE_ARG(&PHI), DOUBLE_ARG(&HA), DOUBLE_ARG(&DEC));
    *ha = HA;
    *dec = DEC;
}


F77_SUBROUTINE(sla_oap) (CHARACTER(TYPE),
                         DOUBLE(OB1),
                         DOUBLE(OB2),
                         DOUBLE(DATE),
                         DOUBLE(DUT),
                         DOUBLE(ELONGM),
                         DOUBLE(PHIM),
                         DOUBLE(HM),
                         DOUBLE(XP),
                         DOUBLE(YP),
                         DOUBLE(TDK),
                         DOUBLE(PMB),
                         DOUBLE(RH),
                         DOUBLE(WL), DOUBLE(TLR), DOUBLE(RAP), DOUBLE(DAP)
                         TRAIL(TYPE));

/* Note that SLA insists that "c" has space for 10 characters + nul
   and "name" has space for 40 characters + nul */

void
slaOap(char *type, double ob1, double ob2, double date,
       double dut, double elongm, double phim, double hm,
       double xp, double yp, double tdk, double pmb,
       double rh, double wl, double tlr, double *rap, double *dap)
{

    DECLARE_CHARACTER(TYPE, 1);
    DECLARE_DOUBLE(OB1);
    DECLARE_DOUBLE(OB2);
    DECLARE_DOUBLE(DATE);
    DECLARE_DOUBLE(DUT);
    DECLARE_DOUBLE(ELONGM);
    DECLARE_DOUBLE(PHIM);
    DECLARE_DOUBLE(HM);
    DECLARE_DOUBLE(XP);
    DECLARE_DOUBLE(YP);
    DECLARE_DOUBLE(TDK);
    DECLARE_DOUBLE(PMB);
    DECLARE_DOUBLE(RH);
    DECLARE_DOUBLE(WL);
    DECLARE_DOUBLE(TLR);
    DECLARE_DOUBLE(RAP);
    DECLARE_DOUBLE(DAP);
    OB1 = ob1;
    OB2 = ob2;
    DATE = date;
    DUT = dut;
    ELONGM = elongm;
    PHIM = phim;
    HM = hm;
    XP = xp;
    YP = yp;
    TDK = tdk;
    PMB = pmb;
    RH = rh;
    WL = wl;
    TLR = tlr;
    slaStringExport(type, TYPE, 1);

    /* call the routine */
    F77_CALL(sla_oap) (CHARACTER_ARG(TYPE),
                       DOUBLE_ARG(&OB1),
                       DOUBLE_ARG(&OB2),
                       DOUBLE_ARG(&DATE),
                       DOUBLE_ARG(&DUT),
                       DOUBLE_ARG(&ELONGM),
                       DOUBLE_ARG(&PHIM),
                       DOUBLE_ARG(&HM),
                       DOUBLE_ARG(&XP),
                       DOUBLE_ARG(&YP),
                       DOUBLE_ARG(&TDK),
                       DOUBLE_ARG(&PMB),
                       DOUBLE_ARG(&RH),
                       DOUBLE_ARG(&WL),
                       DOUBLE_ARG(&TLR), DOUBLE_ARG(&RAP), DOUBLE_ARG(&DAP)
                       TRAIL_ARG(TYPE));
    *rap = RAP;
    *dap = DAP;
}

F77_SUBROUTINE(sla_obs) (INTEGER(I),
                         CHARACTER(C),
                         CHARACTER(NAME), DOUBLE(W), DOUBLE(P), DOUBLE(H)
                         TRAIL(C)
                         TRAIL(NAME));

/* Note that SLA insists that "c" has space for 10 characters + nul
   and "name" has space for 40 characters + nul */

void slaObs(int n, char *c, char *name, double *w, double *p, double *h)
{

    DECLARE_INTEGER(N);
    DECLARE_CHARACTER(C, 10);
    DECLARE_CHARACTER(NAME, 40);
    DECLARE_DOUBLE(W);
    DECLARE_DOUBLE(P);
    DECLARE_DOUBLE(H);

    if (n < 1) {
        /* C needs to be imported */
        slaStringExport(c, C, 10);
    } else {
        /* initialise C */
        slaStringExport("", C, 10);
    }
    F77_EXPORT_INTEGER(n, N);

    /* w, p and h are not touched on error but for consistency this means
       we copy the current values to Fortran so that we can correctly copy
       back the result. */
    F77_EXPORT_DOUBLE(*w, W);
    F77_EXPORT_DOUBLE(*p, P);
    F77_EXPORT_DOUBLE(*h, H);

    /* call the routine */
    F77_CALL(sla_obs) (INTEGER_ARG(&N),
                       CHARACTER_ARG(C),
                       CHARACTER_ARG(NAME),
                       DOUBLE_ARG(&W), DOUBLE_ARG(&P), DOUBLE_ARG(&H)
                       TRAIL_ARG(C)
                       TRAIL_ARG(NAME));

    /* extract results */
    slaStringImport(NAME, 40, name);
    if (n > 0 && name[0] != '?') {
        /* only do this if we know we used a numeric input and if the result
           for the NAME is not '?' (since we are not allowed to alter the string
           in that case). This allows people
           to call slaObs with a string constant */
        slaStringImport(C, 10, c);
    }
    F77_IMPORT_DOUBLE(W, *w);
    F77_IMPORT_DOUBLE(P, *p);
    F77_IMPORT_DOUBLE(H, *h);

}

F77_DOUBLE_FUNCTION(sla_pa) (DOUBLE(HA), DOUBLE(DEC), DOUBLE(PHI));

double slaPa(double ha, double dec, double phi)
{
    DECLARE_DOUBLE(HA);
    DECLARE_DOUBLE(DEC);
    DECLARE_DOUBLE(PHI);
    DECLARE_DOUBLE(RETVAL);
    double retval;

    F77_EXPORT_DOUBLE(ha, HA);
    F77_EXPORT_DOUBLE(dec, DEC);
    F77_EXPORT_DOUBLE(phi, PHI);

    RETVAL = F77_CALL(sla_pa) (DOUBLE_ARG(&HA), DOUBLE_ARG(&DEC), DOUBLE_ARG(&PHI));

    F77_IMPORT_DOUBLE(RETVAL, retval);
    return retval;
}

F77_DOUBLE_FUNCTION(sla_dtt) (DOUBLE(UTC));

double slaDtt(double utc)
{
    DECLARE_DOUBLE(UTC);
    DECLARE_DOUBLE(RETVAL);
    double retval;

    F77_EXPORT_DOUBLE(utc, UTC);
    RETVAL = F77_CALL(sla_dtt) (DOUBLE_ARG(&UTC));

    F77_IMPORT_DOUBLE(RETVAL, retval);
    return retval;
}

F77_DOUBLE_FUNCTION(sla_dat) (DOUBLE(UTC));

double slaDat(double utc)
{
    DECLARE_DOUBLE(UTC);
    DECLARE_DOUBLE(RETVAL);
    double retval;

    F77_EXPORT_DOUBLE(utc, UTC);
    RETVAL = F77_CALL(sla_dat) (DOUBLE_ARG(&UTC));

    F77_IMPORT_DOUBLE(RETVAL, retval);
    return retval;
}

F77_SUBROUTINE(sla_rdplan) (DOUBLE(DATE), INTEGER(I), DOUBLE(ELONG), DOUBLE(PHI),
                            DOUBLE(RA), DOUBLE(DEC), DOUBLE(DIAM));

void
slaRdplan(double date, int i, double elong, double phi,
          double *ra, double *dec, double *diam)
{
    DECLARE_DOUBLE(DATE);
    DECLARE_INTEGER(I);
    DECLARE_DOUBLE(ELONG);
    DECLARE_DOUBLE(PHI);
    DECLARE_DOUBLE(RA);
    DECLARE_DOUBLE(DEC);
    DECLARE_DOUBLE(DIAM);

    F77_EXPORT_DOUBLE(date, DATE);
    F77_EXPORT_INTEGER(i, I);
    F77_EXPORT_DOUBLE(elong, ELONG);
    F77_EXPORT_DOUBLE(phi, PHI);

    F77_CALL(sla_rdplan) (DOUBLE_ARG(&DATE),
                          INTEGER_ARG(&I),
                          DOUBLE_ARG(&ELONG),
                          DOUBLE_ARG(&PHI),
                          DOUBLE_ARG(&RA), DOUBLE_ARG(&DEC), DOUBLE_ARG(&DIAM));

    F77_IMPORT_DOUBLE(RA, *ra);
    F77_IMPORT_DOUBLE(DEC, *dec);
    F77_IMPORT_DOUBLE(DIAM, *diam);
}
