/*
*+
*  Name:
*     f77.h and cnf.h

*  Purpose:
*     C - FORTRAN interace macros and prototypes

*  Language:
*     C (part ANSI, part not)

*  Type of Module:
*     C include file

*  Description:
*     Including this file in SLALIB allows SLALIB to remain free of any 
*     dependencies on the rest of the Starlink Software Collection.
*
*     For historical reasons two files, F77.h and cnf.h are required
*     but the have now been combined and for new code, only one is 
*     necessary.
*
*     This file defines the macros needed to write C functions that are
*     designed to be called from FORTRAN programs, and to do so in a
*     portable way. Arguments are normally passed by reference from a
*     FORTRAN program, and so the F77 macros arrange for a pointer to
*     all arguments to be available. This requires no work on most
*     machines, but will actually generate the pointers on a machine
*     that passes FORTRAN arguments by value.

*  Notes:
*     -  Macros are provided to handle the conversion of logical data
*        values between the way that FORTRAN represents a value and the
*        way that C represents it.
*     -  Macros are provided to convert variables between the FORTRAN and
*        C method of representing them. In most cases there is no
*        conversion required, the macros just arrange for a pointer to
*        the FORTRAN variable to be set appropriately. The possibility that
*        FORTRAN and C might use different ways of representing integer
*        and floating point values is considered remote, the macros are
*        really only there for completeness and to assist in the automatic 
*        generation of C interfaces.
*     -  For character variables the macros convert between
*        the FORTRAN method of representing them (fixed length, blank
*        filled strings) and the C method (variable length, null
*        terminated strings) using calls to the CNF functions.

*  Implementation Deficiencies:
*     -  The macros support the K&R style of function definition, but
*        this file may not work with all K&R compilers as it contains
*        "#if defined" statements. These could be replaced with #ifdef's
*        if necessary. This has not been done as is would make the code
*        less clear and the need for support for K&R sytle definitions
*        should disappear as ANSI compilers become the default.

*  Copyright:
*     Copyright (C) 1991, 1993 Science & Engineering Research Council
*     Copyright (C) 2006 Particle Physics and Astronomy Research Council

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*     
*     This program is distributed in the hope that it will be
*     useful,but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*     
*     You should have received a copy of the GNU General Public License
*     along with this program; if not, write to the Free Software
*     Foundation, Inc., 59 Temple Place,Suite 330, Boston, MA
*     02111-1307, USA

*  Authors:
*     PMA: Peter Allan (Starlink, RAL)
*     AJC: Alan Chipperfield (Starlink, RAL)
*     DSB: David S Berry (JAC)
*     PWD: Peter W. Draper (JAC, Durham University)
*     {enter_new_authors_here}

*  History:
*     13-JUN-2006 (DSB):
*        Original version, copied from AST.
*     22-JUN-2006 (PWD):
*        Updated to include changes that handle the return type of 
*        REAL functions (double, not float, on some platforms).
*     {enter_further_changes_here}
*        

*  Bugs:
*     {note_any_bugs_here}

*-
------------------------------------------------------------------------------
*/
#if !defined(CNF_MACROS)
#define CNF_MACROS

#include <stdlib.h>
/*  This initial sections defines values for all macros. These are the	    */
/*  values that are generally appropriate to an ANSI C compiler on Unix.    */
/*  For macros that have different values on other systems, the macros	    */
/*  should be undefined and then redefined in the system specific sections. */
/*  At the end of this section, some macros are redefined if the compiler   */
/*  is non-ANSI.							    */


#if defined(__STDC__) || defined(VMS)
#define CNF_CONST const
#else
#define CNF_CONST
#endif

/*  -----  Macros common to calling C from FORTRAN and FORTRAN from C  ---- */


/*  ---  External Names  ---						    */

/*  Macro to define the name of a Fortran routine or common block. This	    */
/*  ends in an underscore on many Unix systems.				    */

#define F77_EXTERNAL_NAME(X) X ## _


/*  ---  Logical Values  ---						    */

/*  Define the values that are used to represent the logical values TRUE    */
/*  and FALSE in Fortran.						    */

#define F77_TRUE 1
#define F77_FALSE 0

/*  Define macros that evaluate to C logical values, given a FORTRAN	    */
/*  logical value.							    */

#define F77_ISTRUE(X) ( X )
#define F77_ISFALSE(X) ( !( X ) )


/*  ---  Common Blocks  ---						    */

/*  Macros used in referring to FORTRAN common blocks.			    */

#define F77_BLANK_COMMON                _BLNK__
#define F77_NAMED_COMMON(B)             F77_EXTERNAL_NAME(B)



/*  ------------------  Calling C from FORTRAN  --------------------------- */


/*  ---  Data Types  ---						    */

/*  Define macros for all the Fortran data types (except COMPLEX, which is  */
/*  not handled by this package).					    */

#define F77_INTEGER_TYPE   int 
#define F77_REAL_TYPE      float 
#define F77_REAL_FUNCTION_TYPE      float
#define F77_DOUBLE_TYPE    double 
#define F77_LOGICAL_TYPE   int
#define F77_CHARACTER_TYPE char
#define F77_BYTE_TYPE      signed char
#define F77_WORD_TYPE      short int
#define F77_UBYTE_TYPE     unsigned char
#define F77_UWORD_TYPE     unsigned short int

/*  Define macros for the type of a CHARACTER and CHARACTER_ARRAY argument  */
#define F77_CHARACTER_ARG_TYPE char
#define F77_CHARACTER_ARRAY_ARG_TYPE char

/*  Define a macro to use when passing arguments that STARLINK FORTRAN	    */
/*  treats as a pointer. From the point of view of C, this type should be   */
/*  (void *), but it is declared as type unsigned int as we actually pass   */
/*  an INTEGER from the FORTRAN routine. The distinction is important for   */
/*  architectures where the size of an INTEGER is not the same as the size  */
/*  of a pointer.							    */

#define F77_POINTER_TYPE   unsigned int


/*  ---  Subroutine Names  ---						    */

/* This declares that the C function returns a value of void.		    */

#define F77_SUBROUTINE(X)  void F77_EXTERNAL_NAME(X)


/*  ---  Function Names  ---						    */

/*  Macros to define the types and names of functions that return values.   */
/*  Due the the different ways that function return values could be	    */
/*  implemented, it is better not to use functions, but to stick to using   */
/*  subroutines.							    */

/*  Character functions are implemented, but in a way that cannot be	    */
/*  guaranteed to be portable although it will work on VMS, SunOS, Ultrix   */
/*  and DEC OSF/1. It would be better to return the character value as a    */
/*  subroutine argument where possible, rather than use a character	    */
/*  function.								    */

#define F77_INTEGER_FUNCTION(X)   F77_INTEGER_TYPE F77_EXTERNAL_NAME(X)
#define F77_REAL_FUNCTION(X)      F77_REAL_FUNCTION_TYPE F77_EXTERNAL_NAME(X)
#define F77_DOUBLE_FUNCTION(X)    F77_DOUBLE_TYPE F77_EXTERNAL_NAME(X)
#define F77_LOGICAL_FUNCTION(X)   F77_LOGICAL_TYPE F77_EXTERNAL_NAME(X)
#define F77_CHARACTER_FUNCTION(X) void F77_EXTERNAL_NAME(X)
#define F77_BYTE_FUNCTION(X)      F77_BYTE_TYPE F77_EXTERNAL_NAME(X)
#define F77_WORD_FUNCTION(X)      F77_WORD_TYPE F77_EXTERNAL_NAME(X)
#define F77_UBYTE_FUNCTION(X)     F77_UBYTE_TYPE F77_EXTERNAL_NAME(X)
#define F77_UWORD_FUNCTION(X)     F77_UWORD_TYPE F77_EXTERNAL_NAME(X)
#define F77_POINTER_FUNCTION(X)   F77_POINTER_TYPE F77_EXTERNAL_NAME(X)


/*  ---  Character return value for a function  ---			    */

#define CHARACTER_RETURN_VALUE(X) CHARACTER(X) TRAIL(X)
#define CHARACTER_RETURN_ARG(X) CHARACTER_ARG(X) TRAIL_ARG(X)

/*  ---  Dummy Arguments  ---						    */

/*  Macros for defining subroutine arguments. All these macros take a	    */
/*  single argument; the name of the parameter. On most systems, a numeric  */
/*  argument is passed as a pointer.					    */

#define INTEGER(X)     F77_INTEGER_TYPE *CNF_CONST X
#define REAL(X)        F77_REAL_TYPE    *CNF_CONST X
#define DOUBLE(X)      F77_DOUBLE_TYPE  *CNF_CONST X
#define LOGICAL(X)     F77_LOGICAL_TYPE *CNF_CONST X
#define BYTE(X)        F77_BYTE_TYPE    *CNF_CONST X
#define WORD(X)        F77_WORD_TYPE    *CNF_CONST X
#define UBYTE(X)       F77_UBYTE_TYPE   *CNF_CONST X
#define UWORD(X)       F77_UWORD_TYPE   *CNF_CONST X

/*  Pointer arguments. Define a pointer type for passing pointer values	    */
/*  between subroutines.						    */

#define POINTER(X)     F77_POINTER_TYPE *CNF_CONST X

/*  EXTERNAL arguments. Define a passed subroutine or function name */
#define SUBROUTINE(X)  void (*X)()
#define INTEGER_FUNCTION(X)  F77_INTEGER_TYPE (*X)()
#define REAL_FUNCTION(X)  F77_REAL_TYPE (*X)()
#define DOUBLE_FUNCTION(X)  F77_DOUBLE_TYPE (*X)()
#define LOGICAL_FUNCTION(X)  F77_LOGICAL_TYPE (*X)()
#define CHARACTER_FUNCTION(X)  F77_CHARACTER_TYPE (*X)()
#define BYTE_FUNCTION(X)  F77_BYTE_TYPE (*X)()
#define WORD_FUNCTION(X)  F77_WORD_TYPE (*X)()
#define UBYTE_FUNCTION(X)  F77_UBYTE_TYPE (*X)()
#define UWORD_FUNCTION(X)  F77_UWORD_TYPE (*X)()
#define POINTER_FUNCTION(X)  F77_POINTER_TYPE (*X)()

/*  Array arguments.							    */

#define INTEGER_ARRAY(X)     F77_INTEGER_TYPE *CNF_CONST X
#define REAL_ARRAY(X)        F77_REAL_TYPE    *CNF_CONST X
#define DOUBLE_ARRAY(X)      F77_DOUBLE_TYPE  *CNF_CONST X
#define LOGICAL_ARRAY(X)     F77_LOGICAL_TYPE *CNF_CONST X
#define BYTE_ARRAY(X)        F77_BYTE_TYPE    *CNF_CONST X
#define WORD_ARRAY(X)        F77_WORD_TYPE    *CNF_CONST X
#define UBYTE_ARRAY(X)       F77_UBYTE_TYPE   *CNF_CONST X
#define UWORD_ARRAY(X)       F77_UWORD_TYPE   *CNF_CONST X

#define POINTER_ARRAY(X)     F77_POINTER_TYPE *CNF_CONST X

/*  Macros to handle character arguments.				    */

/*  Character arguments can be passed in many ways. The purpose of these    */
/*  macros and the GENPTR_CHARACTER macro (defined in the next section) is  */
/*  to generate a pointer to a character variable called ARG and an integer */
/*  ARG_length containing the length of ARG. If these two variables are	    */
/*  available directly from the argument list of the routine, then the	    */
/*  GENPTR_CHARACTER macro is null, otherwise it works on intermediate	    */
/*  variables.								    */

#define CHARACTER(X)             F77_CHARACTER_TYPE *CNF_CONST X
#define TRAIL(X)                 ,int X ## _length
#define CHARACTER_ARRAY(X)       F77_CHARACTER_TYPE *CNF_CONST X


/*  ---  Getting Pointers to Arguments  ---				    */

/*  Macros that ensure that a pointer to each argument is available for the */
/*  programmer to use. Usually this means that these macros are null. On    */
/*  VMS, a pointer to a character variable has to be generated. If a	    */
/*  particular machine were to pass arguments by reference, rather than by  */
/*  value, then these macros would construct the appropriate pointers.	    */

#define GENPTR_INTEGER(X)
#define GENPTR_REAL(X)
#define GENPTR_DOUBLE(X)
#define GENPTR_CHARACTER(X)
#define GENPTR_LOGICAL(X)
#define GENPTR_BYTE(X)
#define GENPTR_WORD(X)
#define GENPTR_UBYTE(X)
#define GENPTR_UWORD(X)
#define GENPTR_POINTER(X)

#define GENPTR_INTEGER_ARRAY(X)
#define GENPTR_REAL_ARRAY(X)
#define GENPTR_DOUBLE_ARRAY(X)
#define GENPTR_CHARACTER_ARRAY(X)
#define GENPTR_LOGICAL_ARRAY(X)
#define GENPTR_BYTE_ARRAY(X)
#define GENPTR_WORD_ARRAY(X)
#define GENPTR_UBYTE_ARRAY(X)
#define GENPTR_UWORD_ARRAY(X)
#define GENPTR_POINTER_ARRAY(X)

#define GENPTR_SUBROUTINE(X)
#define GENPTR_INTEGER_FUNCTION(X)
#define GENPTR_REAL_FUNCTION(X)
#define GENPTR_DOUBLE_FUNCTION(X)
#define GENPTR_CHARACTER_FUNCTION(X)
#define GENPTR_LOGICAL_FUNCTION(X)
#define GENPTR_BYTE_FUNCTION(X)
#define GENPTR_WORD_FUNCTION(X)
#define GENPTR_UBYTE_FUNCTION(X)
#define GENPTR_UWORD_FUNCTION(X)
#define GENPTR_POINTER_FUNCTION(X)



/*  ------------------  Calling FORTRAN from C  --------------------------- */


/*  ---  Declare variables  ---						    */

#define DECLARE_INTEGER(X) F77_INTEGER_TYPE X
#define DECLARE_REAL(X)    F77_REAL_TYPE X
#define DECLARE_DOUBLE(X)  F77_DOUBLE_TYPE X
#define DECLARE_LOGICAL(X) F77_LOGICAL_TYPE X
#define DECLARE_BYTE(X)    F77_BYTE_TYPE X
#define DECLARE_WORD(X)    F77_WORD_TYPE X
#define DECLARE_UBYTE(X)   F77_UBYTE_TYPE X
#define DECLARE_UWORD(X)   F77_UWORD_TYPE X

#define DECLARE_POINTER(X) F77_POINTER_TYPE X

#define DECLARE_CHARACTER(X,L) F77_CHARACTER_TYPE X[L]; \
   const int X##_length = L


/*  ---  Declare arrays ---						    */

#define DECLARE_INTEGER_ARRAY(X,D) F77_INTEGER_TYPE X[D]
#define DECLARE_REAL_ARRAY(X,D)    F77_REAL_TYPE X[D]
#define DECLARE_DOUBLE_ARRAY(X,D)  F77_DOUBLE_TYPE X[D]
#define DECLARE_LOGICAL_ARRAY(X,D) F77_LOGICAL_TYPE X[D]
#define DECLARE_BYTE_ARRAY(X,D)    F77_BYTE_TYPE X[D]
#define DECLARE_WORD_ARRAY(X,D)    F77_WORD_TYPE X[D]
#define DECLARE_UBYTE_ARRAY(X,D)   F77_UBYTE_TYPE X[D]
#define DECLARE_UWORD_ARRAY(X,D)   F77_UWORD_TYPE X[D]
#define DECLARE_POINTER_ARRAY(X,D) F77_POINTER_TYPE X[D]
#define DECLARE_CHARACTER_ARRAY(X,L,D) F77_CHARACTER_TYPE X[D][L]; \
   const int X##_length = L

/*  ---  Declare and construct dynamic CHARACTER arguments ---                      */
#define DECLARE_CHARACTER_DYN(X)   F77_CHARACTER_TYPE *X;\
   int X##_length
#define F77_CREATE_CHARACTER(X,L)  X=cnfCref(L);\
   X##_length = L

/* Declare Dynamic Fortran arrays */
#define DECLARE_INTEGER_ARRAY_DYN(X) F77_INTEGER_TYPE *X
#define DECLARE_REAL_ARRAY_DYN(X)    F77_REAL_TYPE *X
#define DECLARE_DOUBLE_ARRAY_DYN(X)  F77_DOUBLE_TYPE *X
#define DECLARE_LOGICAL_ARRAY_DYN(X) F77_LOGICAL_TYPE *X
#define DECLARE_BYTE_ARRAY_DYN(X)    F77_BYTE_TYPE *X
#define DECLARE_WORD_ARRAY_DYN(X)    F77_WORD_TYPE *X
#define DECLARE_UBYTE_ARRAY_DYN(X)   F77_UBYTE_TYPE *X
#define DECLARE_UWORD_ARRAY_DYN(X)   F77_UWORD_TYPE *X
#define DECLARE_POINTER_ARRAY_DYN(X) F77_POINTER_TYPE *X
#define DECLARE_CHARACTER_ARRAY_DYN(X)   F77_CHARACTER_TYPE *X;\
   int X##_length

/* Create arrays dynamic Fortran arrays for those types which require */
/* Separate space for Fortran and C arrays                            */
/* Character and logical are already defined                          */
/* For most types there is nothing to do                              */
#define F77_CREATE_CHARACTER_ARRAY(X,L,N) \
        {int f77dims[1];f77dims[0]=N;X=cnfCrefa(L,1,f77dims);X##_length=L;}
#define F77_CREATE_CHARACTER_ARRAY_M(X,L,N,D)  X=cnfCrefa(L,N,D);\
   X##_length = L
#define F77_CREATE_LOGICAL_ARRAY(X,N) \
        {int f77dims[1];f77dims[0]=N;X=cnfCrela(1,f77dims);}
#define F77_CREATE_LOGICAL_ARRAY_M(X,N,D) X=cnfCrela(N,D)
#define F77_CREATE_INTEGER_ARRAY(X,N)
#define F77_CREATE_REAL_ARRAY(X,N)
#define F77_CREATE_DOUBLE_ARRAY(X,N)
#define F77_CREATE_BYTE_ARRAY(X,N)
#define F77_CREATE_UBYTE_ARRAY(X,N)
#define F77_CREATE_WORD_ARRAY(X,N)
#define F77_CREATE_UWORD_ARRAY(X,N)
#define F77_CREATE_POINTER_ARRAY(X,N)\
        X=(F77_POINTER_TYPE *) malloc(N*sizeof(F77_POINTER_TYPE))

/* Associate Fortran arrays with C arrays                             */
/* These macros ensure that there is space somewhere for the Fortran  */
/* array. They are complemetary to the CREATE_type_ARRAY macros       */
#define F77_ASSOC_CHARACTER_ARRAY(F,C)
#define F77_ASSOC_LOGICAL_ARRAY(F,C)
#define F77_ASSOC_INTEGER_ARRAY(F,C) F=C
#define F77_ASSOC_REAL_ARRAY(F,C) F=C
#define F77_ASSOC_DOUBLE_ARRAY(F,C) F=C
#define F77_ASSOC_BYTE_ARRAY(F,C) F=C
#define F77_ASSOC_UBYTE_ARRAY(F,C) F=C
#define F77_ASSOC_WORD_ARRAY(F,C) F=C
#define F77_ASSOC_UWORD_ARRAY(F,C) F=C
#define F77_ASSOC_POINTER_ARRAY(F,C)

/* Free created dynamic arrays */
/* Character and logical are already defined */
/* For most types there is nothing to do */
#define F77_FREE_INTEGER(X)
#define F77_FREE_REAL(X)
#define F77_FREE_DOUBLE(X)
#define F77_FREE_BYTE(X)
#define F77_FREE_UBYTE(X)
#define F77_FREE_WORD(X)
#define F77_FREE_UWORD(X)
#define F77_FREE_POINTER(X) cnfFree((void *)X);
#define F77_FREE_CHARACTER(X) cnfFreef( X )
#define F77_FREE_LOGICAL(X) cnfFree( (char *)X )

/*  ---  IMPORT and EXPORT of values  --- */
/* Export C variables to Fortran variables */
#define F77_EXPORT_CHARACTER(C,F,L) cnfExprt(C,F,L)
#define F77_EXPORT_DOUBLE(C,F) F=C
#define F77_EXPORT_INTEGER(C,F) F=C
#define F77_EXPORT_LOGICAL(C,F) F=C?F77_TRUE:F77_FALSE
#define F77_EXPORT_REAL(C,F) F=C
#define F77_EXPORT_BYTE(C,F) F=C
#define F77_EXPORT_WORD(C,F) F=C
#define F77_EXPORT_UBYTE(C,F) F=C
#define F77_EXPORT_UWORD(C,F) F=C
#define F77_EXPORT_POINTER(C,F) F=cnfFptr(C)
#define F77_EXPORT_LOCATOR(C,F) cnfExpch(C,F,DAT__SZLOC)

/* Export C arrays to Fortran */
/* Arrays are assumed to be 1-d so just the number of elements is given  */
/* This may be OK for n-d arrays also */
/* CHARACTER arrays may be represented in C as arrays of arrays of char or */
/* as arrays of pointers to char (the _P variant) */
#define F77_EXPORT_CHARACTER_ARRAY(C,LC,F,LF,N) \
        {int f77dims[1];f77dims[0]=N;cnfExprta(C,LC,F,LF,1,f77dims);}
#define F77_EXPORT_CHARACTER_ARRAY_P(C,F,LF,N) \
        {int f77dims[1];f77dims[0]=N;cnfExprtap(C,F,LF,1,f77dims);}
#define F77_EXPORT_DOUBLE_ARRAY(C,F,N) F=(F77_DOUBLE_TYPE *)C
#define F77_EXPORT_INTEGER_ARRAY(C,F,N) F=(F77_INTEGER_TYPE *)C
#define F77_EXPORT_LOGICAL_ARRAY(C,F,N) \
        {int f77dims[1];f77dims[0]=N;cnfExpla(C,F,1,f77dims);}
#define F77_EXPORT_REAL_ARRAY(C,F,N) F=(F77_REAL_TYPE *)C
#define F77_EXPORT_BYTE_ARRAY(C,F,N) F=(F77_BYTE_TYPE *)C
#define F77_EXPORT_WORD_ARRAY(C,F,N) F=(F77_WORD_TYPE *)C
#define F77_EXPORT_UBYTE_ARRAY(C,F,N) F=(F77_UBYTE_TYPE *)C
#define F77_EXPORT_UWORD_ARRAY(C,F,N) F=(F77_UWORD_TYPE * )C
#define F77_EXPORT_POINTER_ARRAY(C,F,N) \
     {int f77i;for (f77i=0;f77i<N;f77i++)F[f77i]=cnfFptr(C[f77i]);}
#define F77_EXPORT_LOCATOR_ARRAY(C,F,N) \
        {int f77i;for (f77i=0;f77i<N;f77i++)cnfExpch(C,F,DAT__SZLOC);}

/* Import Fortran variables to C */
#define F77_IMPORT_CHARACTER(F,L,C) cnfImprt(F,L,C)
#define F77_IMPORT_DOUBLE(F,C) C=F
#define F77_IMPORT_INTEGER(F,C) C=F
#define F77_IMPORT_LOGICAL(F,C) C=F77_ISTRUE(F)
#define F77_IMPORT_REAL(F,C) C=F
#define F77_IMPORT_BYTE(F,C) C=F
#define F77_IMPORT_WORD(F,C) C=F
#define F77_IMPORT_UBYTE(F,C) C=F
#define F77_IMPORT_UWORD(F,C) C=F
#define F77_IMPORT_POINTER(F,C) C=cnfCptr(F)
#define F77_IMPORT_LOCATOR(F,C) cnfImpch(F,DAT__SZLOC,C)

/* Import Fortran arrays to C */
/* Arrays are assumed to be 1-d so just the number of elements is given  */
/* This may be OK for n-d arrays also */
/* CHARACTER arrays may be represented in C as arrays of arrays of char or */
/* as arrays of pointers to char (the _P variant) */
#define F77_IMPORT_CHARACTER_ARRAY(F,LF,C,LC,N) \
        {int f77dims[1];f77dims[0]=N;cnfImprta(F,LF,C,LC,1,f77dims);}
#define F77_IMPORT_CHARACTER_ARRAY_P(F,LF,C,LC,N) \
        {int f77dims[1];f77dims[0]=N;cnfImprtap(F,LF,C,LC,1,f77dims);}
#define F77_IMPORT_DOUBLE_ARRAY(F,C,N)
#define F77_IMPORT_INTEGER_ARRAY(F,C,N)
#define F77_IMPORT_LOGICAL_ARRAY(F,C,N) \
        {int f77dims[1];f77dims[0]=N;cnfImpla(F,C,1,f77dims);}
#define F77_IMPORT_REAL_ARRAY(F,C,N)
#define F77_IMPORT_BYTE_ARRAY(F,C,N)
#define F77_IMPORT_WORD_ARRAY(F,C,N)
#define F77_IMPORT_UBYTE_ARRAY(F,C,N)
#define F77_IMPORT_UWORD_ARRAY(F,C,N)
#define F77_IMPORT_POINTER_ARRAY(F,C,N) \
        {int f77i;for (f77i=0;f77i<N;f77i++)C[f77i]=cnfCptr(F[f77i]);}
#define F77_IMPORT_LOCATOR_ARRAY(F,C,N) \
        {int f77i;for (f77i=0;f77i<N;f77i++)cnfImpch(F,DAT__SZLOC,C);}

/*  ---  Call a FORTRAN routine  ---					    */

#define F77_CALL(X)  F77_EXTERNAL_NAME(X)


/*  ---  Pass arguments to a FORTRAN routine  ---			    */

#define INTEGER_ARG(X)   X
#define REAL_ARG(X)      X
#define DOUBLE_ARG(X)    X
#define LOGICAL_ARG(X)   X
#define BYTE_ARG(X)      X
#define WORD_ARG(X)      X
#define UBYTE_ARG(X)     X
#define UWORD_ARG(X)     X
#define POINTER_ARG(X)   X
#define CHARACTER_ARG(X) X
#define TRAIL_ARG(X)     ,X##_length

#define SUBROUTINE_ARG(X)  X
#define INTEGER_FUNCTION_ARG(X)  X
#define REAL_FUNCTION_ARG(X)  X
#define DOUBLE_FUNCTION_ARG(X)  X
#define LOGICAL_FUNCTION_ARG(X)  X
#define CHARACTER_FUNCTION_ARG(X)  X
#define BYTE_FUNCTION_ARG(X)  X
#define WORD_FUNCTION_ARG(X)  X
#define UBYTE_FUNCTION_ARG(X)  X
#define UWORD_FUNCTION_ARG(X)  X
#define POINTER_FUNCTION_ARG(X)  X

#define INTEGER_ARRAY_ARG(X)   (F77_INTEGER_TYPE *)X
#define REAL_ARRAY_ARG(X)      (F77_REAL_TYPE *)X
#define DOUBLE_ARRAY_ARG(X)    (F77_DOUBLE_TYPE *)X
#define LOGICAL_ARRAY_ARG(X)   (F77_LOGICAL_TYPE *)X
#define BYTE_ARRAY_ARG(X)      (F77_BYTE_TYPE *)X
#define WORD_ARRAY_ARG(X)      (F77_WORD_TYPE *)X
#define UBYTE_ARRAY_ARG(X)     (F77_UBYTE_TYPE *)X
#define UWORD_ARRAY_ARG(X)     (F77_UWORD_TYPE *)X
#define POINTER_ARRAY_ARG(X)   (F77_POINTER_TYPE *)X
#define CHARACTER_ARRAY_ARG(X) (F77_CHARACTER_ARRAY_ARG_TYPE *)X


/* ------------------------ Non-ansi section ------------------------------ */

/*  The difference between ANSI and non-ANSI compilers, as far as macro	    */
/*  definition is concerned, is that non-ANSI compilers do not support the  */
/*  token concatenation operator (##). To work around this, we use the fact */
/*  that the null comment is preprocessed to produce no characters at all   */
/*  by our non-ANSI compilers.						    */
/*  This section does not deal with the fact that some non-ANSI compilers   */
/*  cannot handle function prototypes. That is handled in the machine 	    */
/*  specific sections.							    */

#if !defined(__STDC__)

/*  ---  External Name  ---						    */

/*  Macro to define the name of a Fortran routine or common block. This	    */
/*  ends in an underscore on many Unix systems.				    */

#undef  F77_EXTERNAL_NAME
#define F77_EXTERNAL_NAME(X) X/**/_


/*  ---  Dummy Arguments  ---						    */

/*  Macros to handle character dummy arguments.				    */

#undef  TRAIL
#define TRAIL(X) ,int X/**/_length


/*  ---  Declare variables  ---						    */

#undef  DECLARE_CHARACTER
#define DECLARE_CHARACTER(X,L)         F77_CHARACTER_TYPE X[L]; \
   const int X/**/_length = L
#undef  DECLARE_CHARACTER_ARRAY
#define DECLARE_CHARACTER_ARRAY(X,L,D) F77_CHARACTER_TYPE X[D][L]; \
   const int X/**/_length = L
#undef DECLARE_CHARACTER_DYN
#define DECLARE_CHARACTER_DYN(X)   F77_CHARACTER_TYPE *X;\
   int X/**/_length
#undef DECLARE_CHARACTER_ARRAY_DYN
#define DECLARE_CHARACTER_ARRAY_DYN(X)   F77_CHARACTER_TYPE *X;\
   int X/**/_length
#undef F77_CREATE_CHARACTER
#define F77_CREATE_CHARACTER(X,L)  X=cnfCref(L);\
   X/**/_length = L
#undef F77_CREATE_CHARACTER_ARRAY
#define F77_CREATE_CHARACTER_ARRAY(X,L,N) \
        {int f77dims[1];f77dims[0]=N;X=cnfCrefa(L,1,f77dims);X/**/_length=L;}

/*  ---  Pass arguments to a FORTRAN routine  ---			    */

#undef  TRAIL_ARG
#define TRAIL_ARG(X)     ,X/**/_length


#endif  /* of non ANSI redefinitions					    */


/* -----------------------------------------------------------------------  */

/* The standard macros defined above are known to work with the following   */
/* systems:								    */
 
/*--------
|   Sun   |
---------*/

/*  On SunOS, the ANSI definitions work with the acc and gcc compilers.     */
/*  The cc compiler uses the non ANSI definitions. It also needs the K&R    */
/*  definitions in the file kr.h.					    */
/*  On Solaris, the standard definitions work with the cc compiler.	    */

#if defined(sun)

#if !defined(__STDC__)
#if !defined(_F77_KR)
#define _F77_KR
#endif
#endif

#endif	/* Sun								    */

/* -------------------- System dependent sections ------------------------- */

/*------------
|   VAX/VMS   |
-------------*/

/* Many macros need to be changed due to the way that VMS handles external  */
/* names, passes character arguments and handles logical values.	    */


#if defined(VMS)

/*  ---  Data Types  ---						    */

/*  Redefine the macro for the byte data type as signed is not valid syntax */
/*  as the VMS compiler is not ANSI compliant.				    */

#undef  F77_BYTE_TYPE
#define F77_BYTE_TYPE      char


/*  ---  External Names  ---						    */

/*  Macro to define the name of a Fortran routine or common block.	    */
/*  Fortran and C routines names are the same on VMS.			    */

#undef  F77_EXTERNAL_NAME
#define F77_EXTERNAL_NAME(X) X


/*  ---  Dummy Arguments  ---						    */

/*  Macros to handle character arguments.				    */
/*  Character string arguments are pointers to character string descriptors */
/*  and there are no trailing arguments.				    */

#if( VMS != 0 )
#include <descrip.h>
#endif


#undef  F77_CHARACTER_ARG_TYPE
#define F77_CHARACTER_ARG_TYPE struct dsc$descriptor_s
#undef  F77_CHARACTER_ARRAY_ARG_TYPE
#define F77_CHARACTER_ARRAY_ARG_TYPE struct dsc$descriptor_a
#undef  CHARACTER
#define CHARACTER(X) F77_CHARACTER_ARG_TYPE *CNF_CONST X/**/_arg
#undef  TRAIL
#define TRAIL(X)
#undef  CHARACTER_ARRAY
#define CHARACTER_ARRAY(X)  F77_CHARACTER_ARRAY_ARG_TYPE *CNF_CONST X/**/_arg
#undef  GENPTR_CHARACTER
#define GENPTR_CHARACTER(X) \
   F77_CHARACTER_TYPE *X = X/**/_arg->dsc$a_pointer; \
   int X/**/_length = X/**/_arg->dsc$w_length;
#undef  GENPTR_CHARACTER_ARRAY
#define GENPTR_CHARACTER_ARRAY(X)   GENPTR_CHARACTER(X)


/*  ---  Logical Values  ---						    */

#undef  F77_TRUE
#define F77_TRUE -1
#undef  F77_ISTRUE
#define F77_ISTRUE(X) ( (X)&1 )
#undef  F77_ISFALSE
#define F77_ISFALSE(X) ( ! ( (X)&1 ) )


/*  ---  Common Blocks  ---						    */

#undef  F77_BLANK_COMMON
#define F77_BLANK_COMMON  $BLANK


/*  ---  Declare Variables  ---						    */

#undef  DECLARE_CHARACTER
#define DECLARE_CHARACTER(X,L) \
   F77_CHARACTER_TYPE X[L];    const int X/**/_length = L; \
   F77_CHARACTER_ARG_TYPE X/**/_descr = \
      { L, DSC$K_DTYPE_T, DSC$K_CLASS_S, X }; \
   F77_CHARACTER_ARG_TYPE *X/**/_arg = &X/**/_descr
#undef  DECLARE_CHARACTER_ARRAY
#define DECLARE_CHARACTER_ARRAY(X,L,D) \
   F77_CHARACTER_TYPE X[D][L]; const int X/**/_length = L; \
   F77_CHARACTER_ARRAY_ARG_TYPE X/**/_descr = \
      { L, DSC$K_DTYPE_T, DSC$K_CLASS_S, X }; \
   F77_CHARACTER_ARRAY_ARG_TYPE *X/**/_arg = &X/**/_descr


/*  ---  The dynamic allocation of character arguments  ---                 */
#undef DECLARE_CHARACTER_DYN
#define DECLARE_CHARACTER_DYN(X) int X/**/_length;\
                                  F77_CHARACTER_ARG_TYPE *X/**/_arg;\
                                  F77_CHARACTER_TYPE *X
#undef DECLARE_CHARACTER_ARRAY_DYN
#define DECLARE_CHARACTER_ARRAY_DYN(X) int X/**/_length;\
                                  F77_CHARACTER_ARRAY_ARG_TYPE *X/**/_arg;\
                                  F77_CHARACTER_TYPE *X
#undef F77_CREATE_CHARACTER
#define F77_CREATE_CHARACTER(X,L) X/**/_arg = cnfCref(L);\
                                  X = X/**/_arg->dsc$a_pointer; \
                                  X/**/_length = X/**/_arg->dsc$w_length
#undef F77_CREATE_CHARACTER_ARRAY
#define F77_CREATE_CHARACTER_ARRAY(X,L,N) \
  {int f77dims[1];f77dims[0]=N;X/**/_arg=cnfCrefa(L,1,f77dims);X/**/_length=L;}
#define F77_CREATE_CHARACTER_ARRAY_M(X,L,N,D) X/**/_arg = cnfCrefa(L,N,D);\
                                  X = X/**/_arg->dsc$a_pointer; \
                                  X/**/_length = X/**/_arg->dsc$w_length
#undef F77_FREE_CHARACTER
#define F77_FREE_CHARACTER(X)     cnfFreef( X/**/_arg )

/*  ---  Pass arguments to a FORTRAN routine  ---			    */

#undef  CHARACTER_ARG
#define CHARACTER_ARG(X) X/**/_arg
#undef  CHARACTER_ARRAY_ARG
#define CHARACTER_ARRAY_ARG(X) X/**/_arg
#undef  TRAIL_ARG
#define TRAIL_ARG(X)

#endif  /* VMS								    */

/* -----------------------------------------------------------------------  */

/*--------------------------
|   DECstation Ultrix (cc)  |
|   DECstation Ultrix (c89) |
|   DECstation OSF/1        |
|   Alpha OSF/1             |
 --------------------------*/

/* Do this complicated set of definitions as a single #if cannot be	    */
/* continued across multiple lines.					    */

#if defined(mips) && defined(ultrix)
#define _dec_unix 1
#endif
#if defined(__mips) && defined(__ultrix)
#define _dec_unix 1
#endif
#if defined(__mips__) && defined(__osf__)
#define _dec_unix 1
#endif
#if defined(__alpha) && defined(__osf__)
#define _dec_unix 1
#endif

#if _dec_unix

/*  The macros for Ultrix are the same as the standard ones except for ones */
/*  dealing with logical values. The ANSI definitions work with the c89	    */
/*  compiler, and the non ANSI definitions work with the cc compiler.	    */
/*  The same applies to DEC OSF/1, except that its cc compiler is ANSI	    */
/*  compliant.								    */


/*  ---  Logical Values  ---						    */

/*  Redefine macros that evaluate to a C logical value, given a FORTRAN	    */
/*  logical value. These definitions are only valid when used with the DEC  */
/*  FORTRAN for RISC compiler. If you are using the earlier FORTRAN for	    */
/*  RISC compiler from MIPS, then these macros should be deleted.	    */

#undef  F77_TRUE
#define F77_TRUE -1
#undef  F77_ISTRUE
#define F77_ISTRUE(X) ( (X)&1 )
#undef  F77_ISFALSE
#define F77_ISFALSE(X) ( ! ( (X)&1 ) )


#endif  /* DEC Unix							    */

/*
*+
*  Name:
*     cnf.h

*  Purpose:
*     Function prototypes for cnf routines

*  Language:
*     ANSI C

*  Type of Module:
*     C include file

*  Description:
*     These are the prototype definitions for the functions in the CNF
*     library. They are used used in mixing C and FORTRAN programs.

*  Copyright:
*     Copyright (C) 1991 Science & Engineering Research Council

*  Authors:
*     PMA: Peter Allan (Starlink, RAL)
*     AJC: Alan Chipperfield (Starlink, RAL)
*     {enter_new_authors_here}

*  History:
*     23-MAY-1991 (PMA):
*        Original version.
*     12-JAN-1996 (AJC):
*        Add cnf_cref and cnf_freef
*     14-JUN-1996 (AJC):
*        Add cnf_crefa, imprta, exprta
*                crela, impla, expla
*     18-JUL-1996 (AJC):
*        Add impch and expch
*     17-MAR-1998 (AJC):
*        Add imprtap and exprtap
*     {enter_changes_here}

*  Bugs:
*     {note_any_bugs_here}

*-
------------------------------------------------------------------------------
*/
void *cnfCalloc( size_t, size_t );
void cnfCopyf( const char *source_f, int source_len, char *dest_f,
                int dest_len );
void *cnfCptr( F77_POINTER_TYPE );
char *cnfCreat( int length );
F77_CHARACTER_ARG_TYPE *cnfCref( int length );
F77_CHARACTER_ARG_TYPE *cnfCrefa( int length, int ndims, const int *dims );
char *cnfCreib( const char *source_f, int source_len );
char *cnfCreim( const char *source_f, int source_len );
F77_LOGICAL_TYPE *cnfCrela( int ndims, const int *dims );
void cnfExpch( const char *source_c, char *dest_f, int nchars );
void cnfExpla( const int *source_c, F77_LOGICAL_TYPE *dest_f, int ndims,
                const int *dims );
void cnfExpn( const char *source_c, int max, char *dest_f, int dest_len );
void cnfExprt( const char *source_c, char *dest_f, int dest_len );
void cnfExprta( const char *source_c, int source_len, char *dest_f,
                 int dest_len, int ndims, const int *dims );
void cnfExprtap( char *const *source_c, char *dest_f, int dest_len, 
                  int ndims, const int *dims );
F77_POINTER_TYPE cnfFptr( void *cpointer );
void cnfFree( void * );
void cnfFreef( F77_CHARACTER_ARG_TYPE *temp );
void cnfImpb( const char *source_f, int source_len, char *dest_c );
void cnfImpbn( const char *source_f, int source_len, int max, char *dest_c );
void cnfImpch( const char *source_f, int nchars, char *dest_c );
void cnfImpla( const F77_LOGICAL_TYPE *source_f, int *dest_c,
                int ndims, const int *dims );
void cnfImpn( const char *source_f, int source_len, int max, char *dest_c );
void cnfImprt( const char *source_f, int source_len, char *dest_c );
void cnfImprta( const char *source_f, int source_len, char *dest_c,
                 int dest_len, int ndims, const int *dims );
void cnfImprtap( const char *source_f, int source_len, char *const *dest_c,
                  int dest_len, int ndims, const int *dims );
int cnfLenc( const char *source_c );
int cnfLenf( const char *source_f, int source_len );
void *cnfMalloc( size_t );
int cnfRegp( void * );
void cnfUregp( void * );
#endif

#ifndef CNF_OLD_DEFINED
#define CNF_OLD_DEFINED
/* Define old names to be new names */
#define cnf_calloc cnfCalloc
#define cnf_copyf cnfCopyf
#define cnf_cptr cnfCptr
#define cnf_creat cnfCreat
#define cnf_cref cnfCref
#define cnf_crefa cnfCrefa
#define cnf_creib cnfCreib
#define cnf_creim cnfCreim
#define cnf_crela cnfCrela
#define cnf_expch cnfExpch
#define cnf_expla cnfExpla
#define cnf_expn cnfExpn
#define cnf_exprt cnfExprt
#define cnf_exprta cnfExprta
#define cnf_exprtap cnfExprtap
#define cnf_fptr cnfFptr
#define cnf_free cnfFree
#define cnf_freef cnfFreef
#define cnf_impb cnfImpb
#define cnf_impbn cnfImpbn
#define cnf_impch cnfImpch
#define cnf_impla cnfImpla
#define cnf_impn cnfImpn
#define cnf_imprt cnfImprt
#define cnf_imprta cnfImprta
#define cnf_imprtap cnfImprtap
#define cnf_lenc cnfLenc
#define cnf_lenf cnfLenf
#define cnf_malloc cnfMalloc
#define cnf_regp cnfRegp
#define cnf_uregp cnfUregp

#endif
