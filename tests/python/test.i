%module my_range

%include typemaps.i

%{

#include <stdlib.h>
#include <stdio.h>
#include "arrayobject.h"


/*
Use the following macro to show the properties of a PyArrayObject.
Call it with something like the following:

#ifdef DEBUG_NUMPY
  PRINT_PYARRAY(arr)
#endif

Where arr is a pointer to the PyArrayObject you want to examine.
*/

/* #define DEBUG_NUMPY 1 */

#define PRINT_PYARRAY(obj) \
{ \
    int i; \
    printf("      ob_refcnt = %d\n", (obj)->ob_refcnt); \
    printf("             nd = %d\n", (obj)->nd); \
    for (i=0; i<(obj)->nd; i++){ \
      printf("  dimensions[%d] = %d\n", i, (obj)->dimensions[i]); \
      printf("     strides[%d] = %d\n", i, (obj)->strides[i]); \
    } \
    printf("          flags = %d  (", (obj)->flags); \
    if ((obj)->flags & CONTIGUOUS) printf(" CONTIGUOUS "); \
    if ((obj)->flags & OWN_DIMENSIONS) printf(" OWN_DIMENSIONS "); \
    if ((obj)->flags & OWN_STRIDES) printf(" OWN_STRIDES "); \
    if ((obj)->flags & OWN_DATA) printf(" OWN_DATA "); \
    printf(")\n"); \
    printf("  descr.typenum = %d\n", (obj)->descr->type_num); \
    printf("   desrc.elsize = %d\n", (obj)->descr->elsize); \
    printf("      descr.one = %s\n", (obj)->descr->one); \
    printf("     descr.zero = %s\n", (obj)->descr->zero); \
    printf(" descr.chartype = %c\n", (obj)->descr->type); \
}


typedef struct {
   double real;
   double imag;
} double_complex;

typedef struct {
   float real;
   float imag;
} float_complex;

double * my_range(long n);

%}


%init %{
  import_array();
%}

%wrapper %{
static long _output_arraylen = 0;
%}

// Return an array where each value is a double
%typemap(python, out) double* {
  PyArrayObject *arr;
  int n, tempn = 1;
  
  n = _output_arraylen;
  arr = (PyArrayObject *)PyArray_FromDims(1, (int *)&tempn, PyArray_DOUBLE);
  if (arr == NULL) return NULL;
  free(arr->data);
  arr->dimensions[0] = n;
  arr->data = (char *)$source;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}


// Identify the variable that will give the length of the array
//  that will be returned.  It should be an int or a long.
%typemap(python, in) long ARRAYLEN, int ARRAYLEN {
  _output_arraylen = PyInt_AsLong((PyObject *)$source);
  $target = _output_arraylen;
}


%apply long ARRAYLEN { long n }
double * my_range(long n);
/* Here is an example library function that returns an array (a 1D
   vector).  This is just representative -- the actual functions in my
   library are much more complicated.  They are all malloc'd,
   though. */
