//
//        Rough NumPy Typemap for SWIG 
//   Scott M. Ransom (ransom@cfa.harvard.edu)
// 

%include typemaps.i

%{
#include "arrayobject.h"

#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), (m)->descr->type_num, 0,0)))

#ifndef _FCOMPLEX_DECLARED_
typedef struct FCOMPLEX {
    float r, i;
} fcomplex;
#define _FCOMPLEX_DECLARED_
#endif				/* _FCOMPLEX_DECLARED_ */

#ifndef _DCOMPLEX_DECLARED_
typedef struct DCOMPLEX {
    double r, i;
} dcomplex;
#define _DCOMPLEX_DECLARED_
#endif				/* _DCOMPLEX_DECLARED_ */

%}


%init %{
  import_array();
%}

/* The following are used when the C routine outputs a Numeric Array */
/* The variables below describe the sizes of the NumPy arrays.       */
%wrapper %{
static long _output_arraylen = 0;
static long _output_matrixrows = 0;
static long _output_matrixcols = 0;
%}


// Identify the variable that will give the length of the array
//  that will be returned.  It should be an int or a long.

%typemap(python, in) long ARRAYLEN, int ARRAYLEN {
  _output_arraylen = PyInt_AsLong((PyObject *)$source);
  $target = _output_arraylen;
}
%typemap(python, in) long MATRIXROWS, int MATRIXROWS {
  _output_matrixrows = PyInt_AsLong((PyObject *)$source);
  $target = _output_matrixrows;
}
%typemap(python, in) long MATRIXCOLS, int MATRIXCOLS {
  _output_matrixcols = PyInt_AsLong((PyObject *)$source);
  $target = _output_matrixcols;
}

// Type mapping for grabbing a FILE * from Python

%typemap(python,in) FILE * {
  if (!PyFile_Check($source)) {
    PyErr_SetString(PyExc_TypeError, "Need a file!");
    return NULL;
  }		
  $target = PyFile_AsFile($source);
}

// The following are used to convert various input
// 1-D PyArrays into vectors to be used by the C routine.

%typemap(python, in) float* IN_1D_FLOAT {
  PyArrayObject *arr;
  
  /* Check that obj is really a 1D array of bytes */
  
  if (!PyArray_Check($source)) {
    PyErr_SetString(PyExc_TypeError, \
		    "First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */

  if (PyArray_ObjectType($source,0) != PyArray_FLOAT) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of FLOAT");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$source);
  if (arr->nd != 1) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 1D array");
    return NULL;
  }
  $target = (float *)arr->data;
  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}


%typemap(python, in) double* IN_1D_DOUBLE {
  PyArrayObject *arr;
  
  /* Check that obj is really a 1D array of bytes */
  
  if (!PyArray_Check($source)) {
    PyErr_SetString(PyExc_TypeError,"First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */
  
  if (PyArray_ObjectType($source,0) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of DOUBLE");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$source);
  if (arr->nd != 1) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 1D array");
    return NULL;
  }
  $target = (double *)arr->data;
  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}


%typemap(python, in) fcomplex* IN_1D_CFLOAT {
  PyArrayObject *arr;
  
  /* Check that obj is really a 1D array of bytes */
  
  if (!PyArray_Check($source)) {
    PyErr_SetString(PyExc_TypeError,"First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */
  
  if (PyArray_ObjectType($source,0) != PyArray_CFLOAT) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of CFLOAT");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$source);
  if (arr->nd != 1) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 1D array");
    return NULL;
  }
  $target = (fcomplex *)arr->data;
  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}


%typemap(python, in) dcomplex* IN_1D_CDOUBLE {
  PyArrayObject *arr;
  
  /* Check that obj is really a 1D array of bytes */
  
  if (!PyArray_Check($source)) {
    PyErr_SetString(PyExc_TypeError,"First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */
  
  if (PyArray_ObjectType($source,0) != PyArray_CDOUBLE) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of CDOUBLE");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$source);
  if (arr->nd != 1) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 1D array");
    return NULL;
  }
  $target = (dcomplex *)arr->data;
  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}


%typemap(python, in) float* IN_2D_FLOAT {
  PyObject *obj;
  PyArrayObject *arr;
  int i, j;
  unsigned char *ptr;
  
  /* Check that obj is really an 2D array of bytes */
  
  if (!PyArray_Check($source)) {
    PyErr_SetString(PyExc_TypeError, \
		    "First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */

  if (PyArray_ObjectType($source,0) != PyArray_FLOAT) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of FLOAT");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$source);
  if (arr->nd != 2) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 2d array");
    return NULL;
  }
  $target = (float *)arr->data;
  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}


%typemap(python, in) float* IN_2D_DOUBLE {
  PyObject *obj;
  PyArrayObject *arr;
  int i, j;
  unsigned char *ptr;
  
  /* Check that obj is really an 2D array of bytes */
  
  if (!PyArray_Check($source)) {
    PyErr_SetString(PyExc_TypeError, \
		    "First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */
  
  if (PyArray_ObjectType($source,0) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of DOUBLE");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$source);
  if (arr->nd != 2) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 2d array");
    return NULL;
  }
  $target = (double *)arr->data;
  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}

//  Add typemaps for complex output values

%typemap(python,ignore) fcomplex *L_OUTPUT(fcomplex temp),
                        dcomplex *L_OUTPUT(dcomplex temp)
{
  $target = &temp;
}

%typemap(python,ignore) fcomplex *T_OUTPUT(fcomplex temp),
                        dcomplex *T_OUTPUT(dcomplex temp)
{
  $target = &temp;
}

%typemap(python,argout) fcomplex *L_OUTPUT,
                        dcomplex *L_OUTPUT
{
  PyObject *o;
  double rl, im;

  rl = (double) $source->r;
  im = (double) $source->i;
  o = PyComplex_FromDoubles(rl, im);
  $target = l_output_helper($target, o);
}

%typemap(python,argout) fcomplex *T_OUTPUT,
                        dcomplex *T_OUTPUT
{
  PyObject *o;
  double rl, im;

  rl = (double) $source->r;
  im = (double) $source->i;
  o = PyComplex_FromDoubles(rl, im);
  $target = t_output_helper($target, o);
}

#ifdef OUTPUT_LIST
%typemap(python,ignore) fcomplex  *OUTPUT = fcomplex *L_OUTPUT;
%typemap(python,ignore) dcomplex  *OUTPUT = dcomplex *L_OUTPUT;

%typemap(python,argout) fcomplex  *OUTPUT = fcomplex *L_OUTPUT;
%typemap(python,argout) dcomplex  *OUTPUT = dcomplex *L_OUTPUT;
#else
%typemap(python,ignore) fcomplex  *OUTPUT = fcomplex *T_OUTPUT;
%typemap(python,ignore) dcomplex  *OUTPUT = dcomplex *T_OUTPUT;

%typemap(python,argout) fcomplex  *OUTPUT = fcomplex *T_OUTPUT;
%typemap(python,argout) dcomplex  *OUTPUT = dcomplex *T_OUTPUT;
#endif



//  Functions to let us retrieve arguments as list output.

%typemap(python,argout) double *OUTDOUBLE {
  PyObject *o;
  o = PyFloat_FromDouble(*$source);
  if ((!$target) || ($target == Py_None)) {
    $target = o;
  } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    PyList_Append($target,o);
    Py_XDECREF(o);
  }
}

//  Functions to let us retrieve arguments as list output.

%typemap(python,argout) float *OUTFLOAT {
  PyObject *o;
  o = PyFloat_FromFloat(*$source);
  if ((!$target) || ($target == Py_None)) {
    $target = o;
  } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    PyList_Append($target,o);
    Py_XDECREF(o);
  }
}

//  Functions to let us retrieve arguments as list output.

%typemap(python,argout) long *OUTLONG {
  PyObject *o;
  o = PyInt_FromLong(*$source);
  if ((!$target) || ($target == Py_None)) {
    $target = o;
  } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    PyList_Append($target,o);
    Py_XDECREF(o);
  }
}

//  Functions to let us retrieve arguments as list output.

%typemap(python,argout) int *OUTINT {
  PyObject *o;
  o = PyInt_FromInt(*$source);
  if ((!$target) || ($target == Py_None)) {
    $target = o;
  } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    PyList_Append($target,o);
    Py_XDECREF(o);
  }
}

//  Functions to help us return PyArrays from a C function

// Return an array where each value is a float
%typemap(python, out) float * {
  PyArrayObject *arr;
  int n;
  
  n = _output_arraylen;
  _output_arraylen = 0;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(1, (int *)&n, PyArray_FLOAT, (char *)$source);
  if (arr == NULL) return NULL;
  arr->flags |= OWN_DATA;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}


// Return a complex array where each value is a double
%typemap(python, out) double * {
  PyArrayObject *arr;
  int n;
  
  n = _output_arraylen;
  _output_arraylen = 0;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(1, (int *)&n, PyArray_DOUBLE, (char *)$source);
  if (arr == NULL) return NULL;
  arr->flags |= OWN_DATA;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}


// Return a complex array where each value is a float
%typemap(python, out) fcomplex * {
  PyArrayObject *arr;
  int n;
  
  n = _output_arraylen;
  _output_arraylen = 0;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(1, (int *)&n, PyArray_CFLOAT, (char *)$source);
  if (arr == NULL) return NULL;
  arr->flags |= OWN_DATA;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}


// Return a complex matrix where each value is a float
%typemap(python, out) fcomplex ** {
  PyArrayObject *arr;
  int n[2];

  n[0] = _output_matrixcols;
  n[1] = _output_matrixrows;
  _output_matrixrows = 0;
  _output_matrixcols = 0;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(2, n, PyArray_CFLOAT, (char *)$source[0]);
  free($source);
  if (arr == NULL) return NULL;
  arr->flags |= OWN_DATA;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}


// Return a complex array where each value is a double
%typemap(python, out) dcomplex * {
  PyArrayObject *arr;
  int n;
  
  n = _output_arraylen;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(1, (int *)&n, PyArray_CDOUBLE, (char *)$source);
  if (arr == NULL) return NULL;
  arr->flags |= OWN_DATA;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}








