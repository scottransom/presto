#include <stdlib.h>
#include <stdio.h>
#include "Python.h"
#include "arrayobject.h"

/* Here is an example library function that returns an array (a 1D
   vector).  This is just representative -- the actual functions in my
   library are much more complicated.  They are all malloc'd,
   though. */

double * my_range(long n)
{
  double *v;
  long i;

  v = (double *) malloc((size_t) (sizeof(double) * n));
  if (!v) {
    printf("\nAllocation error in my_range()\n");
    exit(1);
  }
  for (i=0; i<n; i++) v[i] = (double) i;
  return v;
}

/* Here is a representative wrapper function that I am using to
   interface the above routine with Python (using Numeric arrays). */

static PyObject *wrap_my_range(PyObject *self, PyObject *args)
{
  PyObject *obj;
  PyArrayObject *arr;
  double *result;
  long n;
  
  self = self;
  if(!PyArg_ParseTuple(args,"O:wrap_my_range",&obj)) 
    return NULL;
  n = PyInt_AsLong((PyObject *)obj);

  result = my_range(n);

  arr = (PyArrayObject *)PyArray_FromDims(1, (int *)&n, PyArray_DOUBLE);
  if (arr == NULL) return NULL;

  /* Have to copy the data, until I figure out a better way... */

  memcpy(arr->data, (char *)result, n*arr->descr->elsize);

  /* Free the malloced array from the 'C' routine... */

  free(result);

  PyArray_INCREF(arr);
  return (PyObject *)arr;
}

/* Here are the module initialization functions */

static PyObject *ErrorObject;
static PyMethodDef my_range_method[] = {
  { "my_range", wrap_my_range, 1 },
  { NULL, NULL }
};

void initmy_range()
{
  PyObject *m, *d;

  m = Py_InitModule("my_range", my_range_method);
  d = PyModule_GetDict(m);
  import_array();
  ErrorObject = PyString_FromString("my_range.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module my_range");
}
