/* Note:  Define rawtype to be the type you need to transpose */

#ifndef _RAWTYPE_PART_DECLARED_
/* Use the following if you want to transform doubles */
/* typedef double rawtype_part; */
typedef float rawtype_part;
#define _RAWTYPE_PART_DECLARED_
#endif				/* _RAWTYPE_PART_DECLARED_ */

/* Declare a floating point complex type */
#ifndef _FCOMPLEX_DECLARED_
typedef struct FCOMPLEX {
  rawtype_part r, i;
} fcomplex;
#define _FCOMPLEX_DECLARED_
#endif				/* _FCOMPLEX_DECLARED_ */

/* Declare rawtype as fcomplex */
#ifndef _RAWTYPE_DECLARED_
typedef fcomplex rawtype;
typedef rawtype TOMS_el_type;
#define MPI_rawtype fcomplextype
#define _RAWTYPE_DECLARED_
#endif				/* _RAWTYPE_DECLARED_ */


