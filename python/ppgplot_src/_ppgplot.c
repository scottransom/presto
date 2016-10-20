/*
 * FILE:
 *    _ppgplot.c
 * DESCRIPTION:
 *    ppgplot: Python / Numeric Python interface to the PGPLOT graphical 
 *    library.
 *    Tested with PGPLOT 5.2, Python 2.3.1 and Numeric 23.1, on
 *    Linux (v2.4.x).
 * AUTHOR(S):
 *    Nick Patavalis (npat@efault.net)
 * NOTES:
 *    - A few ppgplot functions have not been interfaced yet.
 *    - The pythonic calling conventions of some functions are *not*
 *      identical to the original PGPLOT ones.
 */

#include <Python.h>

#include <assert.h>

#include <cpgplot.h>
#include "numpy/arrayobject.h"

/************************************************************************/

/*
 * Debuging Switches
 */
/*
  #define DEBUG_CONT_S
  #define DEBUG_PGCIR
  #define DEBUG_TOARRAY
*/
/************************************************************************/

/*
 * Default values and stuff
 */

#define STD_DEVICE "/XWINDOW"
#define DEF_XLABEL "x"
#define DEF_YLABEL "y"
#define DEF_PLOTLABEL "x = f(y)"

/************************************************************************/

/*
 * Handy shortcuts.
 */

#define PYF(_fn) static PyObject *(_fn)(PyObject *self, PyObject *args)
#define PYRN Py_INCREF(Py_None); return(Py_None)

/************************************************************************/

/*
 * Globals
 */

static PyObject *PpgIOErr;
static PyObject *PpgTYPEErr;
static PyObject *PpgMEMErr;

/**************************************************************************/
/*                          support functions                             */
/**************************************************************************/

static PyObject *
tofloatvector (PyObject *o, float **v, int *vsz)
{
    PyArrayObject *a1, *af1, *af2;
    PyArray_Descr *descr;
    npy_intp dims;
    int ownedaf1=0;

    /* Check if args are arrays. */
    if (!PyArray_Check(o)) {
		PyErr_SetString(PpgTYPEErr,"object is not an array");
		return(NULL);
    }
    a1 = (PyArrayObject *)o;
    /* Check if args are vectors. */
    if (a1->nd != 1) {
		PyErr_SetString(PpgTYPEErr,"object is not a vector");
		return(NULL);
    }
    
#ifdef DEBUG_TOARRAY
    fprintf(stderr,"(tofloatvector): array type = %d\n",a1->descr->type_num);
#endif

    switch (a1->descr->type_num) {
    case PyArray_FLOAT:
		af1 = a1;
		break;
    case PyArray_CHAR:
#ifndef USE_NUMARRAY 
    case PyArray_UBYTE: 
#endif
//    case PyArray_SBYTE:
    case PyArray_SHORT:
    case PyArray_INT:
#ifndef USE_NUMARRAY
    case PyArray_LONG:
#endif
    case PyArray_DOUBLE:
		if (!(af1 = (PyArrayObject *)PyArray_Cast(a1,PyArray_FLOAT))) {
			PyErr_SetString(PpgTYPEErr,"cannot cast vector to floats");
			return(NULL);
		}
		ownedaf1 = 1;
		break;
    default:
		PyErr_SetString(PpgTYPEErr,"cannot cast vector to floats");
		return(NULL);
		break;
    }
    
#ifdef DEBUG_TOARRAY
    fprintf(stderr,"(tofloatvector): array type = %d\n",a1->descr->type_num);
#endif
    
    af2 = af1;
    descr = PyArray_DescrFromType(PyArray_FLOAT); 
    if (PyArray_AsCArray((PyObject **)&af2, (void *)v, &dims, 1,
                         descr) == -1) {
		af2 = NULL;
    }
    *vsz = dims;

    if (ownedaf1) { Py_DECREF(af1); }

    return((PyObject *)af2);
}

/*************************************************************************/

static PyObject *
tofloatmat(PyObject *o, float **m, int *nr, int *nc)
{
    PyArrayObject *a1, *af1, *af2;
    PyArray_Descr *descr;
    npy_intp dims[2];
    int ownedaf1=0;
    char **tmpdat;
    
    /* Check if args are arrays. */
    if (!PyArray_Check(o)) {
		PyErr_SetString(PpgTYPEErr,"object is not and array");
		return(NULL);
    }
    a1 = (PyArrayObject *)o;
    /* Check if args are matrices. */
    if (a1->nd != 2) {
		PyErr_SetString(PpgTYPEErr,"object is not a matrix");
		return(NULL);
    }
    
#ifdef DEBUG_TOARRAY
    fprintf(stderr,"(tofloatmat): array type = %d\n",a1->descr->type_num);
#endif
    
    switch (a1->descr->type_num) {
    case PyArray_FLOAT:
		af1 = a1;
		break;
    case PyArray_CHAR: 
#ifndef USE_NUMARRAY
    case PyArray_UBYTE: 
#endif
//    case PyArray_SBYTE:
    case PyArray_SHORT: 
    case PyArray_INT: 
#ifndef USE_NUMARRAY
    case PyArray_LONG:
#endif
    case PyArray_DOUBLE:
		if (!(af1 = (PyArrayObject *)PyArray_Cast(a1,PyArray_FLOAT))) {
			PyErr_SetString(PpgTYPEErr,"cannot cast matrix to floats");
			return(NULL);
		}
		ownedaf1 = 1;
		break;
    default:
		PyErr_SetString(PpgTYPEErr,"cannot cast matrix to floats");
		return(NULL);
		break;
    }
    
#ifdef DEBUG_TOARRAY
    fprintf(stderr,"(tofloatmat): array type = %d\n",a1->descr->type_num);
#endif
    
    af2 = af1;
    descr = PyArray_DescrFromType(PyArray_FLOAT); 
    if (PyArray_AsCArray((PyObject **)&af2, (void *)&tmpdat, dims, 2,
                         descr) == -1) {
		af2 = NULL;
		goto bailout;
    }
    *nr = dims[0];
    *nc = dims[1];
    
    /* WARNING: What follows is a little tricky and I dunno if I'm 
       really allowed to do this. On the other hand it really conserves 
       time and memory! So this assert statement will make sure that 
       the program *will* blow in your face if what I'm doing here 
       turns-out be bogus. */
    assert((af2->dimensions[1] * af2->descr->elsize) == af2->strides[0]);
    
    /* Phew! we 're clear! */
    *m = (float *)(*tmpdat);
    /* tmpdat was malloc'ed inside PyArray_As2D. We must free it.
       Look at the code of PyArray_As2D for details... */
    free(tmpdat);

bailout:
    if (ownedaf1) { Py_DECREF(af1); }
    return((PyObject *)af2);
}

/**************************************************************************/

#ifdef DEBUG_TOARRAY

PYF(tstmat)
{
    PyObject *o=NULL;
    PyArrayObject *af=NULL;
    float *v;
    int i=0,j=0, nc=0, nr=0;
    
    if(!PyArg_ParseTuple(args,"O",&o)) return(NULL);
    
    if (!(af =(PyArrayObject *)tofloatmat(o,&v,&nr,&nc))) goto fail;
    
    for (i=0; i<nr; i++) {
		for(j=0; j<nc; j++)
			fprintf(stderr,"%f\t",v[i*nc+j]);
		fprintf(stderr,"\n");
    }

    Py_DECREF(af);
    PYRN;

fail:
    if (af) Py_DECREF(af);
    return(NULL);
}

#endif

/**************************************************************************/

void
minmax (float *v, int nsz, float *min, float *max)
{
    register float *e;
    register float mn, mx;

    for (mn=mx=*v, e=v+nsz; v < e; v++)
		if (*v > mx) mx = *v;
		else if (*v < mn) mn = *v;
    *min = mn;
    *max = mx;
}


void 
lininterp (float min, float max, int npts, float *v)
{
    register int i;
    register float step;
    register float lev;

    step = (max-min) / (npts-1);
    lev = min;
    for (i=0; i<npts; i++) {
		v[i] = lev;
		lev += step;
    }
}


static void
autocal2d(float *a, int rn, int cn,
		  float *fg, float *bg, int nlevels, float *levels,
		  float *x1, float *x2, float *y1, float *y2,
		  float *tr)
{
    float dx1, dx2, dy1, dy2;

    /* autocalibrate intensity-range. */
    if (*fg == *bg) {
		minmax(a,rn*cn,bg,fg);
/* 	fprintf(stderr,"Intensity range:\n  fg=%f\n  bg=%f\n",*fg,*bg); */
    }
    
    if ((nlevels >= 2) && (levels))
		lininterp(*bg, *fg, nlevels, levels);
    
    /* autocalibrate x-y range. */
    if ((*x1 == *x2) || (*y1 == *y2)) cpgqwin(&dx1,&dx2,&dy1,&dy2);
    if (*x1 == *x2) {*x1=dx1; *x2=dx2;}
    if (*y1 == *y2) {*y1=dy1; *y2=dy2;}
/*     fprintf(stderr,"Xrange: [%f, %f]\nYrange[%f, %f]\n",*x1,*x2,*y1,*y2); */

    /* calculate transformation vector. */
    tr[2] = tr[4] = 0.0;
    tr[1] = (*x2 - *x1) / cn;
    tr[0] = *x1 - (tr[1] / 2);
    tr[5] = (*y2 - *y1) / rn;
    tr[3] = *y1 - (tr[5] / 2);

/*     fprintf(stderr,"Tansformation vector:\n"); */
/*     for (i=0; i<6; fprintf(stderr,"  tr[%d]=%f\n",i,tr[i]),i++); */
}


/* ###### ADDED ###### */
/*
  pgqinf(item : string)
  return string
*/
PYF(pgqinf)
{
	char *item;
	char value[512];
	int length;
	memset(value, 0, 511);
	value[511] = 0;
	//printf("blaat\n");

	if(!PyArg_ParseTuple(args, "s", &item))
		return NULL;
	//printf(item);
	length = 512;
	cpgqinf(item, value, &length);
	//printf("length = %i, value = %s\n", length, value);
	return(Py_BuildValue("s", value));
}

/*
  pgpoly(x, y : array)
*/
PYF(pgpoly)
{
	int n;
	int xlen;
	int ylen;
	float* xpoints = NULL;
	float* ypoints = NULL;
	PyObject *xarray = NULL;
	PyObject *yarray = NULL;
	PyArrayObject *xa = NULL;
	PyArrayObject *ya = NULL;

	if(!PyArg_ParseTuple(args, "OO", &xarray, &yarray))
		return NULL;
	if(!(xa = (PyArrayObject*)tofloatvector(xarray, &xpoints, &xlen))) goto fail;
	if(!(ya = (PyArrayObject*)tofloatvector(yarray, &ypoints, &ylen))) goto fail;
	/*for(int i = 0; i < xlen; i++)
	  {
	  printf("%i: %f,%f\n", i, xpoints
	  }*/
	n = xlen;
	if(ylen < n)
		n = ylen;
	cpgpoly(n, xpoints, ypoints);

	Py_DECREF(xa);
	Py_DECREF(ya);

	PYRN;

fail:
	if(xa) { Py_DECREF(xa); }
	if(ya) { Py_DECREF(ya); }
	return NULL;
}

/*
  pgqvp(units : int)
  return x1, x2, y1, y2
*/
PYF(pgqvp)
{
	int units;
	float x1, x2, y1, y2;

	if(!PyArg_ParseTuple(args, "i", &units))
		return NULL;
	cpgqvp(units, &x1, &x2, &y1, &y2);

	return Py_BuildValue("ffff", x1, x2, y1, y2);
}

/*
  pgqvsz(units : int)
  return x1, x2, y1, y2
*/
PYF(pgqvsz)
{
	int units;
	float x1, x2, y1, y2;

	if(!PyArg_ParseTuple(args, "i", &units))
		return NULL;
	cpgqvsz(units, &x1, &x2, &y1, &y2);
	return Py_BuildValue("ffff", x1, x2, y1, y2);
}

PYF(pgsclp)
{
	int state;
	if(!PyArg_ParseTuple(args, "i", &state))
		return NULL;
	cpgsclp(state);
	PYRN;
}

PYF(pgqclp)
{
	int state;
	if(!PyArg_ParseTuple(args, ""))
		return NULL;
	cpgqclp(&state);
	return Py_BuildValue("i", state);
}

PYF(pgqtxt)
{
	float x, y, angle, fjust;
	char* text;
	float b1[4], b2[4];
	PyObject* t[4];
	PyObject* ret;
	int i;

	if(!PyArg_ParseTuple(args, "ffffs", &x, &y, &angle, &fjust, &text))
		return NULL;
	cpgqtxt(x, y, angle, fjust, text, b1, b2);
	for(i = 0; i < 4; i++)
		t[i] = Py_BuildValue("ff", b1[i], b2[i]);
	ret = Py_BuildValue("OOOO", t[0], t[1], t[2], t[3]);
	for(i = 0; i < 4; i++)
		Py_DECREF(t[i]);
	return ret;
}

PYF(pgconf)
{
    PyObject
		*oa=NULL, *otr=NULL;
    PyArrayObject
		*aa=NULL, *atr=NULL;
    float *a=NULL, *tr=NULL, c_1 = 0.0, c_2 = 0.0;
    int cd=0 ,rd=0,
		c1=0,c2=0,r1=0,r2=0,
		trsz=0;

    if (!PyArg_ParseTuple(args,"OiiiiiiffO:pgconl",
						  &oa, &cd, &rd, &c1, &c2, &r1, &r2,
						  &c_1, &c_2, &otr))
		return(NULL);

    if (!(aa = (PyArrayObject *)tofloatmat(oa, &a, &rd, &cd)))
		goto fail;
    if (!(atr = (PyArrayObject *)tofloatvector(otr, &tr, &trsz)))
		goto fail;

    if (trsz < 6) {
		PyErr_SetString(PpgTYPEErr,"contour: invalid transform. vector");
		goto fail;
    }

    cpgconf(a,cd,rd,c1+1,c2+1,r1+1,r2+1,c_1, c_2,tr);

    Py_DECREF(aa);
    Py_DECREF(atr);
    PYRN;

fail:
    if (aa) { Py_DECREF(aa); }
    if (atr) { Py_DECREF(atr); }
    return(NULL);
}


/* #### END ADDED #####*/

/***************************************************************************/
/*                         management functions                            */
/***************************************************************************/

PYF(pgbeg)
{
    char *device = STD_DEVICE;
    int xnsub=1, ynsub=1;

    if (!PyArg_ParseTuple(args, "|sii:pgbeg", &device, &xnsub, &ynsub))
		return NULL;

    if ( cpgbeg(0,device,xnsub,ynsub) != 1 ) {
		PyErr_SetString(PpgIOErr, "Failed to open plot device.");
		return(NULL);
    }
    
    PYRN;
}


PYF(pgopen)
{
    char *dev = NULL, did = 0;

    if (!PyArg_ParseTuple(args,"|z:pgopen",&dev))
		return(NULL);

    if (!dev) dev = STD_DEVICE;

    if ( (did = cpgopen(dev)) <= 0 ) {
		PyErr_SetString(PpgIOErr, "Failed to open plot device.");
		return(NULL);
    }

    return(Py_BuildValue("i",did));
}    

PYF(pgslct)
{
    int did = 0;

    if (!PyArg_ParseTuple(args,"|i:pgslct",&did))
		return(NULL);

    cpgslct(did);

    PYRN;
}


PYF(pgask)
{
    int i;

    if (!PyArg_ParseTuple(args,"i:pgask",&i)) return(NULL);
    
    cpgask(i);

    PYRN;
}

PYF(pgbbuf)
{
    cpgbbuf();

    PYRN;
}

PYF(pgebuf)
{
    cpgebuf();

    PYRN;
}

PYF(pgclos)
{
    cpgclos();

    PYRN;
}

PYF(pgend)
{
    cpgend();

    PYRN;
}

PYF(pgldev)
{
    cpgldev();

    PYRN;
}

PYF(pgpage)
{
    cpgpage();

    PYRN;
}

PYF(pgpanl)
{
    int x = 1, y = 1;

    if (!PyArg_ParseTuple(args,"|ii:pgpanl",&x,&y))
		return(NULL);

    cpgpanl(x,y);

    PYRN;
}

PYF(pgsubp)
{
    int xnsub = 1, ynsub = 1;

    if (!PyArg_ParseTuple(args,"|ii:pgsubp",&xnsub, &ynsub))
		return(NULL);

    cpgsubp(xnsub, ynsub);

    PYRN;
}

PYF(pgupdt)
{
    cpgupdt();

    PYRN;
}

/**************************************************************************/
/*                           Environment functions.                       */
/**************************************************************************/

/*
 * Coordinate systems, Windows, Viewports and stuff
 */

PYF(pgenv)
{
    float xn,xx,yn,yx;
    int just=0, type=0;

    if (!PyArg_ParseTuple(args,"ffff|ii:pgenv",&xn,&xx,&yn,&yx,&just,&type))
		return(NULL);

    cpgenv(xn,xx,yn,yx,just,type);

    PYRN;
}

PYF(pgswin)
{
    float x1 = 0.0, x2 = 1.0, y1 = 0.0, y2 = 1.0;
    
    if (!PyArg_ParseTuple(args,"|ffff:pgswin", 
						  &x1, &x2, &y1, &y2))
		return(NULL);
    
    cpgswin(x1,x2,y1,y2);
    
    PYRN;
}

PYF(pgvstd)
{
    cpgvstd();

    PYRN;
}


PYF(pgsvp)
{
    float xleft = 0.0, xright = 1.0, ybot = 0.0, ytop = 1.0;
    
    if (!PyArg_ParseTuple(args,"|ffff:pgsvp", 
						  &xleft, &xright, &ybot, &ytop))
		return(NULL);
    
    cpgsvp(xleft,xright,ybot,ytop);
    
    PYRN;
}

PYF(pgvsiz)
{
    float xleft = 0.0, xright = 0.0, ybot = 0.0, ytop = 0.0;
    
    if (!PyArg_ParseTuple(args,"ffff:pgvsiz", 
						  &xleft, &xright, &ybot, &ytop))
		return(NULL);
    
    cpgvsiz(xleft,xright,ybot,ytop);
    
    PYRN;
}

PYF(pgwnad)
{
    float x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
    
    if (!PyArg_ParseTuple(args,"ffff:pgwand", 
						  &x1, &x2, &y1, &y2))
		return(NULL);
    
    cpgwnad(x1,x2,y1,y2);

    PYRN;
}

PYF(pgbox)
{
    char *xopt=NULL, *yopt=NULL;
    float xtick=0.0, ytick=0.0;
    int nxsub=0, nysub=0;

    if (!PyArg_ParseTuple(args,"|zfizfi:pgbox",
						  &xopt, &xtick, &nxsub,
						  &yopt, &ytick, &nysub))
		return(NULL);
    
    if (!xopt) xopt = "ABCGNTS";
    if (!yopt) yopt = "ABCGNTS";

    cpgbox(xopt,xtick,nxsub,yopt,ytick,nysub);

    PYRN;
}

PYF(pgtbox)
{
    char *xopt=NULL, *yopt=NULL;
    float xtick=0.0, ytick=0.0;
    int nxsub=0, nysub=0;

    if (!PyArg_ParseTuple(args,"|zfizfi:pgtbox",
						  &xopt, &xtick, &nxsub,
						  &yopt, &ytick, &nysub))
		return(NULL);

    if (!xopt) xopt = "ABCGNTSYXH";
    if (!yopt) yopt = "ABCGNTS";

    cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub);

    PYRN;
}

PYF(pgpap)
{
    float width, aspect = 1;

    if (!PyArg_ParseTuple(args,"f|f:pgpap",&width, &aspect))
		return(NULL);

    cpgpap(width, aspect);

    PYRN;
}

/*
 * Drawing and misc. parametrs
 */

PYF(pgsci)
{
    int ci = 1;

    if (!PyArg_ParseTuple(args,"|i:pgsci", &ci)) return(NULL);

    cpgsci(ci);
    
    PYRN;
}

PYF(pgsfs)
{
    int fs = 1;
    
    if (!PyArg_ParseTuple(args,"|i:pgsfs", &fs)) return(NULL);

    cpgsfs(fs);
    
    PYRN;
}

PYF(pgsitf)
{
    int itf = 1;
    
    if (!PyArg_ParseTuple(args,"|i:pgsitf", &itf)) return(NULL);

    cpgsitf(itf);
    
    PYRN;
}

PYF(pgstbg)
{
    int tbg = 0;
    
    if (!PyArg_ParseTuple(args,"|i:pgstbg", &tbg)) return(NULL);

    cpgstbg(tbg);
    
    PYRN;
}

PYF(pgshs)
{
    float angle = 30.0, sepn = 1.0, phase = 0.0;
    
    if (!PyArg_ParseTuple(args,"|fff:pgshs", 
						  &angle, &sepn, &phase))
		return(NULL);
    
    cpgshs(angle, sepn, phase);
    
    PYRN;
}

PYF(pgsls)
{
    int ls = 0;

    if (!PyArg_ParseTuple(args,"|i:pgsls",&ls))
		return(NULL);

    cpgsls(ls);

    PYRN;
}

PYF(pgscf)
{
    int cf = 1;

    if (!PyArg_ParseTuple(args,"|i:pgscf",&cf))
		return(NULL);

    cpgscf(cf);

    PYRN;
}

PYF(pgsch)
{
    float ch = 1.0;

    if (!PyArg_ParseTuple(args,"|f:pgsch",&ch))
		return(NULL);

    cpgsch(ch);

    PYRN;
}


PYF(pgsah)
{
    int fs = 1;
    float angle=45.0, vent=0.3;

    if (!PyArg_ParseTuple(args,"|iff:pgsah",&fs, &angle, &vent))
		return(NULL);
    
    cpgsah(fs, angle, vent);

    PYRN;
}

PYF(pgslw)
{
    int lw = 1;

    if (!PyArg_ParseTuple(args,"|i:pgslw",&lw))
		return(NULL);

    cpgslw(lw);

    PYRN;
}

/*
 * Colors
 */

PYF(pgscir)
{
    int cir1=0, cir2=0;
	
    if (!PyArg_ParseTuple(args,"|ii:pgscir", &cir1, &cir2))
		return(NULL);
	
    if ((cir1 == 0) && (cir2 == 0)) {
		cpgqcol(&cir1, &cir2);
#ifdef DEBUG_PGCIR
		fprintf(stderr,"Device: lowcol = %d, hicol = %d\n",cir1,cir2);
#endif
		if (cir2 > 15)
			cir1 = 16;
		else
			cir1 = 1;
    }
	
#ifdef DEBUG_PGCIR
    fprintf(stderr,"2D Maps: lowcol = %d, hicol = %d\n",cir1,cir2);
#endif

    cpgscir(cir1,cir2);
	
    PYRN;
}

PYF(pgctab)
{
    PyObject *lo = NULL, *ro = NULL, *go = NULL, *bo = NULL;
    PyArrayObject *la = NULL, *ra = NULL, *ga = NULL, *ba = NULL;
    float *l = NULL, *r = NULL, *g = NULL, *b = NULL,
		contra = 1.0, bright = 0.5;
    int nc = 0, nr=0, ng=0, nb=0;

    if (!PyArg_ParseTuple(args,"OOOO|iff:pgctab",
						  &lo, &ro, &go, &bo, &nc, &contra, &bright))
		return(NULL);

    if (!(la =(PyArrayObject *)tofloatvector(lo, &l, &nc))) goto fail;
    if (!(ra =(PyArrayObject *)tofloatvector(ro, &r, &nr))) goto fail;
    if (!(ga =(PyArrayObject *)tofloatvector(go, &g, &ng))) goto fail;
    if (!(ba =(PyArrayObject *)tofloatvector(bo, &b, &nb))) goto fail;

    if ((nr < nc) || (ng < nc) || (nb < nc)) {
		PyErr_SetString(PpgTYPEErr,"pgtab: invalid color tables");
		goto fail;
    }

    cpgctab(l,r,g,b,nc,contra,bright);

    Py_DECREF(la);
    Py_DECREF(ra);
    Py_DECREF(ga);
    Py_DECREF(ba);
    PYRN;

fail:
    if (la) { Py_DECREF(la); }
    if (ra) { Py_DECREF(ra); }
    if (ga) { Py_DECREF(ga); }
    if (ba) { Py_DECREF(ba); }
    return(NULL);
}

PYF(pgscr)
{
    float cr = 0.0, cg = 1.0, cb = 0.0;
    int ci = 0;
    
    if (!PyArg_ParseTuple(args,"ifff:pgscr", 
						  &ci, &cr, &cg, &cb))
		return(NULL);
    
    cpgscr(ci,cr,cg,cb);
    
    PYRN;
}

PYF(pgshls)
{
    float ch = 0.0, cl = 1.0, cs = 0.0;
    int ci = 0;
    
    if (!PyArg_ParseTuple(args,"ifff:pgshls", 
						  &ci, &ch, &cl, &cs))
		return(NULL);
    
    cpgshls(ci, ch, cl, cs);
    
    PYRN;
}


/*
 * Save and retrieve env.
 */

PYF(pgsave)
{
    cpgsave();

    PYRN;
}

PYF(pgunsa)
{
    cpgunsa();

    PYRN;
}

/************************************************************************/
/*                             Query functions.                         */
/************************************************************************/

/*
 * Interactive Stuff
 */

PYF(pgcurs)
{
    float x=0.0, y=0.0;
    char ch = '\0';

    if (!PyArg_ParseTuple(args,"|ff:pgcurs",&x, &y))
		return(NULL);

    cpgcurs(&x,&y,&ch);

    return(Py_BuildValue("ffc",x,y,ch));
}


PYF(pgband)
{
    int mode=7, i=0;
    float xref = 0.0, yref = 0.0, x=0.0, y=0.0;
    char ch = '\0';

    if (!PyArg_ParseTuple(args,"i|iff:pgband",
						  &mode, &i, &xref, &yref))
		return(NULL);


    cpgband(mode,i,xref,yref,&x,&y,&ch);

    return(Py_BuildValue("ffc",x,y,ch));
}

PYF(pgqcol)
{
    int ci1, ci2;

    cpgqcol(&ci1, &ci2);

    return(Py_BuildValue("ii",ci1, ci2));
}

PYF(pgqcir)
{
    int ci1, ci2;

    cpgqcir(&ci1, &ci2);

    return(Py_BuildValue("ii",ci1, ci2));
}

PYF(pgqcr)
{
    int ci=0;
    float cr, cg, cb;

    if (!PyArg_ParseTuple(args,"|i",&ci))
		return(NULL);

    cpgqcr(ci, &cr,&cg,&cb);

    return(Py_BuildValue("fff",cr,cg,cb));
}


PYF(pgqah)
{
    int fs=0;
    float angle=0.0, vent=0.0;

    cpgqah(&fs, &angle, &vent);

    return(Py_BuildValue("iff", fs, angle, vent));
}

PYF(pgqhs)
{
    float angle=0.0, sepn=0.0, phase=0.0;

    cpgqhs(&angle, &sepn, &phase);

    return(Py_BuildValue("fff", angle, sepn, phase));
}

PYF(pgqpos)
{
    float x=0.0, y=0.0;

    cpgqpos(&x, &y);

    return(Py_BuildValue("ff", x, y));
}

PYF(pgqcs)
{
    int units=4;
    float xch, ych;

    if (!PyArg_ParseTuple(args,"|i",&units))
		return(NULL);

    cpgqcs(units, &xch, &ych);

    return(Py_BuildValue("ff",xch, ych));
}

PYF(pgqitf)
{
    int itf;

    cpgqitf(&itf);

    return(Py_BuildValue("i",itf));
}

PYF(pgqls)
{
    int ls;

    cpgqls(&ls);

    return(Py_BuildValue("i",ls));
}

PYF(pgqlw)
{
    int lw;

    cpgqlw(&lw);

    return(Py_BuildValue("i",lw));
}


PYF(pgqcf)
{
    int cf;

    cpgqcf(&cf);

    return(Py_BuildValue("i",cf));
}

PYF(pgqch)
{
    float ch;

    cpgqch(&ch);

    return(Py_BuildValue("f",ch));
}

PYF(pgqci)
{
    int ci;

    cpgqci(&ci);

    return(Py_BuildValue("i",ci));
}


PYF(pgqfs)
{
    int fs;

    cpgqfs(&fs);

    return(Py_BuildValue("i",fs));
}

PYF(pgqid)
{
    int did;

    cpgqid(&did);

    return(Py_BuildValue("i",did));
}

PYF(pglen)
{
    char *text;
    int units=4;
    float xltxt, yltxt;

    if (!PyArg_ParseTuple(args,"s|i:pglen", &text, &units))
		return(NULL);

    cpglen(units, text, &xltxt, &yltxt);

    return(Py_BuildValue("ff",xltxt, yltxt));
}


/************************************************************************/
/*                       Anotation functions.                           */
/************************************************************************/

PYF(pgiden)
{
    cpgiden();

    PYRN;
}

PYF(pglab)
{
    char *xl=DEF_XLABEL, *yl=DEF_YLABEL, *pl=DEF_PLOTLABEL;

    if (!PyArg_ParseTuple(args,"|sss:pglab",&xl,&yl,&pl))
		return(NULL);

    cpglab(xl,yl,pl);

    PYRN;
}

PYF(pgmtxt)
{
    char
		*side = "T",
		*text;
    float
		disp = 1,
		coord = 0.5,
		fjust = 0.5;

    if (!PyArg_ParseTuple(args,"sfffs:pglab",&side,&disp,&coord,&fjust,
						  &text))
		return(NULL);
    
    cpgmtxt(side, disp, coord, fjust, text);

    PYRN;
}

PYF(pgtext)
{
    float x = 0.0, y=0.0;
    char *string = "ppgplot";
    
    if (!PyArg_ParseTuple(args,"ffs:pgtext", &x, &y, &string)) return(NULL);

    cpgtext(x,y,string);

    PYRN;
}

PYF(pgptxt)
{
    float x = 0.0, y=0.0, angle = 0.0, fjust = 0.0;
    char *string = "ppgplot";
    
    if (!PyArg_ParseTuple(args,"ffffs:pgptxt", &x, &y, 
						  &angle, &fjust, &string)) 
		return(NULL);
	
    cpgptxt(x, y, angle, fjust, string);
    
    PYRN;
}


/************************************************************************/
/*                      Simple Drawing                                  */
/************************************************************************/

PYF(pgeras)
{
    cpgeras();

    PYRN;
}

PYF(pgetxt)
{
    cpgetxt();

    PYRN;
}

PYF(pgarro)
{
    float x1,y1,x2,y2;

    if (!PyArg_ParseTuple(args,"ffff:pgarro",&x1,&y1,&x2,&y2))
		return(NULL);

    cpgarro(x1,y1,x2,y2);

    PYRN;
}

PYF(pgcirc)
{
    float x, y, r;

    if (!PyArg_ParseTuple(args,"fff:pgcirc",&x,&y,&r)) return(NULL);
    
    cpgcirc(x,y,r);

    PYRN;
}

PYF(pgmove)
{
    float x, y;

    if (!PyArg_ParseTuple(args,"ff:pgmove",&x,&y)) return(NULL);
    
    cpgmove(x,y);

    PYRN;
}


PYF(pgdraw)
{
    float x, y;

    if (!PyArg_ParseTuple(args,"ff:pgdraw",&x,&y)) return(NULL);
    
    cpgdraw(x,y);

    PYRN;
}

PYF(pgrect)
{
    float x1, x2, y1, y2;

    if (!PyArg_ParseTuple(args,"ffff:pgrect",&x1, &x2, &y1, &y2)) 
		return(NULL);
    
    cpgrect(x1,x2,y1,y2);

    PYRN;
}


PYF(pgline)
{
    int xsz, ysz;
    PyObject *op1=NULL, *op2=NULL, *o1=NULL, *o2=NULL;
    float *xpts, *ypts;

    if(!PyArg_ParseTuple(args, "OO", &op1, &op2)) goto fail;
    if (!(o1 = tofloatvector(op1, &xpts, &xsz))) goto fail;
    if (!(o2 = tofloatvector(op2, &ypts, &ysz))) goto fail;

    if (xsz > ysz) xsz = ysz;
    cpgline(xsz,xpts,ypts);

    Py_DECREF(o1);
    Py_DECREF(o2);
    PYRN;

fail:
    if (o1) { Py_DECREF(o1); }
    if (o2) { Py_DECREF(o2); }
    return(NULL);
}

PYF(pgpt)
{
    int xsz, ysz;
    PyObject *op1=NULL, *op2=NULL, *o1=NULL, *o2=NULL;
    float *xpts, *ypts;
    int symbol=0;

    if(!PyArg_ParseTuple(args, "OO|i", &op1, &op2, &symbol)) goto fail;
    if (!(o1 = tofloatvector(op1, &xpts, &xsz))) goto fail;
    if (!(o2 = tofloatvector(op2, &ypts, &ysz))) goto fail;

    if (xsz > ysz) xsz = ysz;
    cpgpt(xsz,xpts,ypts, symbol);

    Py_DECREF(o1);
    Py_DECREF(o2);
    PYRN;

fail:
    if (o1) { Py_DECREF(o1); }
    if (o2) { Py_DECREF(o2); }
    return(NULL);
}



/************************************************************************/
/*                           2D maps.                                   */
/************************************************************************/


/*
 * ----------------------------------------------------------
 * NOTE: Array syntax and storage in python (C) and Fortran:
 * ----------------------------------------------------------
 *
 * Memory is visualized as a notebook, where you write left to right, and
 * when the line (row) ends, you go to the leftmost part of the next
 * one and start again, and so on. This paradigm is used throughout
 * the comments in this file, and throughout all the associated 
 * documentation. If you like to visualize memory as writting top
 * to bottom, then just replace the word "column" with the word 
 * "row" (and vice-versa) everywhere.
 *
 * in python (and C) a matrix looks like: a[i,j] i=rows j=cols. 
 * In mem it is: 
 *              a[1,0],    a[1,1],    a[1,2]   ....a[1,nc-1]
 *              c[2,0],    a[2,1],    a[2,2]   ....a[2,nc-1]
 *		       ..................
 *              a[nr-1,1], a[nr-1,1], a[nr-1,2]....a[nr-1,nc-1]
 * A.k.a the major (faster running) index is written (spelled) last.
 *
 * The same matrix in fortran syntax looks like: a(i,j) i=cols 
 * j=rows. 
 * In mem it is:
 *              a(1,1),    a(2,1),    a(3,1),  ....a(nc,1)
 *              a(1,2),    a(2,2),    a(3,2),  ....a(nc,2)
 *		       ..................
 *              a(1,nr),   a(2,nr),   a(3,nr), ....a(nc,nr)
 * A.k.a the major (faster running) index is written (spelled) first.
 *
 * In this file (and all associated documentation) whenever i'm using
 * the () as indexing operators, I'm talking fortran-lingo, and whenever 
 * I'm using the [] as an indexing operator, I'm talking C/Python-lingo.
 *
 * In the pgplot documentation the FORTRAN notation is always used.
 * So translating the function prototype found in the PGPLOT docs, 
 * to the variable naming sceme used in the following function, yields 
 * something like this:
 *
 * idim ---> cd (columns dimension / # of columns)
 * jdim ---> rd (rows dimension / # of rows)
 * i1   ---> c1 (first column)
 * i2   ---> c2 (last column)
 * j1   ---> r1 (first row)
 * j2   ---> r2 (last row)
 *
 * I repeat that columns run faster.
 * 
 * Also remember than when you use python (even if calling PGPLOT)
 * all indeces start from 0, *not* 1.
 * 
 * ----------------------------------------------------------
 * NOTE: A word about translations:
 * ----------------------------------------------------------
 *
 * In order to plot (map) a matrix 
 *            a[r1:r2,c1:c2] 
 * in a viewport scaled as:
 *            x-axis form x1 to x2
 *            y-axis from y1 to y2
 * having rows run on x and columns on run y, you have to use the
 * following transformation matrix:
 *
 * tr = array([t0,t1,t2,t3,t4,t5])
 * where
 *        t0 = x1 - (t1 * c1) / 2
 *        t1 = (x2-x1) / (c2-c1)
 *        t2 = 0.0
 *        t3 = y1 - (t1 * r1) / 2
 *        t4 = (y2 - y1) / (c2 - c1)
 *        t5 = 0.0
 *
 * The 2 denominator in t0 and t3 is justified by the fact that
 * cordinates pin-point the CENTER of the pixels, and not their
 * SW corner.
 * t2 and t5 are used only if you want shearing or rotations.
 *
 * ---------------------------------------------------------
 * NOTE:
 * ---------------------------------------------------------
 *
 * If all these sound just too-much for you, you should 
 * consider using pggray_s which is a simplified allternative 
 * to pggray. Not the complete but enough of pgray's functionality
 * is wrapped around an easy-to-use interface. The basic limitation
 * of pggray_s is that it cannot do rotations and shearing
 * (it can do flips though!)
 * The same is true for color image-maps (pgimag).
 *
 */

static PyObject *
ImageMap(int color, PyObject *args)
{
    float fg=0.0, bg=0.0, *a=NULL, *tr=NULL;
    int cd=0, rd=0, c1=0, c2=0, r1=0, r2=0, rn=0, cn=0, sz=0;
    PyObject *oa = NULL, *ot = NULL;
    PyArrayObject *aa = NULL, *at = NULL;


    if (!PyArg_ParseTuple(args,"OiiiiiiffO:pggray",
						  &oa, &cd, &rd, &c1, &c2, &r1, &r2,
						  &fg, &bg, &ot))
		return(NULL);

    if (!(aa =(PyArrayObject *)tofloatmat(oa, &a, &rn, &cn))) goto fail;
    if (!(at =(PyArrayObject *)tofloatvector(ot, &tr, &sz))) goto fail;

    if (sz < 6) {
		PyErr_SetString(PpgTYPEErr,"pggray: invalid transform. vactor");
		goto fail;
    }

    if (color)
		cpgimag(a, cn, rn, c1+1, c2+1, r1+1, r2+1, bg, fg, tr);
    else
		cpggray(a, cn, rn, c1+1, c2+1, r1+1, r2+1, fg, bg, tr);

    Py_DECREF(aa);
    Py_DECREF(at);
    PYRN;

fail:
    if (aa) { Py_DECREF(aa); }
    if (at) { Py_DECREF(at); }
    return(NULL);
}

PyObject *
ImageMap_s (int color, PyObject *args)
{
    PyObject *oa=NULL;
    PyArrayObject *aa=NULL;
    float 
		*a=NULL, 
        tr[6], levels[10],
		x1=0.0, y1=0.0, x2=0.0, y2=0.0, 
		fg=0.0, bg=0.0;
    int 
		rn=0, cn=0;
		    
    if (!PyArg_ParseTuple(args,"O|ffffff:imagemap_s",
						  &oa, &fg, &bg, &x1, &y1, &x2, &y2))
		return(NULL);
    
    if (!(aa =(PyArrayObject *)tofloatmat(oa, &a, &rn, &cn))) return(NULL);

    /* Perform autocalibrations as nesecairy. */
    autocal2d(a, rn, cn, &fg, &bg, 5, levels, &x1, &x2, &y1, &y2, tr);

    if (color)
		cpgimag(a, cn, rn, 0+1, cn, 0+1, rn, bg, fg, tr);
    else
		cpggray(a, cn, rn, 0+1, cn, 0+1, rn, fg, bg, tr);
    
    Py_DECREF(aa);
    PYRN;
}

PYF(pggray)
{
    return(ImageMap(0,args));
}

PYF(pggray_s)
{
    return(ImageMap_s(0,args));
}

PYF(pgimag)
{
    return(ImageMap(1,args));
}

PYF(pgimag_s)
{
    return(ImageMap_s(1,args));
}

PYF(pgwedg)
{
    char *side = NULL, *label = NULL;
    float disp=0.0, width=0.0, fg=0.0, bg=0.0;

    if (!PyArg_ParseTuple(args,"sffff|s:pgwdg",
						  &side, &disp, &width, &fg, &bg, &label))
		return(NULL);

    if (!label) label = " ";
    
    cpgwedg(side, disp, width, fg, bg, label);

    PYRN;
}

PYF(pgwedg_s)
{
    char *side = NULL, *label = NULL;
    float disp=1.0, width=4.0, fg=0.0, bg=0.0;

    if (!PyArg_ParseTuple(args,"ff|zzff:pgwdg",
						  &fg, &bg, &side, &label, &disp, &width))
		return(NULL);

    if (!side) side = "RG";
    if (!label) label = " ";

    cpgwedg(side, disp, width, fg, bg, label);

    PYRN;
}

/*********************************************************************/
/*                       Histrograms.                                */
/*********************************************************************/

PYF(pgerrb)
{
    PyObject *ox=NULL, *oy=NULL, *oe=NULL;
    PyArrayObject *ax=NULL, *ay=NULL, *ae=NULL;
    float *x=NULL, *y=NULL, *e=NULL, t=1.0;
    int szx=0, szy=0, sze =0, dir = 0, n1;

    if (!PyArg_ParseTuple(args,"iOOO|f:pgerrb", &dir, &ox, &oy, &oe, &t))
		return(NULL);

    if (!(ax = (PyArrayObject *)tofloatvector(ox, &x, &szx))) goto fail;
    if (!(ay = (PyArrayObject *)tofloatvector(oy, &y, &szy))) goto fail;
    if (!(ae = (PyArrayObject *)tofloatvector(oe, &e, &sze))) goto fail;

    /* this is n1 = min(szx, szy, sze) */
    n1=(szx<szy)?((szx<sze)?szx:sze):((szy<sze)?szy:sze);

    cpgerrb(dir, n1, x, y, e, t);

    Py_DECREF(ax);
    Py_DECREF(ay);
    Py_DECREF(ae);
    PYRN;

fail:
    if (ax) { Py_DECREF(ax); }
    if (ay) { Py_DECREF(ay); }
    if (ae) { Py_DECREF(ae); }
    return(NULL);
}

PYF(pgerrx)
{
    PyObject *oy=NULL, *ox1=NULL, *ox2=NULL;
    PyArrayObject *ay=NULL, *ax1=NULL, *ax2=NULL;
    float *y=NULL, *x1=NULL, *x2=NULL, t=1.0;
    int szy=0, szx1=0, szx2 =0, n1;

    if (!PyArg_ParseTuple(args,"OOO|f:pgerrx", &ox1, &ox2, &oy, &t))
        return(NULL);

    if (!(ay = (PyArrayObject *)tofloatvector(oy, &y, &szy))) goto fail;
    if (!(ax1 = (PyArrayObject *)tofloatvector(ox1, &x1, &szx1))) goto fail;
    if (!(ax2 = (PyArrayObject *)tofloatvector(ox2, &x2, &szx2))) goto fail;

    /* this is n1 = min(szx, szx1, szx2) */
    n1=(szy<szx1)?((szy<szx2)?szy:szx2):((szx1<szx2)?szx1:szx2);

    cpgerrx(n1, x1, x2, y, t);

    Py_DECREF(ay);
    Py_DECREF(ax1);
    Py_DECREF(ax2);
    PYRN;

fail:
    if (ay) { Py_DECREF(ay); }
    if (ax1) { Py_DECREF(ax1); }
    if (ax2) { Py_DECREF(ax2); }
    return(NULL);
}


PYF(pgerry)
{
    PyObject *ox=NULL, *oy1=NULL, *oy2=NULL;
    PyArrayObject *ax=NULL, *ay1=NULL, *ay2=NULL;
    float *x=NULL, *y1=NULL, *y2=NULL, t=1.0;
    int szx=0, szy1=0, szy2 =0, n1;

    if (!PyArg_ParseTuple(args,"OOO|f:pgerry", &ox, &oy1, &oy2, &t))
        return(NULL);

    if (!(ax = (PyArrayObject *)tofloatvector(ox, &x, &szx))) goto fail;
    if (!(ay1 = (PyArrayObject *)tofloatvector(oy1, &y1, &szy1))) goto fail;
    if (!(ay2 = (PyArrayObject *)tofloatvector(oy2, &y2, &szy2))) goto fail;

    /* this is n1 = min(szx, szy1, szy2) */
    n1=(szx<szy1)?((szx<szy2)?szx:szy2):((szy1<szy2)?szy1:szy2);

    cpgerry(n1, x, y1, y2, t);

    Py_DECREF(ax);
    Py_DECREF(ay1);
    Py_DECREF(ay2);
    PYRN;

fail:
    if (ax) { Py_DECREF(ax); }
    if (ay1) { Py_DECREF(ay1); }
    if (ay2) { Py_DECREF(ay2); }
    return(NULL);
}




PYF(pghist)
{
    PyObject *od=NULL;
    PyArrayObject *ad=NULL;
    int 
		nusrpts = 0,
		npts = 0,
		nbin = 0,
		pgflag = 0;
    float
		datamin=0,
		datamax=0,
		*data=NULL;
    
    if (!PyArg_ParseTuple(args,"iOffi|i:pghist", &nusrpts, &od, &datamin, 
						  &datamax, &nbin, &pgflag))
		return(NULL);
    
    if (!(ad = (PyArrayObject *)tofloatvector(od,&data,&npts))) goto fail;

    if (nusrpts > npts) nusrpts = npts;
    
    cpghist(nusrpts, data, datamin, datamax, nbin, pgflag);
    
    Py_DECREF(ad);
    PYRN;

fail:
    if (ad) { Py_DECREF(ad); }
    return(NULL);
}

PYF(pghist_s)
{
    PyObject *od=NULL;
    PyArrayObject *ad=NULL;
    int 
		npts = 0,
		nbin = 0,
		pgflag = 0;
    float
		datamin=0,
		datamax=0,
		*data=NULL;
    
    if (!PyArg_ParseTuple(args,"Oi|iff:pghist_s", &od, &nbin, &pgflag, 
						  &datamin, &datamax))
		return(NULL);
    
    if (!(ad = (PyArrayObject *)tofloatvector(od,&data,&npts))) goto fail;

    if (datamin == datamax)
		minmax(data, npts, &datamin, &datamax);

    cpghist(npts, data, datamin, datamax, nbin, pgflag);

    Py_DECREF(ad);
    PYRN;

fail:
    if (ad) { Py_DECREF(ad); }
    return(NULL);
}


PYF(pgbin)
{
    PyObject *ox=NULL, *od=NULL;
    PyArrayObject *ax=NULL, *ad=NULL;
    int 
		szx=0, szd=0,
		center = 1;
    float
		*x,
		*data;

    if (!PyArg_ParseTuple(args,"OO|i:pgbin", &ox, &od, &center))
		return(NULL);

    if (!(ax = (PyArrayObject *)tofloatvector(ox,&x,&szx))) goto fail;
    if (!(ad = (PyArrayObject *)tofloatvector(od,&data,&szd))) goto fail;

    if (szx > szd) szx = szd;

    cpgbin(szx,x,data,center);

    Py_DECREF(ax);
    Py_DECREF(ad);
    PYRN;

fail:
    if (ax) { Py_DECREF(ax); }
    if (ad) { Py_DECREF(ad); }
    return(NULL);
}

PYF(pgbin_s)
{
    PyObject  *od=NULL;
    PyArrayObject *ad=NULL;
    int nbin = 0, center = 1;
    float x1 = 0.0, x2 = 0.0, *x, *data,
		dummy1, dummy2;

    if (!PyArg_ParseTuple(args,"O|ffi:pgbin_s",
						  &od, &x1, &x2, &center))
		return(NULL);

    if (!(ad = (PyArrayObject *)tofloatvector(od,&data,&nbin)))
		return(NULL);

/*     fprintf(stderr,"x1=%f\nx2=%f\n",x1,x2); */

    if (!(x = (float *)malloc(nbin * sizeof(*x)))) {
		PyErr_SetString(PpgMEMErr,"pgbin_s: out of memory!");
		Py_DECREF(ad);
		return(NULL);
    }

    if (x1==x2) cpgqwin(&x1,&x2,&dummy1, &dummy2);
    lininterp(x1,x2,nbin,x);

/*     fprintf(stderr,"bin=%d\n",nbin); */
/*     fprintf(stderr,"x1=%f\nx2=%f\n",x1,x2); */
/*     for (i=0; i<nbin; i++) */
/* 	fprintf(stderr,"  x[%d]=%f\n",i,x[i]); */
    
    cpgbin(nbin,x,data,center);

    free(x);
    Py_DECREF(ad);
    PYRN;
}


PYF(pghi2d)
{
    PyObject 
		*od=NULL, *ox=NULL, *oyl=NULL;
    PyArrayObject 
		*ad=NULL, *ax=NULL, *ayl=NULL;
    int 
		cd=0, rd=0, c1=0, c2=0, r1=0, r2=0,
		ioff, center, vxsz=0, vylsz=0;
    float *md = NULL, *vx=NULL, *vyl,
		bias=0;

    if (!PyArg_ParseTuple(args,"OiiiiiiOifiO:pghi2d", 
						  &od,&cd,&rd,&c1,&c2,&r1,&r2,&ox,
						  &ioff, &bias, &center, &oyl))
		return(NULL);

    if (!(ad = (PyArrayObject *)tofloatmat(od,&md, &rd, &cd))) goto fail;
    if (!(ax = (PyArrayObject *)tofloatvector(ox,&vx, &vxsz))) goto fail;
    if (!(ayl = (PyArrayObject *)tofloatvector(oyl,&vyl, &vylsz))) goto fail;

    if ((vxsz != vylsz) || (vxsz != cd)) {
		PyErr_SetString(PpgTYPEErr,"pghi2d: it must be: " 
						"x size == y-lims size == data-columns");
		goto fail;
    }

    cpghi2d(md, cd, rd, c1+1, c2+1, r1+1, r2+2, 
			vx, ioff, bias, center, vyl);

    Py_DECREF(ad);
    Py_DECREF(ax);
    Py_DECREF(ayl);
    PYRN;
    
fail:
    if (ad) { Py_DECREF(ad); }
    if (ax) { Py_DECREF(ax); }
    if (ayl) { Py_DECREF(ayl); }
    return(NULL);
}


/*
 * This is not good! More work needs to be done.
 * Atomatic bias / offset selection gives awfull results,
 * you have to do it manually. 
 * The best way would be to heve
 * it *resize* the graph (axis) so that the histogram allways
 * fits in the page. Also a way to automatically calculate
 * a reasonable bias would be nice.
 */ 
PYF(pghi2d_s)
{
    PyObject *od = NULL;
    PyArrayObject *ad = NULL;
    float x1=0.0, x2=0.0, 
		dx1=0.0, dx2=0.0, dy1=0.0, dy2=0.0,
		miny=0.0, maxy=0.0,
		bias=0.0,
		*x=NULL, *yl=NULL, *md=NULL;
    int rd=0, cd=0, ioff=1, center = 1; 


    if (!PyArg_ParseTuple(args,"Off|ifi",
						  &od, &x1, &x2, &ioff, &bias, &center)) 
		return(NULL);

    if (!(ad = (PyArrayObject *)tofloatmat(od,&md, &rd, &cd))) 
		goto fail;

    if (!(x = (float *)malloc(cd * sizeof(*x)))) {
		PyErr_SetString(PpgMEMErr,"pghi2d: Out of memory!");
		goto fail;
    }

    if (!(yl = (float *)malloc(cd * sizeof(*yl)))) {
		PyErr_SetString(PpgMEMErr,"pghi2d: Out of memory!");
		goto fail;
    }

    if (bias == 0) {
		cpgqwin(&dx1,&dx2,&dy1,&dy2);
		minmax(md, cd*rd, &miny, &maxy);
		bias = ((dy2 - maxy) / rd) * 0.8;
    }

    lininterp(x1, x2, cd, x);

    cpghi2d(md,cd,rd,0+1,cd,0+1,rd,x,ioff,bias,center,yl);

    Py_DECREF(ad);
    PYRN;
	
fail:
    if (x) free(x);
    if (yl) free(yl);
    if (ad) { Py_DECREF(ad); }
    return(NULL);
}

/************************************************************************/
/*                            Contours.                                 */
/************************************************************************/

enum pp_contour_funcs {FUN_PGCONB, FUN_PGCONS, FUN_PGCONT};

static PyObject *
genContours (enum pp_contour_funcs ft, PyObject *args)
{
    PyObject
		*oa=NULL, *oc=NULL, *otr=NULL;
    PyArrayObject 
		*aa=NULL, *ac=NULL, *atr=NULL;
    float *a=NULL, *c=NULL, *tr=NULL, blank=0.0;
    int cd=0 ,rd=0,
		c1=0,c2=0,r1=0,r2=0, 
		csz=0, trsz=0,
		nc=0;

    if (!PyArg_ParseTuple(args,"OiiiiiiOiO|f:contour",
						  &oa, &cd, &rd, &c1, &c2, &r1, &r2,
						  &oc, &nc, &otr, &blank))
		return(NULL);
    
    if (!(aa = (PyArrayObject *)tofloatmat(oa, &a, &rd, &cd)))
		goto fail;
    if (!(ac = (PyArrayObject *)tofloatvector(oc, &c, &csz)))
		goto fail;
    if (!(atr = (PyArrayObject *)tofloatvector(otr, &tr, &trsz)))
		goto fail;

    if (abs(nc) > csz) {
		PyErr_SetString(PpgTYPEErr,"contour: size of cont vec < than the "
						"req. contours number");
		goto fail;
    }
    if (trsz < 6) {
		PyErr_SetString(PpgTYPEErr,"contour: invalid transform. vector");
		goto fail;
    }

    switch (ft) {
    case FUN_PGCONB:
		cpgconb(a,cd,rd,c1+1,c2+1,r1+1,r2+1,c,nc,tr,blank);
		break;
    case FUN_PGCONS:
		cpgcons(a,cd,rd,c1+1,c2+1,r1+1,r2+1,c,nc,tr);
		break;
    case FUN_PGCONT:
		cpgcont(a,cd,rd,c1+1,c2+1,r1+1,r2+1,c,nc,tr);
		break;
    default:
		assert(0);
		break;
    }

    Py_DECREF(aa);
    Py_DECREF(ac);
    Py_DECREF(atr);
    PYRN;

fail:
    if (aa) { Py_DECREF(aa); }
    if (ac) { Py_DECREF(ac); }
    if (atr) { Py_DECREF(atr); }
    return(NULL);
}

PYF(pgconb)
{
    return(genContours(FUN_PGCONB,args));
}

PYF(pgcons)
{
    return(genContours(FUN_PGCONS,args));
}

PYF(pgcont)
{
    return(genContours(FUN_PGCONT,args));
}

/********************************************************************/

static PyObject *
genContours_s (enum pp_contour_funcs ft, PyObject *args)
{
    PyObject *oa=NULL, *oc=NULL;
    PyArrayObject *aa=NULL, *ac=NULL;
    float *a = NULL, *c = NULL, tr[6],
		x1=0.0,y1=0.0,x2=0.0,y2=0.0,blank=0.0,
		mn = 0.0, mx = 0.0;
    int rd=0, cd=0, csz=0, nc=0, ncont=0;

    if (!PyArg_ParseTuple(args,"Oi|Offfff:contour_s",
						  &oa,&nc,&oc,&x1,&y1,&x2,&y2,&blank))
		return(NULL);

    if (abs(nc)<1) {
		PyErr_SetString(PpgTYPEErr,"_ppgplot.error: Number of contours is 0");
		return(NULL);
    }
    if (!(aa = (PyArrayObject *)tofloatmat(oa, &a, &rd, &cd)))
		goto fail;
    if (oc) {
		if (!(ac = (PyArrayObject *)tofloatvector(oc, &c, &csz)))
			goto fail;
    } else {
		if (!(c = malloc(abs(nc)*sizeof(*c)))) {
			PyErr_SetString(PpgTYPEErr,"_ppgplot.error: Out of mem!");
			goto fail;
		}
		ncont = abs(nc);
    }

    /* Perform autocalibrations as nesecairy. */
    autocal2d(a, rd, cd, &mx, &mn, ncont, c, &x1, &x2, &y1, &y2, tr);

#ifdef DEBUG_CONT_S
    {
		int i;
		fprintf(stderr,"ncontours = %d = %d\n",nc,ncont);
		fprintf(stderr,"Contours:\n");
		for (i=0; i<abs(nc); i++)
			fprintf(stderr,"   cont[%d] = %f\n",i,c[i]);
		fprintf(stderr,"blank = %f\n",blank);
    }
#endif


    switch (ft) {
    case FUN_PGCONB:
		cpgconb(a,cd,rd,0+1,cd,0+1,rd,c,nc,tr,blank);
		break;
    case FUN_PGCONS:
		cpgcons(a,cd,rd,0+1,cd,0+1,rd,c,nc,tr);
		break;
    case FUN_PGCONT:
		cpgcont(a,cd,rd,0+1,cd,0+1,rd,c,nc,tr);
		break;
    default:
		assert(0);
		break;
    }

    Py_DECREF(aa);
    if (ac)
		Py_DECREF(ac);
    else if (c)
		free(c);
    PYRN;

fail:
    if (aa) { Py_DECREF(aa); }
    if (ac) { 
		Py_DECREF(ac);
    } else if (c) {
		free(c);
	}

    return(NULL);
}

PYF(pgconb_s)
{
    return(genContours_s(FUN_PGCONB,args));
}

PYF(pgcons_s)
{
    return(genContours_s(FUN_PGCONS,args));
}

PYF(pgcont_s)
{
    return(genContours_s(FUN_PGCONT,args));
}

PYF(pgconl)
{
    PyObject
		*oa=NULL, *otr=NULL;
    PyArrayObject
		*aa=NULL, *atr=NULL;
    float *a=NULL, *tr=NULL, c = 0.0;
    int cd=0 ,rd=0,
		c1=0,c2=0,r1=0,r2=0,
		trsz=0,
		intval = 20, minint = 10;
    char *label = NULL;

    if (!PyArg_ParseTuple(args,"OiiiiiifOs|ii:pgconl",
						  &oa, &cd, &rd, &c1, &c2, &r1, &r2,
						  &c, &otr, &label, &intval, &minint))
		return(NULL);

    if (!(aa = (PyArrayObject *)tofloatmat(oa, &a, &rd, &cd)))
		goto fail;
    if (!(atr = (PyArrayObject *)tofloatvector(otr, &tr, &trsz)))
		goto fail;

    if (trsz < 6) {
		PyErr_SetString(PpgTYPEErr,"contour: invalid transform. vector");
		goto fail;
    }

    cpgconl(a,cd,rd,c1+1,c2+1,r1+1,r2+1,c,tr,label,intval,minint);

    Py_DECREF(aa);
    Py_DECREF(atr);
    PYRN;

fail:
    if (aa) { Py_DECREF(aa); }
    if (atr) { Py_DECREF(atr); }
    return(NULL);
}


PYF(pgconl_s)
{
    PyObject
		*oa=NULL;
    PyArrayObject
		*aa=NULL;
    float *a=NULL, tr[6], c = 0.0,
		/* just to avoid eval. of min and max during autocal. */
		mn = 0.0, mx = 1.0,
		x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
    int cd=0 ,rd=0,
		intval = 20, minint = 10;
    char *label = NULL;

    if (!PyArg_ParseTuple(args,"Ofs|iiffff:pgconl",
						  &oa, &c, &label, &intval, &minint, &x1, &x2, 
						  &y1, &y2))
		return(NULL);
        
    if (!(aa = (PyArrayObject *)tofloatmat(oa, &a, &rd, &cd)))
		return(NULL);

    /* Perform autocalibrations as nesecairy. */
    autocal2d(a, rd, cd, &mx, &mn, 0, NULL, &x1, &x2, &y1, &y2, tr);    

    cpgconl(a,cd,rd,0+1,cd,0+1,rd,c,tr,label,intval,minint);

    Py_DECREF(aa);
    PYRN;
}

/***************************************************************************/

static PyMethodDef PpgMethods[] = {
    {"pgarro", pgarro, 1},
    {"pgask", pgask,1},
    {"pgband", pgband, 1},
    {"pgbbuf", pgbbuf, 1},
    {"pgbeg",  pgbeg, 1},
    {"pgbin", pgbin, 1},
    {"pgbin_s", pgbin_s, 1},
    {"pgbox", pgbox, 1},
    {"pgcirc", pgcirc, 1},
    {"pgclos", pgclos, 1},
    {"pgconb",pgconb,1},
    {"pgconb_s",pgconb_s,1},
    {"pgconl", pgconl, 1},
    {"pgconl_s", pgconl_s, 1},
    {"pgcons",pgcons,1},
    {"pgcons_s",pgcons_s,1},
    {"pgcont",pgcont,1},    
    {"pgcont_s",pgcont_s,1},
    {"pgctab",pgctab, 1},
    {"pgcurs", pgcurs, 1},
    {"pgdraw", pgdraw, 1},
    {"pgebuf", pgebuf, 1},
    {"pgend", pgend, 1},
    {"pgenv", pgenv, 1},
    {"pgeras", pgeras, 1},
    {"pgerrb", pgerrb, 1},
    {"pgerrx", pgerrx, 1},
    {"pgerry", pgerry, 1},
    {"pgetxt", pgetxt, 1},
    {"pggray", pggray, 1},
    {"pggray_s", pggray_s, 1},
    {"pghi2d", pghi2d, 1},
    {"pghi2d_s", pghi2d_s, 1},
    {"pghist", pghist, 1},
    {"pghist_s", pghist_s, 1},
    {"pgiden", pgiden, 1},
    {"pgimag", pgimag, 1},
    {"pgimag_s", pgimag_s, 1},
    {"pglab", pglab, 1},
    {"pgldev", pgldev, 1},
    {"pglen", pglen, 1},
    {"pgline", pgline, 1},
    {"pgmove", pgmove, 1},
    {"pgmtxt", pgmtxt, 1},
    {"pgopen",  pgopen, 1},
    {"pgpage", pgpage, 1},
    {"pgpanl", pgpanl, 1},
    {"pgpap", pgpap, 1},
    {"pgpt", pgpt, 1},
    {"pgptxt", pgptxt, 1},
    {"pgqah", pgqah, 1},
    {"pgqcf", pgqcf, 1},
    {"pgqch", pgqch, 1},
    {"pgqci", pgqci, 1},
    {"pgqcir", pgqcir, 1},
    {"pgqcol", pgqcol, 1},
    {"pgqcr", pgqcr, 1},
    {"pgqcs", pgqcs, 1},
    {"pgqfs", pgqfs, 1},
    {"pgqhs", pgqhs, 1},
    {"pgqid", pgqid, 1},
    {"pgqitf", pgqitf, 1},
    {"pgqls", pgqls, 1},
    {"pgqlw", pgqlw, 1},
    {"pgqpos", pgqpos, 1},
    {"pgrect", pgrect, 1},
    {"pgsah", pgsah, 1},
    {"pgsave", pgsave, 1},
    {"pgscf", pgscf, 1},
    {"pgsch", pgsch, 1},
    {"pgsci", pgsci, 1},
    {"pgscir", pgscir,1},
    {"pgscr", pgscr, 1},
    {"pgsfs", pgsfs, 1},
    {"pgshls", pgshls, 1},
    {"pgshs", pgshs, 1},
    {"pgsitf", pgsitf, 1},
    {"pgslct", pgslct, 1},
    {"pgsls", pgsls,1},
    {"pgslw", pgslw,1},
    {"pgstbg", pgstbg, 1},
    {"pgsubp", pgsubp, 1},
    {"pgsvp",pgsvp,1},
    {"pgswin",pgswin,1},
    {"pgtbox", pgtbox, 1},
    {"pgtext", pgtext, 1},
    {"pgunsa", pgunsa, 1},
    {"pgupdt", pgupdt, 1},
    {"pgvsiz", pgvsiz, 1},
    {"pgvstd", pgvstd, 1},
    {"pgwedg", pgwedg, 1},
    {"pgwedg_s", pgwedg_s, 1},
    {"pgwnad", pgwnad, 1},
	/* ADDED */
    {"pgqinf", pgqinf, 1},
    {"pgpoly", pgpoly, 1},
    {"pgqvp", pgqvp, 1},
    {"pgqvsz", pgqvsz, 1},
    {"pgsclp", pgsclp, 1},
    {"pgqclp", pgqclp, 1},
    {"pgconf", pgconf, 1},
    {"pgqtxt", pgqtxt, 1},

/* END ADDED */
#ifdef DEBUG
    {"tstmat", tstmat, 1},
#endif

    {NULL,      NULL}        /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef ppgplotdef = {
        PyModuleDef_HEAD_INIT,
        "_ppgplot",     /* m_name */
        "PPGPLOT Module",  /* m_doc */
        -1,                  /* m_size */
        PpgMethods,    /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
#endif

/************************************************************************/

static PyObject *
moduleinit(void)
{
    PyObject *m, *d;
#if PY_MAJOR_VERSION <= 2
    m = Py_InitModule("_ppgplot", PpgMethods);
#else
    m = PyModule_Create(&ppgplotdef);
#endif
    d = PyModule_GetDict(m);
    import_array();
#if PY_MAJOR_VERSION <= 2
    PpgIOErr = PyString_FromString("_ppgplot.ioerror");
    PpgTYPEErr = PyString_FromString("_ppgplot.typeerror");
    PpgMEMErr = PyString_FromString("_ppgplot.memerror");
#else
    PpgIOErr = PyBytes_FromString("_ppgplot.ioerror");
    PpgTYPEErr = PyBytes_FromString("_ppgplot.typeerror");
    PpgMEMErr = PyBytes_FromString("_ppgplot.memerror");
#endif
    PyDict_SetItemString(d, "ioerror", PpgIOErr);
    PyDict_SetItemString(d, "typeerror", PpgTYPEErr);
    PyDict_SetItemString(d, "memerror", PpgMEMErr);
    return m;
}

#if PY_MAJOR_VERSION < 3
    void
    init_ppgplot(void)
    {
        moduleinit();
    }
#else
    PyMODINIT_FUNC
    PyInit__ppgplot(void)
    {
        return moduleinit();
    }
#endif
/************************************************************************/
/* End of _ppgplot.c */
/************************************************************************/


