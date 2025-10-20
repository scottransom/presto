/*** imio.h  memory access subroutines
 *** September 27, 1999
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

#ifndef imio_h_
#define imio_h_

/* Image pixel access subroutines in imio.c */
extern double getpix(char *image, int bitpix, int w, int h, double bzero, double bscale, int x, int y); /* Read one pixel from any data type 2-D array (0,0)*/
extern double getpix1(char *image, int bitpix, int w, int h, double bzero, double bscale, int x, int y); /* Read one pixel from any data type 2-D array (1,1)*/
extern void putpix(char *image, int bitpix, int w, int h, double bzero, double bscale, int x, int y, double dpix);   /* Write one pixel to any data type 2-D array (0,0)*/
extern void putpix1(char *image, int bitpix, int w, int h, double bzero, double bscale, int x, int y, double dpix);  /* Write one pixel to any data type 2-D array (1,1) */
extern void addpix(char *image, int bitpix, int w, int h, double bzero, double bscale, int x, int y, double dpix);   /* Add to one pixel in any data type 2-D array (0,0)*/
extern void addpix1(char *image, int bitpix, int w, int h, double bzero, double bscale, int x, int y, double dpix);  /* Add to one pixel in any data type 2-D array (1,1)*/
extern void movepix(char *image1, int bitpix1, int w1, int x1, int y1, char *image2, int bitpix2, int w2, int x2, int y2);  /* Move one pixel value between two 2-D arrays (0,0) */
extern void getvec(char *image, int bitpix, double bzero, double bscale, int pix1, int npix, double *dvec0);   /* Read vector from 2-D array */
extern void putvec(char *image, int bitpix, double bzero, double bscale, int pix1, int npix, double *dvec);   /* Write vector into 2-D array */
extern void imswap(int bitpix, char *string, int nbytes);   /* Swap alternating bytes in a vector */
extern void imswap2(char *string, int nbytes);  /* Swap bytes in a vector of 2-byte (short) integers */
extern void imswap4(char *string, int nbytes);  /* Reverse bytes in a vector of 4-byte numbers */
extern void imswap8(char *string, int nbytes);  /* Reverse bytes in a vector of 8-byte numbers */
extern int imswapped(void); /* Return 1 if machine byte order is not FITS order */

#endif	/* imio_h_ */

/* May 31 1996	Use stream I/O for reading as well as writing
 * Jun 12 1996	Add byte-swapping subroutines
 * Aug  6 1996	Add MOVEPIX, HDEL and HCHANGE declarations
 *
 * May 27 1998	Split off imio subroutines to imio.h

 * Sep 27 1999	Add Fortran-indexed (1,1), not (0,0) image access *1()
 * Sep 28 1999	Add addpix()
 */
