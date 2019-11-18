/*** File wcslib/imio.c
 *** May 20, 2003
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 1996-2003
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

 * Module:      imio.c (image pixel manipulation)
 * Purpose:     Read and write pixels from arbitrary data type 2D arrays
 * Subroutine:	getpix (image, bitpix, w, h, bz, bs, x, y)
 *		Read pixel from 2D image of any numeric type (0,0 lower left)
 * Subroutine:	getpix1 (image, bitpix, w, h, bz, bs, x, y)
 *		Read pixel from 2D image of any numeric type (1,1 lower left)
 * Subroutine:	putpix (image, bitpix, w, h, bz, bs, x, y, dpix)
 *		Write pixel into 2D image of any numeric type (0,0 lower left)
 * Subroutine:	putpix1 (image, bitpix, w, h, bz, bs, x, y, dpix)
 *		Write pixel into 2D image of any numeric type (1,1 lower left)
 * Subroutine:	addpix (image, bitpix, w, h, bz, bs, x, y, dpix)
 *		Copy pixel into 2D image of any numeric type (0,0 lower left)
 * Subroutine:	addpix1 (image, bitpix, w, h, bz, bs, x, y, dpix)
 *		Add pixel into 2D image of any numeric type (1,1 lower left)
 * Subroutine:	getvec (image, bitpix, bz, bs, pix1, npix, dpix)
 *		Get vector from 2D image of any numeric type
 * Subroutine:	putvec (image, bitpix, bz, bs, pix1, npix, dpix)
 *		Copy pixel vector into 2D image of any numeric type
 * Subroutine:	movepix (image1, bitpix, w1, x1, y1, image2, w2, x2, y2)
 *		Copy pixel from one image location to another
 * Subroutine:	imswap (bitpix,string,nbytes)
 *		Swap bytes in string in place, with FITS bits/pixel code
 * Subroutine:	imswap2 (string,nbytes)
 *		Swap bytes in string in place
 * Subroutine	imswap4 (string,nbytes)
 *		Reverse bytes of Integer*4 or Real*4 vector in place
 * Subroutine	imswap8 (string,nbytes)
 *		Reverse bytes of Real*8 vector in place
 * Subroutine	imswapped ()
 *		Return 1 if PC/DEC byte order, else 0
 */

#include <stdio.h>
#include <stdlib.h>
#include "imio.h"

static int scale = 1;           /* If 0, skip scaling step */
void setscale(scale0)
int scale0;
{
    scale = scale0;
    return;
}

/* GETPIX1 -- Get pixel from 2D FITS image of any numeric type */

double getpix1(image, bitpix, w, h, bzero, bscale, x, y)

char *image;                    /* Image array as 1-D vector */
int bitpix;                     /* FITS bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w;                          /* Image width in pixels */
int h;                          /* Image height in pixels */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int x;                          /* One-based horizontal pixel number */
int y;                          /* One-based vertical pixel number */

{
    return (getpix(image, bitpix, w, h, bzero, bscale, x - 1, y - 1));
}


/* GETPIX -- Get pixel from 2D image of any numeric type */

double getpix(image, bitpix, w, h, bzero, bscale, x, y)

char *image;                    /* Image array as 1-D vector */
int bitpix;                     /* FITS bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w;                          /* Image width in pixels */
int h;                          /* Image height in pixels */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int x;                          /* Zero-based horizontal pixel number */
int y;                          /* Zero-based vertical pixel number */

{
    short *im2;
    int *im4;
    unsigned short *imu;
    float *imr;
    double *imd;
    double dpix;

/* Return 0 if coordinates are not inside image */
    if (x < 0 || x >= w)
        return (0.0);
    if (y < 0 || y >= h)
        return (0.0);

/* Extract pixel from appropriate type of array */
    switch (bitpix) {

    case 8:
        dpix = (double) image[(y * w) + x];
        break;

    case 16:
        im2 = (short *) image;
        dpix = (double) im2[(y * w) + x];
        break;

    case 32:
        im4 = (int *) image;
        dpix = (double) im4[(y * w) + x];
        break;

    case -16:
        imu = (unsigned short *) image;
        dpix = (double) imu[(y * w) + x];
        break;

    case -32:
        imr = (float *) image;
        dpix = (double) imr[(y * w) + x];
        break;

    case -64:
        imd = (double *) image;
        dpix = imd[(y * w) + x];
        break;

    default:
        dpix = 0.0;
    }
    if (scale)
        return (bzero + (bscale * dpix));
    else
        return (dpix);
}


/* PUTPIX1 -- Copy pixel into 2D FITS image of any numeric type */

void putpix1(image, bitpix, w, h, bzero, bscale, x, y, dpix)

char *image;
int bitpix;                     /* Number of bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w;                          /* Image width in pixels */
int h;                          /* Image height in pixels */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int x;                          /* One-based horizontal pixel number */
int y;                          /* One-based vertical pixel number */
double dpix;

{
    putpix(image, bitpix, w, h, bzero, bscale, x - 1, y - 1, dpix);
    return;
}


/* PUTPIX -- Copy pixel into 2D image of any numeric type */

void putpix(image, bitpix, w, h, bzero, bscale, x, y, dpix)

char *image;
int bitpix;                     /* Number of bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w;                          /* Image width in pixels */
int h;                          /* Image height in pixels */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int x;
int y;
double dpix;

{
    short *im2;
    int *im4;
    unsigned short *imu;
    float *imr;
    double *imd;

/* Return if coordinates are not inside image */
    if (x < 0 || x >= w)
        return;
    if (y < 0 || y >= h)
        return;

    if (scale)
        dpix = (dpix - bzero) / bscale;

    switch (bitpix) {

    case 8:
        if (dpix < 0)
            image[(y * w) + x] = (char) (dpix - 0.5);
        else
            image[(y * w) + x] = (char) (dpix + 0.5);
        break;

    case 16:
        im2 = (short *) image;
        if (dpix < 0)
            im2[(y * w) + x] = (short) (dpix - 0.5);
        else
            im2[(y * w) + x] = (short) (dpix + 0.5);
        break;

    case 32:
        im4 = (int *) image;
        if (dpix < 0)
            im4[(y * w) + x] = (int) (dpix - 0.5);
        else
            im4[(y * w) + x] = (int) (dpix + 0.5);
        break;

    case -16:
        imu = (unsigned short *) image;
        if (dpix < 0)
            imu[(y * w) + x] = (unsigned short) 0;
        else
            imu[(y * w) + x] = (unsigned short) (dpix + 0.5);
        break;

    case -32:
        imr = (float *) image;
        imr[(y * w) + x] = (float) dpix;
        break;

    case -64:
        imd = (double *) image;
        imd[(y * w) + x] = dpix;
        break;

    }
    return;
}


/* ADDPIX1 -- Add pixel value into 2D FITS image of any numeric type */

void addpix1(image, bitpix, w, h, bzero, bscale, x, y, dpix)

char *image;
int bitpix;                     /* Number of bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w;                          /* Image width in pixels */
int h;                          /* Image height in pixels */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int x;                          /* One-based horizontal pixel number */
int y;                          /* One-based vertical pixel number */
double dpix;                    /* Value to add to pixel */

{
    addpix(image, bitpix, w, h, bzero, bscale, x - 1, y - 1, dpix);
    return;
}


/* ADDPIX -- Add pixel value into 2D image of any numeric type */

void addpix(image, bitpix, w, h, bzero, bscale, x, y, dpix)

char *image;
int bitpix;                     /* Number of bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w;                          /* Image width in pixels */
int h;                          /* Image height in pixels */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int x;                          /* Zero-based horizontal pixel number */
int y;                          /* Zero-based vertical pixel number */
double dpix;                    /* Value to add to pixel */

{
    short *im2;
    int *im4;
    unsigned short *imu;
    float *imr;
    double *imd;
    int ipix;

/* Return if coordinates are not inside image */
    if (x < 0 || x >= w)
        return;
    if (y < 0 || y >= h)
        return;

    if (scale)
        dpix = (dpix - bzero) / bscale;
    ipix = (y * w) + x;

    switch (bitpix) {

    case 8:
        if (dpix < 0)
            image[ipix] = image[ipix] + (char) (dpix - 0.5);
        else
            image[ipix] = image[ipix] + (char) (dpix + 0.5);
        break;

    case 16:
        im2 = (short *) image;
        if (dpix < 0)
            im2[ipix] = im2[ipix] + (short) (dpix - 0.5);
        else
            im2[ipix] = im2[ipix] + (short) (dpix + 0.5);
        break;

    case 32:
        im4 = (int *) image;
        if (dpix < 0)
            im4[ipix] = im4[ipix] + (int) (dpix - 0.5);
        else
            im4[ipix] = im4[ipix] + (int) (dpix + 0.5);
        break;

    case -16:
        imu = (unsigned short *) image;
        if (dpix > 0)
            imu[ipix] = imu[ipix] + (unsigned short) (dpix + 0.5);
        break;

    case -32:
        imr = (float *) image;
        imr[ipix] = imr[ipix] + (float) dpix;
        break;

    case -64:
        imd = (double *) image;
        imd[ipix] = imd[ipix] + dpix;
        break;

    }
    return;
}


/* MOVEPIX -- Copy pixel between images */

void movepix(image1, bitpix1, w1, x1, y1, image2, bitpix2, w2, x2, y2)

char *image1;                   /* Pointer to first pixel in input image */
int bitpix1;                    /* Bits per input pixel (FITS codes) */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w1;                         /* Number of horizontal pixels in input image */
int x1, y1;                     /* Row and column for input pixel */

char *image2;                   /* Pointer to first pixel in output image */
int bitpix2;                    /* Bits per output pixel (FITS codes) */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
int w2;                         /* Number of horizontal pixels in output image */
int x2, y2;                     /* Row and column for output pixel */

{
    short *ims1, *ims2;
    int *imi1, *imi2;
    unsigned short *imu1, *imu2;
    float rpix, *imr1, *imr2;
    double dpix, *imd1, *imd2;

    switch (bitpix1) {

    case 8:
        switch (bitpix2) {
        case 8:
            image2[(y2 * w2) + x2] = image1[(y1 * w1) + x1];
            break;
        case 16:
            ims2 = (short *) image2;
            ims2[(y2 * w2) + x2] = image1[(y1 * w1) + x1];
            break;
        case 32:
            imi2 = (int *) image2;
            imi2[(y2 * w2) + x2] = (int) image1[(y1 * w1) + x1];
            break;
        case -16:
            imu2 = (unsigned short *) image2;
            imu2[(y2 * w2) + x2] = (unsigned short) image1[(y1 * w1) + x1];
            break;
        case -32:
            imr2 = (float *) image2;
            imr2[(y2 * w2) + x2] = (float) image1[(y1 * w1) + x1];
            break;
        case -64:
            imd2 = (double *) image2;
            imd2[(y2 * w2) + x2] = (double) image1[(y1 * w1) + x1];
            break;
        }
        break;

    case 16:
        switch (bitpix2) {
        case 8:
            ims1 = (short *) image1;
            image2[(y2 * w2) + x2] = (char) ims1[(y1 * w1) + x1];
            break;
        case 16:
            ims1 = (short *) image1;
            ims2 = (short *) image2;
            ims2[(y2 * w2) + x2] = ims1[(y1 * w1) + x1];
            break;
        case 32:
            ims1 = (short *) image1;
            imi2 = (int *) image2;
            imi2[(y2 * w2) + x2] = (int) ims1[(y1 * w1) + x1];
            break;
        case -16:
            ims1 = (short *) image1;
            imu2 = (unsigned short *) image2;
            imu2[(y2 * w2) + x2] = (unsigned short) ims1[(y1 * w1) + x1];
            break;
        case -32:
            ims1 = (short *) image1;
            imr2 = (float *) image2;
            imr2[(y2 * w2) + x2] = (float) ims1[(y1 * w1) + x1];
            break;
        case -64:
            ims1 = (short *) image1;
            imd2 = (double *) image2;
            imd2[(y2 * w2) + x2] = (double) ims1[(y1 * w1) + x1];
            break;
        }
        break;

    case 32:
        switch (bitpix2) {
        case 8:
            imi1 = (int *) image1;
            image2[(y2 * w2) + x2] = (char) imi1[(y1 * w1) + x1];
            break;
        case 16:
            imi1 = (int *) image1;
            ims2 = (short *) image2;
            ims2[(y2 * w2) + x2] = (short) imi1[(y1 * w1) + x1];
            break;
        case 32:
            imi1 = (int *) image1;
            imi2 = (int *) image2;
            imi2[(y2 * w2) + x2] = imi1[(y1 * w1) + x1];
            break;
        case -16:
            imi1 = (int *) image1;
            imu2 = (unsigned short *) image2;
            imu2[(y2 * w2) + x2] = (unsigned short) imi1[(y1 * w1) + x1];
            break;
        case -32:
            imi1 = (int *) image1;
            imr2 = (float *) image2;
            imr2[(y2 * w2) + x2] = (float) imi1[(y1 * w1) + x1];
            break;
        case -64:
            imi1 = (int *) image1;
            imd2 = (double *) image2;
            imd2[(y2 * w2) + x2] = (double) imi1[(y1 * w1) + x1];
            break;
        }
        break;

    case -16:
        switch (bitpix2) {
        case 8:
            imu1 = (unsigned short *) image1;
            image2[(y2 * w2) + x2] = (char) imu1[(y1 * w1) + x1];
            break;
        case 16:
            imu1 = (unsigned short *) image1;
            ims2 = (short *) image2;
            ims2[(y2 * w2) + x2] = (short) imu1[(y1 * w1) + x1];
            break;
        case 32:
            imu1 = (unsigned short *) image1;
            imi2 = (int *) image2;
            imi2[(y2 * w2) + x2] = (int) imu1[(y1 * w1) + x1];
            break;
        case -16:
            imu1 = (unsigned short *) image1;
            imu2 = (unsigned short *) image2;
            imu2[(y2 * w2) + x2] = imu1[(y1 * w1) + x1];
            break;
        case -32:
            imu1 = (unsigned short *) image1;
            imr2 = (float *) image2;
            imr2[(y2 * w2) + x2] = (float) imu1[(y1 * w1) + x1];
            break;
        case -64:
            imu1 = (unsigned short *) image1;
            imd2 = (double *) image2;
            imd2[(y2 * w2) + x2] = (double) imu1[(y1 * w1) + x1];
            break;
        }
        break;

    case -32:
        imr1 = (float *) image1;
        rpix = imr1[(y1 * w1) + x1];
        switch (bitpix2) {
        case 8:
            if (rpix < 0.0)
                image2[(y2 * w2) + x2] = (char) (rpix - 0.5);
            else
                image2[(y2 * w2) + x2] = (char) (rpix + 0.5);
            break;
        case 16:
            ims2 = (short *) image2;
            if (rpix < 0.0)
                ims2[(y2 * w2) + x2] = (short) (rpix - 0.5);
            else
                ims2[(y2 * w2) + x2] = (short) (rpix + 0.5);
            break;
        case 32:
            imi2 = (int *) image2;
            if (rpix < 0.0)
                imi2[(y2 * w2) + x2] = (int) (rpix - 0.5);
            else
                imi2[(y2 * w2) + x2] = (int) (rpix + 0.5);
            break;
        case -16:
            imu2 = (unsigned short *) image2;
            if (rpix < 0.0)
                imu2[(y2 * w2) + x2] = (unsigned short) 0;
            else
                imu2[(y2 * w2) + x2] = (unsigned short) (rpix + 0.5);
            break;
        case -32:
            imr2 = (float *) image2;
            imr2[(y2 * w2) + x2] = rpix;
            break;
        case -64:
            imd2 = (double *) image2;
            imd2[(y2 * w2) + x2] = (double) rpix;
            break;
        }
        break;

    case -64:
        imd1 = (double *) image1;
        dpix = imd1[(y1 * w1) + x1];
        switch (bitpix2) {
        case 8:
            imd1 = (double *) image1;
            if (dpix < 0.0)
                image2[(y2 * w2) + x2] = (char) (dpix - 0.5);
            else
                image2[(y2 * w2) + x2] = (char) (dpix + 0.5);
            break;
        case 16:
            ims2 = (short *) image2;
            if (dpix < 0.0)
                ims2[(y2 * w2) + x2] = (short) (dpix - 0.5);
            else
                ims2[(y2 * w2) + x2] = (short) (dpix + 0.5);
            break;
        case 32:
            imi2 = (int *) image2;
            if (dpix < 0.0)
                imi2[(y2 * w2) + x2] = (int) (dpix - 0.5);
            else
                imi2[(y2 * w2) + x2] = (int) (dpix + 0.5);
            break;
        case -16:
            imu2 = (unsigned short *) image2;
            if (dpix < 0.0)
                imu2[(y2 * w2) + x2] = (unsigned short) 0;
            else
                imu2[(y2 * w2) + x2] = (unsigned short) (dpix + 0.5);
            break;
        case -32:
            imr2 = (float *) image2;
            imr2[(y2 * w2) + x2] = (float) dpix;
            break;
        case -64:
            imd2 = (double *) image2;
            imd2[(y2 * w2) + x2] = dpix;
            break;
        }
        break;
    }
    return;
}


/* GETVEC -- Get vector from 2D image of any numeric type */

void getvec(image, bitpix, bzero, bscale, pix1, npix, dvec0)

char *image;                    /* Image array from which to extract vector */
int bitpix;                     /* Number of bits per pixel in image */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int pix1;                       /* Offset of first pixel to extract */
int npix;                       /* Number of pixels to extract */
double *dvec0;                  /* Vector of pixels (returned) */

{
    short *im2;
    int *im4;
    unsigned short *imu;
    float *imr;
    double *imd;
    double *dvec;
    int ipix, pix2;

    pix2 = pix1 + npix;
    dvec = dvec0;

    switch (bitpix) {

    case 8:
        for (ipix = pix1; ipix < pix2; ipix++)
            *dvec++ = (double) *(image + ipix);
        break;

    case 16:
        im2 = (short *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *dvec++ = (double) *(im2 + ipix);
        break;

    case 32:
        im4 = (int *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *dvec++ = bscale * (double) *(im4 + ipix);
        break;

    case -16:
        imu = (unsigned short *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *dvec++ = (double) *(imu + ipix);
        break;

    case -32:
        imr = (float *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *dvec++ = (double) *(imr + ipix);
        break;

    case -64:
        imd = (double *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *dvec++ = (double) *(imd + ipix);
        break;

    }

    /* Scale data if either BZERO or BSCALE keyword has been set */
    if (scale && (bzero != 0.0 || bscale != 1.0)) {
        dvec = dvec0;
        for (ipix = pix1; ipix < pix2; ipix++) {
            *dvec = (*dvec * bscale) + bzero;
            dvec++;
        }
    }

    return;
}


/* PUTVEC -- Copy pixel vector into 2D image of any numeric type */

void putvec(image, bitpix, bzero, bscale, pix1, npix, dvec)

char *image;                    /* Image into which to copy vector */
int bitpix;                     /* Number of bits per pixel im image */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
double bzero;                   /* Zero point for pixel scaling */
double bscale;                  /* Scale factor for pixel scaling */
int pix1;                       /* Offset of first pixel of vector in image */
int npix;                       /* Number of pixels to copy */
double *dvec;                   /* Vector of pixels to copy */

{
    short *im2;
    int *im4;
    unsigned short *imu;
    float *imr;
    double *imd;
    int ipix, pix2;
    double *dp = dvec;

    pix2 = pix1 + npix;

    /* Scale data if either BZERO or BSCALE keyword has been set */
    if (scale && (bzero != 0.0 || bscale != 1.0)) {
        for (ipix = pix1; ipix < pix2; ipix++) {
            *dp = (*dp - bzero) / bscale;
            dp++;
        }
        dp = dvec;
    }

    switch (bitpix) {

    case 8:
        for (ipix = pix1; ipix < pix2; ipix++)
            *(image + ipix) = (char) *dp++;
        break;

    case 16:
        im2 = (short *) image;
        for (ipix = pix1; ipix < pix2; ipix++) {
            if (*dp < 0.0)
                *(im2 + ipix) = (short) (*dp++ - 0.5);
            else
                *(im2 + ipix) = (short) (*dp++ + 0.5);
        }
        break;

    case 32:
        im4 = (int *) image;
        for (ipix = pix1; ipix < pix2; ipix++) {
            if (*dp < 0.0)
                *(im4 + ipix) = (int) (*dp++ - 0.5);
            else
                *(im4 + ipix) = (int) (*dp++ + 0.5);
        }
        break;

    case -16:
        imu = (unsigned short *) image;
        for (ipix = pix1; ipix < pix2; ipix++) {
            if (*dp < 0.0)
                *(imu + ipix) = (unsigned short) 0;
            else
                *(imu + ipix) = (unsigned short) (*dp++ + 0.5);
        }
        break;

    case -32:
        imr = (float *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *(imr + ipix) = (float) *dp++;
        break;

    case -64:
        imd = (double *) image;
        for (ipix = pix1; ipix < pix2; ipix++)
            *(imd + ipix) = (double) *dp++;
        break;
    }
    return;
}


/* IMSWAP -- Reverse bytes of any type of vector in place */

void imswap(bitpix, string, nbytes)

int bitpix;                     /* Number of bits per pixel */
                        /*  16 = short, -16 = unsigned short, 32 = int */
                        /* -32 = float, -64 = double */
char *string;                   /* Address of starting point of bytes to swap */
int nbytes;                     /* Number of bytes to swap */

{
    switch (bitpix) {

    case 8:
        break;

    case 16:
        if (nbytes < 2)
            return;
        imswap2(string, nbytes);
        break;

    case 32:
        if (nbytes < 4)
            return;
        imswap4(string, nbytes);
        break;

    case -16:
        if (nbytes < 2)
            return;
        imswap2(string, nbytes);
        break;

    case -32:
        if (nbytes < 4)
            return;
        imswap4(string, nbytes);
        break;

    case -64:
        if (nbytes < 8)
            return;
        imswap8(string, nbytes);
        break;

    }
    return;
}


/* IMSWAP2 -- Swap bytes in string in place */

void imswap2(string, nbytes)


char *string;                   /* Address of starting point of bytes to swap */
int nbytes;                     /* Number of bytes to swap */

{
    char *sbyte, temp, *slast;

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
        temp = sbyte[0];
        sbyte[0] = sbyte[1];
        sbyte[1] = temp;
        sbyte = sbyte + 2;
    }
    return;
}


/* IMSWAP4 -- Reverse bytes of Integer*4 or Real*4 vector in place */

void imswap4(string, nbytes)

char *string;                   /* Address of Integer*4 or Real*4 vector */
int nbytes;                     /* Number of bytes to reverse */

{
    char *sbyte, *slast;
    char temp0, temp1, temp2, temp3;

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
        temp3 = sbyte[0];
        temp2 = sbyte[1];
        temp1 = sbyte[2];
        temp0 = sbyte[3];
        sbyte[0] = temp0;
        sbyte[1] = temp1;
        sbyte[2] = temp2;
        sbyte[3] = temp3;
        sbyte = sbyte + 4;
    }

    return;
}


/* IMSWAP8 -- Reverse bytes of Real*8 vector in place */

void imswap8(string, nbytes)

char *string;                   /* Address of Real*8 vector */
int nbytes;                     /* Number of bytes to reverse */

{
    char *sbyte, *slast;
    char temp[8];

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
        temp[7] = sbyte[0];
        temp[6] = sbyte[1];
        temp[5] = sbyte[2];
        temp[4] = sbyte[3];
        temp[3] = sbyte[4];
        temp[2] = sbyte[5];
        temp[1] = sbyte[6];
        temp[0] = sbyte[7];
        sbyte[0] = temp[0];
        sbyte[1] = temp[1];
        sbyte[2] = temp[2];
        sbyte[3] = temp[3];
        sbyte[4] = temp[4];
        sbyte[5] = temp[5];
        sbyte[6] = temp[6];
        sbyte[7] = temp[7];
        sbyte = sbyte + 8;
    }
    return;
}

/* IMSWAPPED -- Returns 0 if big-endian (Sun,Mac),
		1 if little-endian(PC,Alpha) */

int imswapped()
{
    char *ctest;
    int itest;

    itest = 1;
    ctest = (char *) &itest;
    if (*ctest)
        return (1);
    else
        return (0);
}

/* Apr 17 1996	New file
 * May 22 1996	Add H so that PUTPIX and GETPIX can check coordinates
 * Jun 11 1996	Simplify NEWIMAGE subroutine
 * Jun 12 1996	Add byte-swapping subroutines
 *
 * Jul 24 1997	Add 8-bit option to subroutines
 *
 * May 27 1998	Include imio.h instead of fitshead.h
 * Jun 17 1998	Fix bug, changing all unsigned int's to unsigned short's
 *
 * Apr 29 1999	Add scaling to getpix, putpix, getvec, and putvec
 * Apr 29 1999	Fix bug in getvec in dealing with 1-byte data
 * Sep 14 1999	Change dp incrementing so it works on Alpha compiler
 * Sep 27 1999	Add interface for 1-based (FITS) image access
 * Sep 27 1999	Add addpix() and addpix1()
 * Dec 14 1999	In putpix(), addpix(), putvec(), round when output is integer
 *
 * Sep 20 2000	In getvec(), scale only if necessary
 *
 * Nov 27 2001	In movepix(), add char to char move
 *
 * Jan 23 2002	Add global scale switch to turn off scaling
 * Jun  4 2002	In getvec() and putvec(), change dpix to dvec
 * Jun  4 2002	Add addvec() to add to a vector
 * Jul 19 2002	Fix getvec() bug rescaling scaled numbers
 *
 * May 20 2003	Declare scale0 in setscale()
 */
