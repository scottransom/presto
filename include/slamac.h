#ifndef SLAMACHDEF
#define SLAMACHDEF

/*
   **  - - - - - - - - -
   **   s l a m a c . h
   **  - - - - - - - - -
   **
   **  Macros used by slalib library.
   **
   **  Last revision:   5 June 1996
   **
   **  Copyright P.T.Wallace.  All rights reserved.
 */

/* max(A,B) - larger (most +ve) of two numbers (generic) */
#define gmax(A,B) ((A)>(B)?(A):(B))

/* min(A,B) - smaller (least +ve) of two numbers (generic) */
#define gmin(A,B) ((A)<(B)?(A):(B))

/* dint(A) - truncate to nearest whole number towards zero (double) */
#define dint(A) ((A)<0.0?ceil(A):floor(A))

/* aint(A) - truncate to nearest whole number towards zero (float) */
#define aint(A) ((A)<0.0f?(float)ceil((double)(A)):(float)floor((double)(A)))

/* dnint(A) - round to nearest whole number (double) */
#define dnint(A) ((A)<0.0?ceil((A)-0.5):floor((A)+0.5))

/* anint(A) - round to nearest whole number (float) */
#define anint(A) ((float)dnint((double)(A)))

/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-(A):(A))

/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))

/* logicals */
#if !defined(FALSE) || ((FALSE)!=0)
#define FALSE 0
#endif
#if !defined(TRUE) || ((TRUE)!=1)
#define TRUE 1
#endif

/* pi */
#define DPI 3.1415926535897932384626433832795028841971693993751

/* 2pi */
#define D2PI 6.2831853071795864769252867665590057683943387987502

/* 1/(2pi) */
#define D1B2PI 0.15915494309189533576888376337251436203445964574046

/* 4pi */
#define D4PI 12.566370614359172953850573533118011536788677597500

/* 1/(4pi) */
#define D1B4PI 0.079577471545947667884441881686257181017229822870228

/* pi^2 */
#define DPISQ 9.8696044010893586188344909998761511353136994072408

/* sqrt(pi) */
#define DSQRPI 1.7724538509055160272981674833411451827975494561224

/* pi/2:  90 degrees in radians */
#define DPIBY2 1.5707963267948966192313216916397514420985846996876

/* pi/180:  degrees to radians */
#define DD2R 0.017453292519943295769236907684886127134428718885417

/* 180/pi:  radians to degrees */
#define DR2D 57.295779513082320876798154814105170332405472466564

/* pi/(180*3600):  arcseconds to radians */
#define DAS2R 4.8481368110953599358991410235794797595635330237270e-6

/* 180*3600/pi :  radians to arcseconds */
#define DR2AS 2.0626480624709635515647335733077861319665970087963e5

/* pi/12:  hours to radians */
#define DH2R 0.26179938779914943653855361527329190701643078328126

/* 12/pi:  radians to hours */
#define DR2H 3.8197186342054880584532103209403446888270314977709

/* pi/(12*3600):  seconds of time to radians */
#define DS2R 7.2722052166430399038487115353692196393452995355905e-5

/* 12*3600/pi:  radians to seconds of time */
#define DR2S 1.3750987083139757010431557155385240879777313391975e4

/* 15/(2pi):  hours to degrees x radians to turns */
#define D15B2P 2.3873241463784300365332564505877154305168946861068

#endif
