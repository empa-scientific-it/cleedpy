/*********************************************************************
GH/10.07.10

constants and macros generally used
GH/10.07.10: include <stdlib.h> and <string.h> to avoid compiler warnings

*********************************************************************/

#ifndef GH_STD_DEF_H
#define GH_STD_DEF_H

/*********************************************************************
 define machine
*********************************************************************/

/* Alternatives:
#define DEC
#define IBM
*/

#define IBM
/*
#define free(x) free((void *)(x))
*/

/*********************************************************************
*********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*********************************************************************
 output channels
*********************************************************************/


#define STDOUT stdout

#define STDERR stderr
#define STDWAR stderr
#define STDCTR stdout
#define STDCPU stdout

/*********************************************************************
 general mathematical definitions / constants
*********************************************************************/

#ifdef _XOPEN_SOURCE
#define PI  M_PI
#else
#define PI  3.1415926535897932385
#endif
/*
*/

#define DEG_TO_RAD 0.017453293         /* conversion degree to radian */
#define RAD_TO_DEG 57.29578            /* conversion radian to degree */

#define M2_H      0.2631894506957162   /* 2*m/h       [eV^-1   A^-2] */
#define SQRT_M2_H 0.5130199320647456   /* sqrt(2*m/h) [eV^-0.5 A^-1] */

/*********************************************************************
 general other definitions / constants
*********************************************************************/

#define KBYTE 1024
#define MBYTE 1048576

/*********************************************************************
 special definitions
*********************************************************************/

#define STRSZ 256                /* maximum length of strings */

#define I_END_OF_LIST   -9999    /* list terminator (integer)*/
#define F_END_OF_LIST   -9999.   /* list terminator (float)  */

#define IEND_OF_LIST   I_END_OF_LIST    /* alias for list terminator (integer)*/
#define FEND_OF_LIST   F_END_OF_LIST    /* alias for list terminator (float)  */

/*********************************************************************
 macros:
*********************************************************************/

#define MAX(x,y)  ((x)>(y))?(x):(y)
#define MIN(x,y)  ((x)<(y))?(x):(y)
#define SQUARE(x) (x)*(x)

#define ODD(n)    ((n)%2)
#define M1P(n)    (((n)%2)?(-1.):(1.))

#endif
