/*********************************************************************
GH/02.10.95

include file for type definition "real"

*********************************************************************/
/*
  Definitions in connection with the type of "real"
*/

/* define DOUBLE if neither FLOAT nor DOUBLE have been defined before. */
#ifndef REAL_IS_DOUBLE
#ifndef REAL_IS_FLOAT
#define REAL_IS_DOUBLE
#endif
#endif

/* alternatives:
#define REAL_IS_DOUBLE
#define REAL_IS_FLOAT
*/

/************************************************
*************************************************
    REAL IS DOUBLE
*************************************************
************************************************/
#ifdef REAL_IS_DOUBLE

#ifndef TYPE_REAL
#define TYPE_REAL
 typedef double  real;
#endif

/*
  Define macros for math operations of type double
  (standard math library)
*/

#define R_cos cos
#define R_cosd cosd
#define R_sin sin
#define R_sind sind
#define R_asin asin
#define R_acos acos
#define R_atan2 atan2

#define R_cosh cosh
#define R_sinh sinh

#define R_exp exp
#define R_log log

#define R_fabs fabs
#define R_sqrt sqrt

#define R_cabs hypot
#define R_hypot hypot

#define R_nint rint
#define R_rint rint

#endif  /* REAL_IS_DOUBLE */

/************************************************
*************************************************
    REAL IS FLOAT
*************************************************
************************************************/

#ifdef REAL_IS_FLOAT

#ifndef TYPE_REAL
#define TYPE_REAL
 typedef float  real;
#endif

/*
  Define macros for math operations of type float
*/

#ifdef __osf__                                     /* DEC */

#define R_cos cosf
#define R_cosd cosdf
#define R_sin sinf
#define R_sind sindf
#define R_asin asinf
#define R_acos acosf
#define R_atan2 atan2f

#define R_cosh coshf
#define R_sinh sinhf

#define R_exp expf
#define R_log logf

#define R_fabs fabsf
#define R_sqrt sqrtf

#define R_cabs hypotf
#define R_hypot hypotf

#define R_nint rintf

#else                              /* IBM and others (standard math library) */

#define R_cos(x)   (float) cos  ((double)(x))
#define R_cosd(x)  (float) cosd ((double)(x))
#define R_sin(x)   (float) sin  ((double)(x))
#define R_sind(x)  (float) sind ((double)(x))
#define R_atan2(x,y)  (float)atan2((double)(x), (double)(y))

#define R_cosh(x)  (float) cosh ((double)(x))
#define R_sinh(x)  (float) sinh ((double)(x))

#define R_exp(x)   (float) exp  ((double)(x))
#define R_log(x)   (float) log  ((double)(x))

#define R_fabs(x)  (float) fabs ((double)(x))
#define R_sqrt(x)  (float) sqrt ((double)(x))

#define R_cabs(x,y)   (float) hypot ((double)(x), (double)(y))
#define R_hypot(x,y)  (float) hypot ((double)(x), (double)(y))

#define R_nint(x)  (float) rint ((double)(x))

#endif

#endif  /* REAL_IS_FLOAT */
