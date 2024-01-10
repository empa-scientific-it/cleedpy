/*********************************************************************
GH/28.09.00
  file contains function:

  inp_debtemp

!!!
!!! This function is used in CLEED and in CSEARCH
!!!

Changes:

GH/27.09.00 - bug fix in Debye Waller-Temp input
*********************************************************************/

#include <math.h>
#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif
#include <stdio.h>
#include <strings.h>

#include "leed.h"

/*
#define CONTROL_X
#define CONTROL
*/
#define CONTROL
#define WARNING
#define ERROR

#define EXIT_ON_ERROR

#define PREF_DEBWAL 1559.04170632481439  /* prefactor for the evaluation
                                            of <dr^2> from temperature
                                            and Debye temperature */

/********************************************************************/

real inp_debtemp(real deb_temp, real mass, real temp)

/*********************************************************************

  Calculate mean square displacement <dr^2> from Debye temperature,
  mass and temperature.

  INPUT

  real deb_temp  - Debye temperature in K
  real mass      - mass in amu
  real temp      - temperature in K

DESIGN

  Use equations (27/28, p.30 VHT) to calculate <r^2>:

    <r^2> = 1/2 * 9*T/(m*k_B*T_D^2)               for large T (> T_D)
    <r^2> = 1/2 * 9/(m*k_B*T_D) * ( 1/4 + 1.642 * (T/T_D)^2 )
                                                  for small T (< T_D)
    <r^2> = 1/2 * 9/(m*k_B*T_D) * sqrt(1/16 + (T/T_D)^2 )
                                                  for intermediate T
    (PREF_DEBWAL m= 9/(m*k_B) )
    The factor 1/2 (not in VHT) makes <r^2> the mean square displacement
    rather than the amplitude.

  FUNCTION CALLS

  RETURN VALUES

  <dr^2>
   -1. if failed (and EXIT_ON_ERROR is not defined)

*********************************************************************/
{

int i,j, iaux;                /* counter, dummy  variables */

real dr2;
real faux;                    /* dummy variable */


/*********************************************************************
 Check parameters
*********************************************************************/

 if (deb_temp <= 0.)
 {
#ifdef ERROR
   fprintf(STDERR,
           "*** error (inp_debtemp): wrong value for Debye temperature: %f\n",             deb_temp);
#endif
#ifdef EXIT_ON_ERROR
   exit(1);
#else
   return(-1.);
#endif
 }

 if (mass <= 0.)
 {
#ifdef ERROR
   fprintf(STDERR,
           "*** error (inp_debtemp): wrong value for mass: %f\n", mass);
#endif
#ifdef EXIT_ON_ERROR
   exit(1);
#else
   return(-1.);
#endif
 }

 if (temp < 0.)
 {
#ifdef ERROR
   fprintf(STDERR,
           "*** error (inp_debtemp): wrong value for temperature: %f\n", temp);
#endif
#ifdef EXIT_ON_ERROR
   exit(1);
#else
   return(-1.);
#endif
 }

/*********************************************************************
 calculate dr2
*********************************************************************/

#ifdef CONTROL
 fprintf(STDCTR, "(inp_debtemp): Debye = %.1f, Mass = %.1f, Temp. = %.1f\n",
         deb_temp, mass, temp);
#endif

 faux = temp / deb_temp;

 if( faux < 0.125 )       /* small T < T_D / 8 */
 {
   dr2 = 0.5 * PREF_DEBWAL / (mass * deb_temp) * (0.25 + 1.642 * faux * faux);
 }
 else if( faux > 8.0 )    /* large T  > T_D * 8 */
 {
   dr2 = 0.5 * PREF_DEBWAL / (mass * deb_temp) * faux;
 }
 else                     /* intermediate T */
 {
   dr2 = 0.5 * PREF_DEBWAL / (mass * deb_temp) * R_sqrt(0.0625 + faux*faux);
 }

#ifdef CONTROL
 fprintf(STDCTR, "(inp_debtemp): dr2 = %.3f dr1 = %.3f \n",
         dr2 * BOHR * BOHR, R_sqrt(dr2)*BOHR );
#endif

 return(dr2);
}  /* end of function inp_debtemp */
