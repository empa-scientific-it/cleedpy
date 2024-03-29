/*********************************************************************
  GH/01.03.95
  file contains functions:

  ms_lsum_ii            (23.08.94)
     Calculate the lattice sum Llm for a periodic plane of scatterers.

Changes:
 GH/23.08.94 - Creation

*********************************************************************/

#include <math.h>
#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif

#include <stdio.h>

#include "leed.h"

/*
#define CONTROL
*/
#define WARNING
#define ERROR

#define EXIT_ON_ERROR

#ifdef WARNING
#define WARN_LEVEL 1000
#endif

#ifdef NOT      /***************/
static mat Ylm = NULL;         /* contains spherical harmonics Y(0,0) */
static real * expm_r = NULL;   /* storage space */
static real * expm_i = NULL;   /* storage space */
static int l_max_y = -1;       /* maximum l for the above arrays */
#endif

/*======================================================================*/
/*======================================================================*/

mat ms_lsum_ii ( mat Llm, real k_r, real k_i, real *k_in, real *a,
                 int l_max, real epsilon )

/************************************************************************

 Calculate the lattice sum Llm for a periodic plane of scatterers.

 INPUT:

   mat Llm
   real k_r, k_i
   real k_in   Incident k-vector: k_in[1] = k_in_x, k_in[2] = k_in_y.
   real *a[5]  Basis vectors of the real 2-dim unit cell:
               a[1] = a1_x, a[2] = a2_x, a[3] = a1_y, a[4] = a2_y;
   int l_max
   real epsilon - defines the radius of summation:
               if epsilon < 1. : cut off value for the amplitude of the
                                 wave function.
               if epsilon >= 1.: radius.

 DESIGN:

   * General *

   Calculate the factor Llm which is to be multiplied with the Clebsh-Gordan
   coefficients in order to get the matrix elements of the Greens function
   Gii for a single periodic plane of scatterers.

     Llm =  (-1)^m * 4 PI * Ylm (cos(theta) = 0, phi = 0) *
            * sum(R) [ H(1)l(k*|R|) * exp( i(-kin*R + m*phi(R)) ) ]

     k =   k_r (length of the electron wave vector) +
         i*k_i (damping constant, must be > 0.).
     kin = k-vector of the incoming wave.

     H(1)l Hankel function of the first kind.
     R =   lattice vector (runs over all positions in the plane,
           i.e. Rz = 0.).

   Ylm(0,0) is zero for even (l+m), therefore also Llm is zero for these
   (l,m)-pairs.

   * Radius of the summation: *

   The radius up to which the lattice points are summed up is determined
   from the damping constant k_i and a lower boundary for the modulus
   of the Hankle function epsilon:

     r_m = - ln(epsilon) / k_i or epsilon (see above)

   r_max be the square of the maximum radius r_m, then for all lattice
   points within this radius the equation holds:

     r_max > (n1*a1_x + n2*a2_x)^2 + (n1*a1_y + n2*a2_y)^2 =
             n1^2 * f1 + 2n1*n2 * f12 + n2^2 * f2

   where

     f1  = (a1_x^2 + a1_y^2);
     f2  = (a2_x^2 + a2_y^2);
     f12 = (a1_x*a2_x + a1_y*a2_y).

   The solution of the quadratic equation for n2 is:

     n2_min = -n1*f12/f2 - R_sqrt( n1^2*(f12^2 - f2f1)/f2^2 + r_max/f2 )
     n2_max = -n1*f12/f2 - R_sqrt( n1^2*(f12^2 - f2f1)/f2^2 + r_max/f2 )

   The boundaries for n1 are given by the condition that the argument of
   the square root must be positive:

     n1^2   <  r_max*f2 / (f1f2 - f2^2)

 RETURN VALUES:

   NULL if failed (and EXIT_ON_ERROR is not defined)

   Llm (may be different from input parameter). The storage scheme for
   Llm is in the natural order:

   l      0  1  1  1  2  2  2  2  2  3  3  3  3  3  3  3  4  4 ...
   m      0 -1  0  1 -2 -1  0  1  2 -3 -2 -1  0  1  2  3 -4 -3 ...
   index  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 ...

   I.e. index(l,m) = l*(l+1) + m + 1. Note that, like usually for matrices,
   the first array element Llm[0] is not occupied.

*************************************************************************/
{
int l,m;                       /* quantum numbers l,m */
int off;                       /* offset in the arrays Ylm and Llm */
int i, iaux;

int n1, n1_max;                /* counters for lattice vectors */
int n2, n2_min, n2_max;

real f1, f2, f12;              /* factors used in determining n1, n2 */
real fa, fb, fc;

real r0_x, r0_y, r_x, r_y;     /* lattice vectors */
real r_max, r_abs;
real a1_x, a1_y, a2_x, a2_y;   /* basic lattice vectors */

real faux_r, faux_i;
real faux2_r, faux2_i;
real exp_phi_i, exp_phi_r;
real exp_ikr_i, exp_ikr_r;

real *expm_r, *expm_i;         /* storage space */

mat Ylm;                       /* contains spherical harmonics Y(0,0) */
mat Hl;                        /* Hankel function */

 Ylm = NULL;
 Hl = NULL;

/**************************************************************************
  Check arguments: k_i
**************************************************************************/

 if( k_i <= 0.)                /* no convergence */
 {
#ifdef ERROR
   fprintf(STDERR," *** error (ms_lsum_ii): damping too small: k_i = %.2e\n",
           k_i);
#endif

#ifndef EXIT_ON_ERROR
   return(NULL);
#endif

#ifdef EXIT_ON_ERROR
   exit(1);
#endif
 }

/*
  Allocate memory for Llm (and preset all Llm with zero).
*/
 iaux = (l_max + 1)*(l_max + 1);
 Llm = matalloc( Llm, iaux, 1, NUM_COMPLEX );

/*
 -Calculate the spherical harmonics Ylm(0,0), if not done already,
  and multiply with factor 4PI.
  Note that 4PI * Yl-m(0,0) = (-1)^m * 4PI * Ylm(0,0) is stored
  in the array Ylm and not Ylm(0,0).

 -Allocate storage space expm, if necessary.
*/


 expm_r = (real *)calloc ( l_max+1, sizeof(real) );
 expm_i = (real *)calloc ( l_max+1, sizeof(real) );

 Ylm = r_ylm(Ylm, 0., 0., l_max);

 for(l = 0, i = 1; l <= l_max; l ++)
 {
   faux_r = M1P(l)*4.*PI;
   for(m = -l; m <= l; m++, i++)
   {
/*
  Ylm(0,0) is purely real => only multiply Ylm->rel
*/
     Ylm->rel[i] *= faux_r;
     faux_r = -faux_r;
   }
 }

#ifdef NOT

 if (l_max > l_max_y)
 {
   Ylm = r_ylm(Ylm, 0., 0., l_max);
   l_max_y = l_max;

   for(l = 0, i = 1; l <= l_max; l ++)
   {
     faux_r = M1P(l)*4.*PI;
     for(m = -l; m <= l; m++, i++)
     {
   /*
     Ylm(0,0) is purely real => only multiply Ylm->rel
   */
       Ylm->rel[i] *= faux_r;
       faux_r = -faux_r;
     }
   }

   if (expm_r == NULL) expm_r = (real *)calloc (l_max+1, sizeof(real) );
   else         expm_r = (real *)realloc(expm_r, l_max+1*sizeof(real) );

   if (expm_i == NULL) expm_i = (real *)calloc (l_max+1, sizeof(real) );
   else         expm_i = (real *)realloc(expm_i, l_max+1*sizeof(real) );
 }

#endif

/*
  Some often used values
*/

 if (epsilon < 1.) r_max = - log(epsilon) / k_i;
 else              r_max = epsilon;

#ifdef WARNING
 if( r_max > WARN_LEVEL)    /* poor convergence */
   fprintf(STDWAR,
     "* warning (ms_lsum_ii): damping very weak: k_i = %.2e, eps = %.2e\n",
     k_i, epsilon);
#endif

 r_max *= r_max;

 a1_x = a[1]; a1_y = a[3];
 a2_x = a[2]; a2_y = a[4];

 f1  = a1_x*a1_x + a1_y*a1_y;
 f2  = a2_x*a2_x + a2_y*a2_y;
 f12 = a1_x*a2_x + a1_y*a2_y;

 fa  = - f12/f2;                        /* prefactor of n1 */
 fb  = (f12*f12 - f1*f2 ) /(f2*f2);     /* prefactor of n1^2 */
 fc =  r_max/ f2;                 /* constant under the root */

#ifdef CONTROL
 fprintf(STDCTR,
   "(ms_lsum_ii):  a1  = (%.3f,%.3f) A, a2  =  (%.3f,%.3f) A\n",
     a1_x*BOHR, a1_y*BOHR, a2_x*BOHR, a2_y*BOHR);
 fprintf(STDCTR,
   "               k_in = (%7.3f,%7.3f) A^-1\n", k_in[1]/BOHR, k_in[2]/BOHR);
 fprintf(STDCTR,
   "               eps = %7.5f, k_i = %7.4f A^-1, r_max = %7.3f A\n",
    epsilon, k_i/BOHR, R_sqrt(r_max)*BOHR);
#endif


/***********************
  Two cases:
   - k_in != 0,
   - k_in == 0. makes things much more easier.
***********************/

/***********************
  Fist: k_in == 0.
  (Alternatively:
   if (R_hypot(k_in[1], k_in[2]) < Tolerance) )
***********************/
 if( (k_in[1] == 0.) && (k_in[2] == 0.) )
 {
/*
 loop over lattice vectors a1 and a2:
 n1 >= 0: Only one half of the lattice. The other half is added within the
          m - loop: phi(-r) = phi(r) + PI.
*/

   n1_max = (int) R_sqrt(r_max * f2 / (f1*f2 - f12*f12) );

   for ( n1 = 0, r0_x = 0., r0_y = 0.;
         n1 <= n1_max;
         r0_x += a1_x, r0_y += a1_y, n1 ++ )
   {
     faux_r = R_sqrt( fb*n1*n1 + fc );
     n2_min = (int) (fa*n1 - faux_r);
     n2_max = (int) (fa*n1 + faux_r);


     for ( n2 = n2_min, r_x = r0_x + n2_min*a2_x, r_y = r0_y + n2_min*a2_y;
           n2 <= n2_max; r_x += a2_x, r_y += a2_y, n2 ++ )
     {
     /*
       The origin is not included in the summation.
     */
       if ((n1 == 0) && (n2 == 0)) break;

       r_abs = R_hypot(r_x, r_y);
       Hl = c_hank1 ( Hl, k_r*r_abs, k_i*r_abs, l_max);

       exp_phi_r =  r_x/r_abs;                        /*   cos(phi(r)) */
       exp_phi_i = -r_y/r_abs;                        /* - sin(phi(r)) */

/*
 loops over l and m:
*/
     /* First l = 0:  +/- r */

       expm_r[0] = 1.; expm_i[0] = 0.;

       Llm->rel[1] += 2*Hl->rel[1];
       Llm->iel[1] += 2*Hl->iel[1];

/*
        Now all other l,m:

  Summation over +/- r (phi / phi + PI) and +/- m

  for even m's != 0:
     exp (-im (phi+PI)) = exp (-im phi)
  => add contributions from the two half lattices for even m's

  for odd m's:
     exp (-im (phi+PI)) = -exp (-im phi)
  => cancellation of the two half lattices for odd m's

  Because Y(0,0) = 0. for odd (l+m), only even l's have nonzero elements
  in Llm.

  => skip summation for odd l's
*/
       for(l = 2; l <= l_max; l+= 2)
       {
         cri_mul(expm_r+l-1, expm_i+l-1,
                 expm_r[l-2], expm_i[l-2], exp_phi_r, exp_phi_i);
         cri_mul(expm_r+l, expm_i+l,
                 expm_r[l-1], expm_i[l-1], exp_phi_r, exp_phi_i);

         off = l*(l+1) + 1;

     /* First m = 0, exp(-im*phi) = 1. */
         Llm->rel[off] += 2*Hl->rel[l+1];
         Llm->iel[off] += 2*Hl->iel[l+1];

     /* Now all other m's (+ m and - m in the same loop) */
         for(m = 2; m <= l; m +=2 )
         {
         /* + m: */
           cri_mul(&faux_r, &faux_i,
                   expm_r[m], expm_i[m], Hl->rel[l+1], Hl->iel[l+1]);
           Llm->rel[off + m] += 2*faux_r;
           Llm->iel[off + m] += 2*faux_i;

         /* - m: exp(-i(-m) phi) = exp(-im phi)* */
           cri_mul(&faux_r, &faux_i,
                   expm_r[m], -expm_i[m], Hl->rel[l+1], Hl->iel[l+1]);
           Llm->rel[off - m] += 2*faux_r;
           Llm->iel[off - m] += 2*faux_i;
         }   /* m */
       }   /* l */
     }   /* lattice vectors a2 */
   }   /* lattice vectors a1 */
 }    /* end of k_in == 0. */

/***********************
  Now: k_in != 0.
***********************/

 else /* k_in != 0. */
 {
/*
 loop over lattice vectors a1 and a2:
 n1 >= 0: Only one half of the lattice. The other half is added within the
          m - loop: phi(-r) = phi(r) + PI.
*/

   Hl = NULL;
   n1_max = (int) R_sqrt(r_max * f2 / (f1*f2 - f12*f12) );

   for ( n1 = 0, r0_x = 0., r0_y = 0.;
         n1 <= n1_max;
         r0_x += a1_x, r0_y += a1_y, n1 ++ )
   {
     faux_r = R_sqrt( fb*n1*n1 + fc );
     n2_min = (int) (fa*n1 - faux_r);
     n2_max = (int) (fa*n1 + faux_r);

     for ( n2 = n2_min, r_x = r0_x + n2_min*a2_x, r_y = r0_y + n2_min*a2_y;
         n2 <= n2_max; r_x += a2_x, r_y += a2_y, n2 ++ )
     {
     /*
       The origin is not included in the summation.
     */
       if ((n1 == 0) && (n2 == 0)) break;

       r_abs = R_hypot(r_x, r_y);
       Hl = c_hank1 ( Hl, k_r*r_abs, k_i*r_abs, l_max);

       exp_phi_r =  r_x/r_abs;                        /*   cos(phi(r)) */
       exp_phi_i = -r_y/r_abs;                        /* - sin(phi(r)) */

       faux_r = k_in[1]*r_x + k_in[2]*r_y;
       cri_expi(&exp_ikr_r, &exp_ikr_i, faux_r, 0.); /* exp(i*k_in*r) */

/*
  loops over l and m:
*/
     /* first l = 0 */

       expm_r[0] = 1.; expm_i[0] = 0.;

       /* + r */
       cri_mul(&faux_r, &faux_i, exp_ikr_r, exp_ikr_i, Hl->rel[1], Hl->iel[1]);
       Llm->rel[1] += faux_r;
       Llm->iel[1] += faux_i;

       /* - r */
       cri_mul(&faux_r, &faux_i, exp_ikr_r,-exp_ikr_i, Hl->rel[1], Hl->iel[1]);
       Llm->rel[1] += faux_r;
       Llm->iel[1] += faux_i;

     /* Now all other l,m */
       for(l = 1; l <= l_max; l++ )
       {
         off = l*(l+1) + 1;
         cri_mul(expm_r+l, expm_i+l,
                 expm_r[l-1], expm_i[l-1], exp_phi_r, exp_phi_i);

/*
  Summation over +/- r (phi / phi + PI) and +/- m

  for even m's != 0:
     exp (-im (phi+PI)) = exp (-im phi)
     exp (ik(-r) )     = (exp (ikr) )*

  for odd m's:
     exp (-im (phi+PI)) = -exp (-im phi)
     exp (ik(-r) )     = (exp (ikr) )*

  => no cancellation.

  Because Y(0,0) = 0. for odd (l+m), only even l's have nonzero elements
  in Llm.

  => skip summation for odd (l+m)'s
*/
         if(ODD(l))
         {
       /*
         Loop over all odd m's (+ m and - m in the same loop)
       */
           for(m = 1; m <= l; m += 2)
           {
         /* + m: */
             cri_mul(&faux_r, &faux_i,
                     expm_r[m], expm_i[m], Hl->rel[l+1], Hl->iel[l+1]);
           /* + r */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r, exp_ikr_i, faux_r, faux_i);
             Llm->rel[off + m] += faux2_r;
             Llm->iel[off + m] += faux2_i;
           /* - r  '-' for odd m's */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r,-exp_ikr_i,-faux_r,-faux_i);
             Llm->rel[off + m] += faux2_r;
             Llm->iel[off + m] += faux2_i;

         /* - m: exp(-i(-m) phi) = exp(-im phi)* */
             cri_mul(&faux_r, &faux_i,
                     expm_r[m], -expm_i[m], Hl->rel[l+1], Hl->iel[l+1]);
           /* + r */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r, exp_ikr_i, faux_r, faux_i);
             Llm->rel[off - m] += faux2_r;
             Llm->iel[off - m] += faux2_i;
           /* - r  '-' for odd m's */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r,-exp_ikr_i,-faux_r,-faux_i);
             Llm->rel[off - m] += faux2_r;
             Llm->iel[off - m] += faux2_i;
           }   /* m */
         }      /* odd l's */
         else   /* even l's */
         {
       /*
         Loop over all even m's:
         - First m = 0, exp(-im*phi) = 1.
         - Then all other even m's.
       */

       /* First m = 0, exp(-im*phi) = 1. */
         /* + r */
           cri_mul(&faux_r, &faux_i,
                   exp_ikr_r, exp_ikr_i, Hl->rel[l+1], Hl->iel[l+1]);
           Llm->rel[off] += faux_r;
           Llm->iel[off] += faux_i;
         /* - r  */
           cri_mul(&faux_r, &faux_i,
                   exp_ikr_r,-exp_ikr_i, Hl->rel[l+1], Hl->iel[l+1]);
           Llm->rel[off] += faux_r;
           Llm->iel[off] += faux_i;

     /* Now all other even m's (+ m and - m in the same loop) */
           for(m = 2; m <= l; m += 2)
           {
         /* + m: */
             cri_mul(&faux_r, &faux_i,
                     expm_r[m], expm_i[m], Hl->rel[l+1], Hl->iel[l+1]);
           /* + r */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r, exp_ikr_i, faux_r, faux_i);
             Llm->rel[off + m] += faux2_r;
             Llm->iel[off + m] += faux2_i;
           /* - r  '+' for even m's */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r,-exp_ikr_i, faux_r, faux_i);
             Llm->rel[off + m] += faux2_r;
             Llm->iel[off + m] += faux2_i;

         /* - m: exp(-i(-m) phi) = exp(-im phi)* */
             cri_mul(&faux_r, &faux_i,
                     expm_r[m], -expm_i[m], Hl->rel[l+1], Hl->iel[l+1]);
           /* + r */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r, exp_ikr_i, faux_r, faux_i);
             Llm->rel[off - m] += faux2_r;
             Llm->iel[off - m] += faux2_i;
           /* - r  '+' for even m's */
             cri_mul(&faux2_r, &faux2_i, exp_ikr_r,-exp_ikr_i,faux_r,faux_i);
             Llm->rel[off - m] += faux2_r;
             Llm->iel[off - m] += faux2_i;
           }   /* even l's */
         }   /* m */
       }   /* l */
     }   /* lattice vectors a2 */
   }   /* lattice vectors a1 */
 }  /* end of k_in != 0 */


/***********************
 Final loop over (l,m) to multiply Llm with the factor
   (4PI * Yl-m)
 (Multiplication with i instead of Ylm gives the same results as
 FMAT in the VHT package)
***********************/
 iaux = (l_max+1)*(l_max+1);
 for(i = 1; i <= iaux; i ++)
 {
 /*
   Ylm(0,0) is purely real => only multiply with Ylm->rel
 */
    Llm->rel[i] *= Ylm->rel[i];
    Llm->iel[i] *= Ylm->rel[i];
 }

 matfree(Hl);
 matfree(Ylm);
 free(expm_r);
 free(expm_i);

 return(Llm);

} /* end of function ms_lsum_ii */

/*======================================================================*/
/*======================================================================*/
