/*********************************************************************
  GH/15.08.94
  file contains functions:

  ms_ymat            (15.08.94)
     Create the transformation matrix from angular momentum space
     into k-space: Ylm(k).

*********************************************************************/

#include <math.h>
#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif
#include <stdio.h>
#include <string.h>

#include "leed.h"

/*
#define CONTROL
#define WARNING
*/
#define ERROR

#define EXIT_ON_ERROR

/*======================================================================*/
/*======================================================================*/

static mat Ylm = NULL;

mat ms_ymat ( mat Ymat, int l_max, struct beam_str *beams, int n_beams)

/************************************************************************

 Create the transformation matrix Ymat from angular momentum space into
 k-space: Ylm(k).

 INPUT:

   mat Ymat - (output) transformation matrix.
   int l_max - max l quantum number spanning up the (l,m)-space.
   struct beam_str *beams - beams spanning up the k-space.
   n_beams - number of beams in list beams.

 RETURN VALUES:

   NULL if failed (and EXIT_ON_ERROR is not defined)

   Ymat (may be different from input parameter).

   The order of Ymat is:

   k / (l,m) (0,0)     (-1,1)    (0,1)   (1,1)   (-2,2) ....

   k1       Y00(k1)   Y-11(k1)  Y01(k1)
   k2       Y00(k2)   Y-11(k2)  Y01(k2)
   k3       Y00(k3)
   ...                                     Y(l,m)(ki)



*************************************************************************/
{
int l,m;
int ll_max;                    /* total number of (l,m) pairs */
int off;                       /* offset in the arrays Ylm and Llm */
int size;
int i_beams;

real faux_r, faux_i;
real *ptr_r, *ptr_i, *ptr_end;

/*
  Check arguments:
*/

/*
  Allocate memory for Ymat
*/
 ll_max = (l_max + 1)*(l_max + 1);
 Ymat = matalloc( Ymat, n_beams, ll_max, NUM_COMPLEX );

/*
  Loop over k's:
  Calculate the sperical harmonics Ylm(k) and copy into Ymat
*/

 size = ll_max*sizeof(real);

 for (i_beams = 0, off = 1; i_beams < n_beams; i_beams ++, off += ll_max)
 {

#ifdef CONTROL
   fprintf(STDCTR,"(ms_ymat): cos(th): (%.3f %.3f), phi: %.3f\n",
           (beams+i_beams)->cth_r,
           (beams+i_beams)->cth_i, (beams+i_beams)->phi);
#endif

   Ylm = c_ylm(Ylm, (beams+i_beams)->cth_r, (beams+i_beams)->cth_i,
                    (beams+i_beams)->phi, l_max);
   memcpy( Ymat->rel+off, Ylm->rel+1, size);
   memcpy( Ymat->iel+off, Ylm->iel+1, size);
 }

 return(Ymat);
}

/*======================================================================*/
/*======================================================================*/
