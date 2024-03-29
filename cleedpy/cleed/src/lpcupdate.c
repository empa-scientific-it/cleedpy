/*********************************************************************
  GH/20.09.95
  file contains function:

  pc_update(struct var_str *v_par, struct phs_str *phs_shifts,
            real energy)

 Update all parameters, that change during the energy loop.

Changes:

GH/20.01.95 - include function pc_mktl (structure var_str has changed)
GH/20.09.95 - optional variable vi (structure var_str has changed).
            - use pc_mktl_nd

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


#define CONTROL

#define CONTROL
#define WARNING
#define ERROR

#ifndef VI_START            /* should be defined in "leed_def.h" */
#define VI_START (100./HART)
#endif


int pc_update(struct var_str *v_par,
                          struct phs_str *phs_shifts, real energy)

/************************************************************************

 Update all parameters, that change during the energy loop. The

 INPUT:

  struct var_str *v_par - all parameters that change during the
                energy loop (for details see "leed_def.h").
                The parameter structure must exist and must be preset already.
  struct phs_str *phs_shifts - phase shifts (will be handed to function
                pc_mktl_nd)
  real energy - new energy (vacuum energy)

 RETURN VALUES:

  struct var_str *v_par
  (The function returns its first argument.)

 DESIGN:

*Energy*

 Real part of energy (v_par->eng_r) is set to
  E = (vacuum energy) - (real part of opt. potential)

 Imag. part of energy (v_par->eng_i is set to imag. part of opt. potential
 which is defined as:

 v_par->eng_r <  VI_START : prefactor
 v_par->eng_r >= VI_START : prefactor * (E / VI_START) ^ (exponent)

 prefactor = v_par->vi_pre,
 exponent  = v_par->vi_exp.

*k_in*

 |k_in|   = sin(theta_in) * sqrt( 2*(vacuum energy) )
  k_in(x) = cos(phi_in) * |k_in|
  k_in(y) = sin(phi_in) * |k_in|

 FUNCTION CALLS:
  pc_mktl_nd

*************************************************************************/
{
real faux_r, faux_i;

/*********************************************************
  Set new energy
*********************************************************/

 v_par->eng_v = energy;
 v_par->eng_r = energy - v_par->vr;

 if (v_par->eng_r < VI_START)
   v_par->eng_i = v_par->vi_pre;
 else
 {
   faux_r = R_log(v_par->eng_r / VI_START) * v_par->vi_exp;
   v_par->eng_i = v_par->vi_pre * R_exp(faux_r);
 }

/*********************************************************
  Determine k_in,
*********************************************************/

 faux_r = R_sin(v_par->theta) * R_sqrt(2*v_par->eng_v);
 v_par->k_in[0] = faux_r;
 v_par->k_in[1] = faux_r * R_cos(v_par->phi);
 v_par->k_in[2] = faux_r * R_sin(v_par->phi);


#ifdef CONTROL
 fprintf(STDCTR,
  "(pc_update): new energy: Evac = %.2f; (Er, Ei) = (%.2f, %.2f) eV\n",
  v_par->eng_v*HART, v_par->eng_r*HART, v_par->eng_i*HART);
 fprintf(STDCTR,
  "             k_in = (%.3f, %.3f) A-1\n",
  v_par->k_in[1]/BOHR, v_par->k_in[2]/BOHR);
#endif
#ifdef CONTROL_X
 fprintf(STDCTR,"(pc_update): k_in = \t(%.2f, %.2f)\n",
                  v_par->k_in[1], v_par->k_in[2]);
#endif

/*********************************************************
  Update phase shifts (pc_mktl_nd)
*********************************************************/

 v_par->p_tl = pc_mktl_nd(v_par->p_tl, phs_shifts, v_par->l_max, v_par->eng_r);

 return(1);
}  /* end of function pc_update */
