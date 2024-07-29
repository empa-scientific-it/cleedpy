/*********************************************************************
  GH/11.08.95
  file contains function:

  out_int(mat Amp, struct beam_str *beams, struct var_str *par, FILE * outfile)

 Intensity output function

 Changes:

 GH/20.01.95 - Creation
 GH/11.08.95 - Minor changes

*********************************************************************/

#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif
#include <stdio.h>

#include "leed.h"

/*
#define CONTROL_ALL
#define WARNING
*/
//#define CONTROL
#define ERROR

#define EXIT_ON_ERROR

int out_int(mat Amp, struct beam_str *beams_now, struct beam_str *beams_all,
            struct var_str *par, FILE * outfile, real * beam_intensities)

/************************************************************************

 output of beam intensities

 INPUT:

  mat Amp - (input) vector containing the beam amplitudes of all beams
           included at the current energy.
  struct beam_str *beams_now  -  all beams included at the current energy.
  struct beam_str *beams_all  -  all beams included at the highest energy.
  FILE * outfile - pointer to the output file were the intensities are
           written to.

 RETURN VALUES:

  number of intensities written to outfile.
  -1 if failed.

 DESIGN:

*************************************************************************/
{

int i_beams_now, i_beams_all, i_out;

mat Int;
real k_r;

 Int = NULL;
/*********************************************************
   Calculate intensitied as the square of the moduli of the
   amplitudes.
*********************************************************/

 Int = matsqmod(Int, Amp);

 k_r = R_sqrt(2*par->eng_v);

/*********************************************************
   Print intensities for non-evanescent beams
*********************************************************/
 fprintf(outfile,"%.2f ", par->eng_v*HART);

 for(i_beams_all = 0; (beams_all + i_beams_all)->k_par != F_END_OF_LIST;
     i_beams_all ++)
 {
   for(i_beams_now = 0; i_beams_now < Int->rows; i_beams_now ++)
   {
     if( ((beams_all+i_beams_all)->ind_1 == (beams_now+i_beams_now)->ind_1) &&
         ((beams_all+i_beams_all)->ind_2 == (beams_now+i_beams_now)->ind_2)  )
     {
        if (((beams_now+i_beams_now)->k_par <= k_r) && (Int->rel[i_beams_now + 1] > INT_TOLERANCE))
        {
            fprintf(outfile,"%.6e ", Int->rel[i_beams_now + 1]);
            beam_intensities[i_beams_all] = Int->rel[i_beams_now + 1];
        }
        else
        {
            fprintf(outfile,"%.6e ",0.);
            beam_intensities[i_beams_all] = 0.;
        }
       break;
     }  /* if index = index */
   }  /* for i_beams_now */

   if( i_beams_now == Int->rows) fprintf(outfile,"%.6e ",0.);
 } /* for i_beams_all */

 fprintf(outfile,"\n");

 fflush(outfile);
 return(i_beams_now);
}  /* end of function out_int */
