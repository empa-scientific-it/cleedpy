/*********************************************************************
  GH/11.08.95
  file contains function:

  int out_bmlist(struct beam_str ** p_beams_out,
                 struct beam_str *beams_all,
                 sruct eng_str *eng,
                 FILE * outfile)

 Write header information to output file.

 Changes:

 GH/20.07.95 - Creation
 GH/11.08.95 - write only non-evanescent beams to output. Return value
               is a list of nonevanescent beams at eng->fin.

*********************************************************************/

#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif
#include <stdio.h>

#include "leed.h"


int out_bmlist(struct beam_str ** p_beams_out,
               struct beam_str *beams_all,
               struct eng_str *eng,
               real ** index1,
               real ** index2,
               int ** beam_set
               )

/************************************************************************

 write header information to output file.

 INPUT:

  struct beam_str ** p_beams_out  - (input) Pointer to output list (can
            be NULL to indicate that list should be allocated).
            The list will be terminated by "F_END_OF_LIST" in the
            structure element "k_par".

  struct beam_str *beams_all  - (input)
            List of all beams used throughout the energy
            loop. The list must be terminated by "F_END_OF_LIST" in the
            structure element "k_par".

            NOTE:

            All structure elements of this list will be copied if they
            are non-evanescent at energy eng->fin.
            It is assumed, that the elements beams_all->k_par contain the
            SQUARE OF THE PARALLEL MOMENTUM as created in function
            bm_gen.

  sruct eng_str *eng - (input) control parameters for energy loop:
            eng->ini = first energy,
            eng->fin = last energy,
            eng->stp = energy step.

  FILE * outfile - (input) pointer to the output file were the intensities
            are written to.

 DESIGN:

 *Create list beams_out*

  First find the maximum pareallel vector component of non-evanescent wave
  vectors at eng->fin.

  k_max = sqrt(2 * eng->fin)

  All beams with parallel momentum < k_max (i.e. k_par < k_max ^2) are
  included in list beams_out.

  n_beams is then set to number of output beams.

 *Output to file*

  "#en" energies:
        initial energy, final energy, energy step
  "#bn" number of beams
  "#bi" beam indices:
        number (starting from zero), 1st index, 2nd index, beam set.

 RETURN VALUES:

  int n_beams  number of output beams.
  -1           if failed (not implemented).

*************************************************************************/
{

    int n_beams, n_eng;
    int i_bm_all, i_bm_out;

    real k_max;
    real faux;

    struct beam_str *beams_out;

    /************************************************************************
     Find number of beams in list beams_all (n_beams) and
    allocate *p_beams_out = beams_out of the same size.
    ************************************************************************/

    for(n_beams = 0; (beams_all + n_beams)->k_par != F_END_OF_LIST; n_beams ++);


    if (*p_beams_out == NULL)
        *p_beams_out = beams_out = (struct beam_str *) calloc(n_beams + 1, sizeof(struct beam_str));
    else
        *p_beams_out = beams_out = (struct beam_str *) realloc(*p_beams_out, (n_beams+1) * sizeof(struct beam_str));

    if(beams_out == NULL)
    {
        fprintf(STDERR," *** error (out_bmlist): allocation error.\n");
        exit(1);
    }

    /************************************************************************
     Write appropriate beams to beams_out

    k_max is the square of the maximum pareallel vector component of
            non-evanescent wave vectors at eng->fin.
    n_beams is set to number of output beams afterwards.
    ************************************************************************/

    k_max =  2. * eng->fin;

    for(i_bm_all = 0, i_bm_out = 0; i_bm_all < n_beams; i_bm_all ++)
    {
        if( (beams_all + i_bm_all)->k_par <= k_max )
        {
            memcpy(beams_out + i_bm_out, beams_all + i_bm_all, sizeof(struct beam_str) );
            i_bm_out ++;
        }
    }

    /* terminate list beams_out */
    (beams_out + i_bm_out)->k_par = F_END_OF_LIST;
    n_beams = i_bm_out;


    /************************************************************************
     Write energies, number of beams, and beams to output

    "#en" energies:
            number of energies, initial energy, final energy, energy step.
    "#bn" number of beams
    "#bi" beam indices:
            number (starting from zero), 1st index, 2nd index, beam set.

    - Count energy points (Running through the same loop0 as the main
        program does avoids inconsistencies due to rounding errors)
    -
    ************************************************************************/

    /* energies */
    for(faux = eng->ini, n_eng = 0; faux <= eng->fin; faux += eng->stp, n_eng++);


    // Allocate memory for index1 and index2
    (* index1) = (real *) calloc(n_beams, sizeof(real));
    (* index2) = (real *) calloc(n_beams, sizeof(real));
    (* beam_set) = (int *) calloc(n_beams, sizeof(int));


    /* beam indices */
    for(i_bm_out = 0; i_bm_out < n_beams; i_bm_out ++)
    {
        (* index1)[i_bm_out] = (beams_out+i_bm_out)->ind_1;
        (* index2)[i_bm_out] = (beams_out+i_bm_out)->ind_2;
        (* beam_set)[i_bm_out] = (beams_out+i_bm_out)->set;
    }


    /* write beams_out back to pointer */
    *p_beams_out = beams_out;

    return n_beams;
}
