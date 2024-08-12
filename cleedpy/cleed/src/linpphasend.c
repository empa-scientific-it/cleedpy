/**********************************************************************
  GH/15.07.03
  file contains function:

Changes:

  GH/04.07.94 - Creation
  GH/19.01.95 - add eng_max and eng_min
  GH/08.08.95 - update i_phase
  GH/02.09.97 - Set input path by environment variable CLEED_PHASE
                (makes PHASE_PATH obsolete)
  WB/26.02.98 - Correct control output
  GH/03.05.00 - extend list of parameters for inp_phase: t_type
  GH/15.07.03 - fix bug in energy scaleing factor for Ry (was 2./HART, now 2.).

*********************************************************************/

#include <math.h>
#if !defined (__MACH__)
  #include <malloc.h>
#endif
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include "leed.h"
#include "leed_def.h"


int inp_phase_nd( char * filename, real * dr, int t_type, struct phs_str **p_phs_shifts, int *i_phase)
/*********************************************************************

    Read phase shifts from an input file and store them.

    INPUT:

        char * filename (input) path to the phase shift file
        real * dr (input) displacement vector for thermic vibrations

        struct phs_str **p_phs_shifts (output) phase shifts.

    DESIGN:

        The phase shifts in the input file must be for increasing energies.
        The storage scheme is:

*********************************************************************/
{
    FILE *inp_stream;
    char linebuffer[STRSZ];
    char eng_type[STRSZ];

    struct phs_str *phs_shifts;

    int i, i_str, i_eng;
    int neng, lmax, nl;          /* neng = No. of energies to be read
                                    lmax = max. quantum number;
                                    nl   = No. of phase shifts = lmax + 1
                                */

    real eng_scale;
    real faux;

    if((*i_phase) > 0)
    {
        /* Compare filename, dr, and t_type with previous phaseshifts. Return the
        corresponding phase shift number if the same combination has already been read. */
        for(i=0; i< (*i_phase); i++)
            if( (!strcmp( (*p_phs_shifts + i)->input_file, filename) )         &&
                ( R_fabs(dr[1] - (*p_phs_shifts + i)->dr[1]) < GEO_TOLERANCE ) &&
                ( R_fabs(dr[2] - (*p_phs_shifts + i)->dr[2]) < GEO_TOLERANCE ) &&
                ( R_fabs(dr[3] - (*p_phs_shifts + i)->dr[3]) < GEO_TOLERANCE ) &&
                ( t_type == (*p_phs_shifts + i)->t_type )
            )
            {
                return(i);
                break;
            }
        (*i_phase) ++;
        *p_phs_shifts = (struct phs_str *)realloc(*p_phs_shifts, ((*i_phase) + 1) * sizeof(struct phs_str) );
    }
    else
    {
        (*i_phase) ++;
        *p_phs_shifts = (struct phs_str *) malloc( 2 * sizeof(struct phs_str) );
    }

    // Terminate list of phase shifts.

    (*(p_phs_shifts) + (*i_phase))->lmax = I_END_OF_LIST;


    phs_shifts = *(p_phs_shifts) + (*i_phase)-1;

    // Write dr and t_type to phs_shifts.
    for(i=0; i<=3; i++) phs_shifts->dr[i] = dr[i];
    phs_shifts->t_type = t_type;

    /********************************************************************
     Open and Read input file for a new set of phase shifts
    ********************************************************************/

    // Open input file. Copy the full filename into phs_shifts->input_file.
    if( (inp_stream = fopen(filename, "r")) == NULL)
    {
        fprintf(STDERR, " *** error (inp_phase): could not open file \"%s\"\n",filename);
        exit(1);
    }

    phs_shifts->input_file = strdup(filename);

    /********************************************************************
    Read the first line of the input file which contains the number of
    energies to be read in (neng) and the maximum phase shift quantum
    number (lmax).
    ********************************************************************/

    // Read comment lines
    while( *fgets(linebuffer, STRSZ, inp_stream) == '#');

    if ( linebuffer == NULL)     /* EOF found */
    {
        fprintf(STDERR, " *** error (inp_phase): unexpected EOF found while reading file \"%s\"\n", filename);
        exit(1);
    }
    else if( sscanf(linebuffer, "%d %d %s", &neng, &lmax, eng_type) < 2)
    {
        fprintf(STDERR, " *** error (inp_phase): improper input line in file \"%s\":\n%s", filename, linebuffer);
        exit(1);
    }

    // Define energy scale according to eng_type. The default is input in Hartree units (27.18 eV).
    if( !strncmp(eng_type,"eV",2) || !strncmp(eng_type,"EV",2) )
    {
        eng_scale = 1./HART;
    }
    else if( !strncmp(eng_type,"Ry",2) || !strncmp(eng_type,"RY",2) )
    {
        eng_scale = 2.;
    }
    else
    {
        eng_scale = 1.;
    }


    /********************************************************************
     Read energies and phase shifts.
    Find max and min energy
    NB: The search for blank or '-' after reading each number is needed
    because of the FORTRAN format used for the VHT input files which does
    not have any blank character between negative numbers.
    ********************************************************************/

    phs_shifts->lmax = lmax;
    nl = lmax + 1;

    phs_shifts->energy = (real *)calloc( neng, sizeof(real) );
    phs_shifts->pshift = (real *)calloc( neng * nl, sizeof(real) );

    for( i_eng = 0; (i_eng < neng) && (fgets(linebuffer, STRSZ, inp_stream) != NULL); i_eng ++)
    {
        sscanf(linebuffer, "%le", phs_shifts->energy+i_eng);
        phs_shifts->energy[i_eng] *= eng_scale;
        if (i_eng == 0)
            phs_shifts->eng_min = phs_shifts->energy[i_eng];
        else
            phs_shifts->eng_max = phs_shifts->energy[i_eng];

        if( fgets(linebuffer, STRSZ, inp_stream) != NULL)
        {
            for( i_str = 0, i = 0; i<nl; i++)
            {
                sscanf(linebuffer + i_str, "%le", phs_shifts->pshift+i_eng*nl+i);
                while((linebuffer[i_str] == ' ') || (linebuffer[i_str] == '-')) i_str ++;
                while((linebuffer[i_str] != ' ') && (linebuffer[i_str] != '-')) i_str ++;
            }
        }
        else
        {
            phs_shifts->energy[i_eng] = 0.;
            phs_shifts->eng_max = phs_shifts->energy[i_eng-1];
            break;
        }
    }
    phs_shifts->neng = i_eng;

    if(phs_shifts->neng != neng)
    {
        fprintf(STDWAR, " *** warning (inp_phase): EOF found before reading all phase shifts:\n");
        fprintf(STDWAR, "     expected energies: %3d, found: %3d, file: %s\n", neng, i_eng+1, filename);
    }

    return((*i_phase) - 1);
} /* end of function inp_phase */
