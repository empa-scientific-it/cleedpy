/*********************************************************************
GH/05.10.00
  file contains functions:

  main
     Main program for LEED calculations (only for bravaislayer)

Changes:

GH/05.10.00 - early return option (-e).


*********************************************************************/

#include <stdio.h>
#include <string.h>

#include "gh_stddef.h"
#include "leed.h"



/*
#define CONTROL_X
*/
#define CONTROL_FLOW


#define CONTROL_IO
#define CONTROL
#define WARNING
#define ERROR

#define CTR_NORMAL       998
#define CTR_EARLY_RETURN 999

/*======================================================================*/

main(int argc, char *argv[])

/*********************************************************************
 Perform a LEED calculation for anisotropic vibrations a general case

 INPUT:

 DESIGN:

*********************************************************************/
{
struct cryst_str *bulk;
struct cryst_str *over;
struct phs_str *phs_shifts;
struct beam_str *beams_all;
struct beam_str *beams_out;
struct beam_str *beams_now;
struct beam_str *beams_set;
struct var_str *v_par;
struct eng_str *eng;

mat Tpp,   Tmm,   Rpm,   Rmp;
mat Tpp_s, Tmm_s, Rpm_s, Rmp_s;
mat R_bulk, R_tot;
mat Amp;

int ctr_flag;
int i_c, i_arg;
int n_beams_now, n_beams_set;
int i_set, n_set, offset;
int i_layer;

real energy;
real vec[4];

char linebuffer[STRSZ];

char bul_file[STRSZ];                 /* input/output files */
char par_file[STRSZ];
char pro_name[STRSZ];
char res_file[STRSZ];

FILE *pro_stream;
FILE *res_stream;

  Tpp   =  Tmm   =  Rpm   =  Rmp   = NULL;
  Tpp_s =  Tmm_s =  Rpm_s =  Rmp_s = NULL;
  R_bulk = R_tot = NULL;
  Amp = NULL;

  bulk = over = NULL;
  phs_shifts  = NULL;
  beams_all   = NULL;
  beams_now   = NULL;
  beams_set   = NULL;
  beams_out   = NULL;
  v_par = NULL;
  eng   = NULL;


/*********************************************************************
  Preset parameters (file names) set by arguments
*********************************************************************/

  ctr_flag = CTR_NORMAL;

  strncpy(bul_file,"---", STRSZ);

  strncpy(par_file,"---", STRSZ);
  strncpy(res_file,"leed.res", STRSZ);
  strncpy(pro_name,"leed.pro", STRSZ);

/*********************************************************************
  Decode arguments:

    -b <bul_file> - (optional input file) bulk and non-geometrical
                    parameters.
    -i <par_file> - (mandatory input file) overlayer parameters of all
                    parameters (if bul_file does not exist).
    -o <res_file> - (output file) IV output.
*********************************************************************/

  for (i_arg = 1; i_arg < argc; i_arg++)
  {
    if(*argv[i_arg] != '-')
    {
#ifdef ERROR
      fprintf(STDERR,"*** error (LEED_TEMP):\tsyntax error:\n");
      fprintf(STDERR,"\tusage: \tleed -i <par_file> -o <res_file>");
      fprintf(STDERR," [-b <bul_file> -e]\n");
#endif
      exit(1);
    }
    else
    {

/* Read parameter input file */
      if(strncmp(argv[i_arg], "-b", 2) == 0)
      {
        i_arg++;
        strncpy(bul_file, argv[i_arg], STRSZ);
      } /* -b */

/* Read parameter input file */
      if(strncmp(argv[i_arg], "-i", 2) == 0)
      {
        i_arg++;
        strncpy(par_file, argv[i_arg], STRSZ);
      } /* -i */

/* Read and open results file */
      if(strncmp(argv[i_arg], "-o", 2) == 0)
      {
        i_arg++;
        strncpy(res_file, argv[i_arg], STRSZ);
        if ((res_stream = fopen(res_file,"w")) == NULL)
        {
#ifdef ERROR
          fprintf(STDERR,
          "*** error (LEED_TEMP): could not open output file \"%s\"\n",
          res_file);
#endif
          exit(1);
        }
      }  /* -o */

/* Read parameter input file */
      if(strncmp(argv[i_arg], "-e", 2) == 0)
      {
        ctr_flag = CTR_EARLY_RETURN;
      } /* -e */


    }  /* else */
  }  /* for i_arg */

/*********************************************************************
  Check arguments:
  - check existence of par_file.
  - if bul_file is not specified, use par_file instead.
  - check existence of res_file.
*********************************************************************/

  if(strncmp(par_file, "---", 3) == 0)
  {
#ifdef ERROR
    fprintf(STDERR,
    "*** error (LEED_TEMP): no parameter input file (option -i) specified\n");
#endif
    exit(1);
  }

  if(strncmp(bul_file, "---", 3) == 0)
  {
    strncpy(bul_file, par_file, STRSZ);
  }

  if(strncmp(res_file, "leed.res", 8) == 0)
  {
#ifdef WARNING
    fprintf(STDWAR,
            "* warning (LEED_TEMP): no output file (option -o) specified\n");
    fprintf(STDWAR,"\toutput will be written to file \"%s\"\n", res_file);
#endif

    if ((res_stream = fopen(res_file,"w")) == NULL)
    {
#ifdef ERROR
      fprintf(STDERR,
      "*** error (LEED_TEMP): could not open output file \"%s\"\n",
      res_file);
#endif
     exit(1);
    }
  }

/*********************************************************************
  Read input parameters
*********************************************************************/

  inp_rdbul_nd(&bulk, &phs_shifts, bul_file);
  inp_rdpar(&v_par, &eng, bulk, bul_file);
  inp_rdovl_nd(&over, &phs_shifts, bulk, par_file);
  n_set = bm_gen(&beams_all, bulk, v_par, eng->fin);

  inp_showbop(bulk, over, phs_shifts);

  if( ctr_flag == CTR_EARLY_RETURN )
  {
    fprintf(STDCTR, "(LEED_TEMP): EARLY RETURN \n");
    exit(0);
  }

  out_head (bulk, res_stream);
  out_bmlist(&beams_out, beams_all, eng, res_stream);

/*********************************************************************
 Prepare some often used parameters.
*********************************************************************/

  mk_cg_coef (2*v_par->l_max);
  mk_ylm_coef(2*v_par->l_max);

#ifdef CONTROL
  fprintf(STDCTR, "(LEED_TEMP): E_ini = %.1f, E_fin = %.1f, E_stp %.1f\n",
          eng->ini*HART, eng->fin*HART, eng->stp*HART);

  fprintf(STDCTR, "(LEED_TEMP): n_set = %d\n", n_set);
#endif

/*********************************************************************
 Energy Loop
*********************************************************************/

  for( energy = eng->ini; energy < eng->fin + E_TOLERANCE; energy += eng->stp)
  {
    pc_update(v_par, phs_shifts, energy);
    n_beams_now = bm_select(&beams_now, beams_all, v_par, bulk->dmin);

#ifdef CONTROL
    fprintf(STDCTR, "(LEED_TEMP):\n\t => E = %.1f eV (%d beams used) <=\n\n",
              v_par->eng_v*HART, n_beams_now);
#endif


/*********************************************************************
  BULK:
  Loop over beam sets

  Create matrix R_bulk that will eventually contain the bulk
  reflection matrix
*********************************************************************/
    R_bulk = matalloc(R_bulk, n_beams_now, n_beams_now, NUM_COMPLEX);
    for(offset = 1, i_set = 0; i_set < n_set; i_set ++)
    {
      n_beams_set = bm_set(&beams_set, beams_now, i_set);

/*********************************************************************
    Loop over periodic bulk layers
*********************************************************************/

    /**********************************************************
     Compute scattering matrices for bottom-most bulk layer:
     - single Bravais layer or composite layer
    **********************************************************/

#ifdef CONTROL_FLOW
      fprintf(STDCTR, "(LEED_TEMP periodic): bulk layer %d/%d, set %d/%d\n",
                        0, bulk->nlayers - 1, i_set, n_set - 1);
#endif

      if( (bulk->layers + 0)->natoms == 1)
      {
        ms_bravl_nd( &Tpp, &Tmm, &Rpm, &Rmp,
                     v_par, (bulk->layers + 0), beams_set);
      }
      else
      {
        ms_compl_nd( &Tpp, &Tmm, &Rpm, &Rmp,
                     v_par, (bulk->layers + 0), beams_set);
      }

#ifdef CONTROL_X
      fprintf(STDCTR, "(LEED_TEMP): after ms_bravl_nd: Tpp:");
      matshow(Tpp);
      fprintf(STDCTR, "(LEED_TEMP): after ms_bravl_nd: Tmm:");
      matshow(Tmm);
      fprintf(STDCTR, "(LEED_TEMP): after ms_bravl_nd: Rpm:");
      matshow(Rpm);
      fprintf(STDCTR, "(LEED_TEMP): after ms_bravl_nd: Rmp:");
      matshow(Rmp);
#endif

    /**********************************************************
      Loop over the other bulk layers
    **********************************************************/

      for(i_layer = 1;
          ( (bulk->layers+i_layer)->periodic == 1) &&
          (i_layer < bulk->nlayers);
          i_layer ++)
      {
#ifdef CONTROL_FLOW
        fprintf(STDCTR, "(LEED_TEMP periodic): bulk layer %d/%d, set %d/%d\n",
                        i_layer, bulk->nlayers - 1, i_set, n_set - 1);
#endif

    /**************************************************************
      Compute scattering matrices R/T_s for a single bulk layer
       - single Bravais layer or composite layer
    ***************************************************************/

        if( (bulk->layers + i_layer)->natoms == 1)
        {
          ms_bravl_nd ( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                        v_par, (bulk->layers + i_layer), beams_set);
        }
        else
        {
          ms_compl_nd( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                       v_par, (bulk->layers + i_layer), beams_set);
        }

    /***************************************************************************
       Add the single layer matrices to the rest by layer doubling
       - inter layer vector is the vector between layers
         (i_layer - 1) and (i_layer):
         (bulk->layers + i_layer)->vec_from_last
    ****************************************************************************/
#ifdef CONTROL_FLOW
        fprintf(STDCTR,
                "(LEED_TEMP): before ld_2lay vec_from...(%.2f %.2f %.2f)\n",
                       (bulk->layers + i_layer)->vec_from_last[1] * BOHR,
                       (bulk->layers + i_layer)->vec_from_last[2] * BOHR,
                       (bulk->layers + i_layer)->vec_from_last[3] * BOHR);
#endif

        ld_2lay( &Tpp,  &Tmm,  &Rpm,  &Rmp,
                 Tpp,   Tmm,   Rpm,   Rmp,
                 Tpp_s, Tmm_s, Rpm_s, Rmp_s,
                 beams_set, (bulk->layers + i_layer)->vec_from_last);

      } /* for i_layer (bulk) */

   /*********************************************************************
      Layer doubling for all periodic bulk layers until convergence is
      reached:
       - inter layer vector is (bulk->layers + 0)->vec_from_last
   **********************************************************************/
#ifdef CONTROL_FLOW
      fprintf(STDCTR, "(LEED_TEMP): before ld_2n vec_from...(%.2f %.2f %.2f)\n",
                      (bulk->layers + 0)->vec_from_last[1] * BOHR,
                      (bulk->layers + 0)->vec_from_last[2] * BOHR,
                      (bulk->layers + 0)->vec_from_last[3] * BOHR);
#endif

      Rpm = ld_2n( Rpm, Tpp, Tmm, Rpm, Rmp,
                   beams_set, (bulk->layers + 0)->vec_from_last);

   /*******************************************************************
     Compute scattering matrices for top-most bulk layer if it is
     not periodic.
      - single Bravais layer or composite layer
   **********************************************************************/

      if( i_layer == bulk->nlayers - 1 )
      {
#ifdef CONTROL_FLOW
        fprintf(STDCTR,
                "(LEED_TEMP not periodic): bulk layer %d/%d, set %d/%d\n",
                i_layer, bulk->nlayers - 1, i_set, n_set - 1);
#endif

        if( (bulk->layers + i_layer)->natoms == 1)
        {
          ms_bravl_nd( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                    v_par, (bulk->layers + i_layer), beams_set);
        }
        else
        {
          ms_compl_nd( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                       v_par, (bulk->layers + i_layer), beams_set);
        }

    /**************************************************************************
       Add the single layer matrices of the top-most layer to the rest
       by layer doubling:
       - inter layer vector is the vector between layers
         (i_layer - 1) and (i_layer):
         (bulk->layers + i_layer)->vec_from_last
    ***************************************************************************/

        Rpm = ld_2lay_rpm(Rpm, Rpm, Tpp_s, Tmm_s, Rpm_s, Rmp_s,
                          beams_set, (bulk->layers + i_layer)->vec_from_last);

      }  /* if( i_layer == bulk->nlayers - 1 ) */

   /*******************************************************
     Insert reflection matrix for this beam set into R_bulk.
   ********************************************************/

      R_bulk = matins(R_bulk, Rpm, offset, offset);
      offset += n_beams_set;

   /*************************
     Write cpu time to output
   **************************/

      sprintf(linebuffer,"(LEED_TEMP): bulk layers set %d, E = %.1f",
              i_set, energy*HART);
      cpu_time(STDCPU,linebuffer);
    }  /* for i_set */

/*********************************************************************
  OVERLAYER
  Loop over all overlayer layers
*********************************************************************/

    for(i_layer = 0; i_layer < over->nlayers; i_layer ++)
    {
#ifdef CONTROL_FLOW
      fprintf(STDCTR, "(LEED_TEMP): overlayer %d/%d\n", i_layer, over->nlayers - 1);
#endif
   /***********************************************************
     Calculate scattering matrices for a single overlayer layer
      - only single Bravais layer
   ************************************************************/

      if( (over->layers + i_layer)->natoms == 1)
      {
        ms_bravl_nd( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                     v_par, (over->layers + i_layer), beams_now);
      }
      else
      {
        ms_compl_nd( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                     v_par, (over->layers + i_layer), beams_now);
      }

#ifdef CONTROL_X
 fprintf(STDCTR, "\n(LEED_TEMP):overlayer %d  ...\n",i_layer);
 fprintf(STDCTR, "\n(LEED_TEMP): Tpp:\n");
 matshowabs(Tpp_s);
 fprintf(STDCTR, "\n(LEED_TEMP): Tmm:\n");
 matshowabs(Tmm_s);
 fprintf(STDCTR, "\n(LEED_TEMP): Rpm:\n");
 matshowabs(Rpm_s);
 fprintf(STDCTR, "\n(LEED_TEMP): Rmp:\n");
 matshowabs(Rmp_s);
#endif

  /****************************************************************
     Add the single layer matrices to the rest by layer doubling:
     - if the current layer is the bottom-most (i_layer == 0),
       the inter layer vector is calculated from the vectors between
       top-most bulk layer and origin
       ( (bulk->layers + nlayers)->vec_to_next )
       and origin and bottom-most overlayer
       (over->layers + 0)->vec_from_last.

     - inter layer vector is the vector between layers
       (i_layer - 1) and (i_layer): (over->layers + i_layer)->vec_from_last
  **********************************************************************/

      if (i_layer == 0)
      {
        for(i_c = 1; i_c <= 3; i_c ++)
        {
          vec[i_c] = (bulk->layers + bulk->nlayers - 1)->vec_to_next[i_c]
                     + (over->layers + 0)->vec_from_last[i_c];
        }

#ifdef CONTROL_FLOW
        fprintf(STDCTR,
                "(LEED):over0 before ld_2lay_rpm vec..(%.2f %.2f %.2f)\n",
                vec[1] * BOHR,vec[2] * BOHR, vec[3] * BOHR);
#endif

        R_tot = ld_2lay_rpm(R_tot, R_bulk, Tpp_s, Tmm_s, Rpm_s, Rmp_s,
                            beams_now, vec);
      }
      else
      {
#ifdef CONTROL_FLOW
        fprintf(STDCTR,
                "(LEED):over%d  before ld_2lay_rpm vec..(%.2f %.2f %.2f)\n",
                i_layer,(over->layers + i_layer)->vec_from_last[1] * BOHR,
                        (over->layers + i_layer)->vec_from_last[2] * BOHR,
                        (over->layers + i_layer)->vec_from_last[3] * BOHR);
#endif

        R_tot = ld_2lay_rpm(R_tot, R_tot, Tpp_s, Tmm_s, Rpm_s, Rmp_s,
                          beams_now, (over->layers + i_layer)->vec_from_last);
      }

   /**************************
     Write cpu time to output
   **************************/

      sprintf(linebuffer,"(LEED): overlayer %d, E = %.1f",
              i_layer, energy * HART);
      cpu_time(STDCPU,linebuffer);

    }  /* for i_layer (overlayer) */

/*********************************************
   Add propagation towards the potential step.
**********************************************/

    vec[1] = vec[2] = 0.;
    vec[3] = 1.25 / BOHR;

/********************************************
    No scattering at pot. step
********************************************/

    Amp = ld_potstep0(Amp, R_tot, beams_now, v_par->eng_v, vec);
    out_int(Amp, beams_now, beams_out, v_par, res_stream);

/********************************************
    Write cpu time to output
********************************************/

    sprintf(linebuffer,"  %.1f   %d  ",energy * HART,n_beams_now);
    cpu_time(STDWAR,linebuffer);


  } /* end of energy loop */


#ifdef CONTROL_IO
  fprintf(STDCTR, "(LEED): end of energy loop: close files\n");
#endif

  fclose(res_stream);

#ifdef CONTROL
  fprintf(STDCTR, "\n\n(LEED):\tCORRECT TERMINATION");
#endif

/********************************************
    Write cpu time to output
********************************************/

  cpu_time(STDCPU,"");

/********************************************
    set exit status explicitly
********************************************/

  printf("The original code is working fine.\n");

  exit(0);

} /* end of main */
