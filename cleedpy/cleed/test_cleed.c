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

int energy_list_size, energy_index;
real *energy_list;

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

  inp_showbop(bulk, over, phs_shifts);

  if( ctr_flag == CTR_EARLY_RETURN )
  {
    fprintf(STDCTR, "(LEED_TEMP): EARLY RETURN \n");
    exit(0);
  }

  out_head (bulk, res_stream);

/*********************************************************************
 Prepare some often used parameters.
*********************************************************************/

  mk_cg_coef (2*v_par->l_max);
  mk_ylm_coef(2*v_par->l_max);

  // Construct energy list
  energy_list_size = (eng->fin - eng->ini)/eng->stp + 1;
  energy_list = (real *) malloc(energy_list_size * sizeof(real));
  for (energy_index=0; energy_index<energy_list_size; energy_index++)
  {
    energy_list[energy_index] = eng->ini + energy_index * eng->stp;
  }

  leed(bulk, over, phs_shifts, energy_list_size, energy_list, v_par, res_stream);
  printf("Finished leed\n");

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

  printf("The new code is working\n");
  exit(0);

} /* end of main */
