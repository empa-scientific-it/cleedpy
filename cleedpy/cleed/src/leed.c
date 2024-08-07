#include "leed.h"

typedef struct {
    int n_beams;
    real * beam_index1;
    real * beam_index2;
    int * beam_set;
    int n_energies;
    real * energies;
    real * iv_curves;
} CleedResult;

void print_phase_shift(struct phs_str phs_shift)
{
    int i;
    printf("Phase shift:\n");
    printf("  lmax: %d\n", phs_shift.lmax);
    printf("  neng: %d\n", phs_shift.neng);
    printf("  t_type: %d\n", phs_shift.t_type);
    printf("  eng_max: %lf\n", phs_shift.eng_max);
    printf("  eng_min: %lf\n", phs_shift.eng_min);
    printf("  energy: %lf\n", phs_shift.energy[0]);
    printf("  dr: %lf\n", phs_shift.dr[0]);
    printf("  input_file: %s\n", phs_shift.input_file);
    printf("  pshift: ");
    for (i=0; i<phs_shift.neng; i++)
    {
        printf("%lf ", phs_shift.pshift[i]);
        if (i % 10 == 0)
            printf("\n");
    }
    printf("\n");

}


CleedResult leed(char * par_file, char * bul_file, char *phase_path)
{
    struct cryst_str *bulk=NULL;
    struct cryst_str *over=NULL;
    struct phs_str *phs_shifts=NULL;
    struct var_str *v_par=NULL;


    struct beams_str * beams_now=NULL;
    struct beams_str * beams_set=NULL;
    struct beam_str * beams_out=NULL;
    struct beam_str *beams_all=NULL;

    int n_beams_now, n_beams_set;

    CleedResult results;

    int i_c, i_set, offset, i;
    int i_layer;
    int energy_index;
    int n_set;
    real energy;
    real vec[4];
    mat R_bulk=NULL, R_tot=NULL;
    mat Amp=NULL;

    mat Tpp=NULL, Tmm=NULL, Rpm=NULL, Rmp=NULL;
    mat Tpp_s=NULL, Tmm_s=NULL, Rpm_s=NULL, Rmp_s=NULL;

    struct eng_str *eng=NULL;

    // Read input parameters
    inp_rdbul_nd(&bulk, &phs_shifts, bul_file, phase_path);
    inp_rdpar(&v_par, &eng, bulk, bul_file);
    inp_rdovl_nd(&over, &phs_shifts, bulk, par_file, phase_path);
    inp_showbop(bulk, over, phs_shifts);
    for (i=0; (phs_shifts + i)->lmax != I_END_OF_LIST; i++)
        print_phase_shift(phs_shifts[i]);



    // Construct energy list
    results.n_energies = (eng->fin - eng->ini)/eng->stp + 1;
    results.energies = (real *) malloc(results.n_energies * sizeof(real));

    for (energy_index=0; energy_index < results.n_energies; energy_index++)
        results.energies[energy_index] = eng->ini + energy_index * eng->stp;

    eng->fin = results.energies[results.n_energies - 1];

    // Printing stuff
    inp_showbop(bulk, over, phs_shifts);

    mk_cg_coef (2*v_par->l_max); // Setting up Clebsh Gordan coefficients as global variables.
    mk_ylm_coef(2*v_par->l_max); // Setting up spherical harmonics coefficients as global variables.




    /* Generate beams out */
    n_set = bm_gen(&beams_all, bulk, v_par, eng->fin);
    results.n_beams = out_bmlist(&beams_out, beams_all, eng, &results.beam_index1, &results.beam_index2, &results.beam_set);
    results.iv_curves = (real *) calloc(results.n_energies * results.n_beams, sizeof(real));

    /* Main Energy Loop */

    for(energy_index=0; energy_index < results.n_energies; energy_index++){
        pc_update(v_par, phs_shifts, results.energies[energy_index]);
        n_beams_now = bm_select(&beams_now, beams_all, v_par, bulk->dmin);

        /*********************************************************************
        BULK:
        Loop over beam sets

        Create matrix R_bulk that will eventually contain the bulk
        reflection matrix
        *********************************************************************/

        R_bulk = matalloc(R_bulk, n_beams_now, n_beams_now, NUM_COMPLEX);

        /*********************************************************************
            Loop over periodic bulk layers
        *********************************************************************/
        for(offset = 1, i_set = 0; i_set < n_set; i_set ++)
        {
            n_beams_set = bm_set(&beams_set, beams_now, i_set);

            /**********************************************************
             Compute scattering matrices for bottom-most bulk layer:
            - single Bravais layer or composite layer
            **********************************************************/
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

            /**********************************************************
            Loop over the other bulk layers
            **********************************************************/

            for(i_layer = 1;
            ((bulk->layers+i_layer)->periodic == 1) && (i_layer < bulk->nlayers);
            i_layer ++)
            {
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
            Rpm = ld_2n( Rpm, Tpp, Tmm, Rpm, Rmp,
                        beams_set, (bulk->layers + 0)->vec_from_last);

            /*******************************************************************
            Compute scattering matrices for top-most bulk layer if it is
            not periodic.
            - single Bravais layer or composite layer
            **********************************************************************/
            if( i_layer == bulk->nlayers - 1 ){
                if( (bulk->layers + i_layer)->natoms == 1){
                    ms_bravl_nd( &Tpp_s, &Tmm_s, &Rpm_s, &Rmp_s,
                            v_par, (bulk->layers + i_layer), beams_set);
                }
                else{
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

        }  /* for i_set */

        /*********************************************************************
        OVERLAYER
        Loop over all overlayer layers
        *********************************************************************/

        for(i_layer = 0; i_layer < over->nlayers; i_layer ++)
        {
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

                R_tot = ld_2lay_rpm(R_tot, R_bulk, Tpp_s, Tmm_s, Rpm_s, Rmp_s,
                                    beams_now, vec);
            }
            else
            {

                R_tot = ld_2lay_rpm(R_tot, R_tot, Tpp_s, Tmm_s, Rpm_s, Rmp_s,
                                beams_now, (over->layers + i_layer)->vec_from_last);
            }

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
        out_int(Amp, beams_now, beams_out, v_par, &results.iv_curves[energy_index * results.n_beams]);

    }  /* end of energy loop */

    return results;
}
