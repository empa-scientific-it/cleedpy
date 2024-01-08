#include "leed.h"

double * leed(
    struct cryst_str * bulk,
    struct cryst_str *over,
    struct phs_str *phs_shifts,
    struct beam_str *beams_all,
    int n_set,
    struct eng_str *eng){

    struct beams_str *beams_now, *beams_set;
    int n_beams_now, n_beams_set;
    float energy;
    struct var_str *v_par=NULL;
    float **iv_curves=NULL;

    mat R_bulk;

    /* Main Energy Loop */
    for(energy=eng->ini; energy<eng->fin+E_TOLERANCE; energy+= eng->stp){
        pc_update(v_par, phs_shifts, energy);
        n_beams_now = bm_select(&beams_now, beams_all, v_par, bulk->dmin);

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
    }

    return iv_curves;
}
