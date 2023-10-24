/*********************************************************************
GH/29.09.00

functions for the leed programs (sym and temp)

Changes:
GH/30.08.97

version SYM 1.1 / TEMP 0.5
GH/27.09.00 - version SYM 1.1 / TEMP 0.5
GH/29.09.00 - remove unused functions.
*********************************************************************/

/*********************************************************************
 Input
*********************************************************************/
   /* Decide which atoms belong to which layer; file linplayer.c */
int inp_bul_layer(struct cryst_str *, struct atom_str *, real *);
int inpbullay_sym(struct cryst_str *, struct atom_str *, real *);
real inp_debtemp(real , real , real );
int inp_ovl_layer(struct cryst_str *, struct atom_str *);
int inpovlay_sym(struct cryst_str *, struct atom_str *);

   /* read a matrix with specified l1 m1 l2 m2 */
mat inp_mat_lm(mat , int , const char *);

   /* read phase shifts; file linpphase.c */
int inp_phase   ( char * , real * , struct phs_str **);
int inp_phase_nd( char * , real * , int , struct phs_str **);
int upd_phase( int );

   /* read bulk parameters; file linprdbul.c */
int inp_rdbul   (struct cryst_str ** , struct phs_str ** , char *);
int inp_rdbul_nd(struct cryst_str ** , struct phs_str ** , char *);
int inp_rdbulsym(struct cryst_str ** , struct phs_str ** , char *);

   /* read overlayer parameters; file linprdovl.c */
int inp_rdovl   (struct cryst_str ** , struct phs_str ** , struct cryst_str * , char *);
int inp_rdovl_nd(struct cryst_str **, struct phs_str **, struct cryst_str *, char *);
int inp_rdovlsym(struct cryst_str ** , struct phs_str ** , struct cryst_str * , char *);
   /* read other parameters; file linprdpar.c */
int inp_rdpar(struct var_str **, struct eng_str **, struct cryst_str * , char *);
   /* show all parameters; file linpshowbop.c */
int inp_showbop(struct cryst_str *, struct cryst_str *, struct phs_str *);
   /* read and write parameters */
int write_par(struct cryst_str *, struct phs_str *, struct var_str *,
              struct eng_str *, struct beam_str *, FILE * );
int read_par(struct cryst_str **, struct phs_str **, struct var_str **,
             struct eng_str **, struct beam_str **, FILE * );

int ceckrotsym(struct cryst_str *);
int ceckmirsym(struct cryst_str *);

/*********************************************************************
 beams (bm) and parameter control (pc) and output (out)
*********************************************************************/
    /* Free rotation matrices (lbmrotmat.c) */
int bm_freerotmat(real ** );
    /* Find the beams to be included (lbmgen.c) */
int bm_gen(struct beam_str **, struct cryst_str *,
             struct var_str *, real);
    /* Find the beams to be included (lbmgenrot.c) */
int bm_gen_sym(struct beam_str **, struct cryst_str *, struct cryst_str *,
             struct var_str *, real);

    /* create rotation matrices */
real ** bm_rotmat(int );
    /* Find the beams to be included at the current energy (lbmselect.c) */
int bm_select(struct beam_str **, struct beam_str *,
             struct var_str *, real);
    /* Find the beams of a particular beam set (lbmset.c) */
int bm_set(struct beam_str **, struct beam_str *, int);

/*********************************************************************
 Parameter control
*********************************************************************/
    /* reset parameter structure (lpcreset.c) */
int pc_reset(struct var_str *, struct cryst_str *);

    /* update energy (lpcupdate.c) and tl (lpcmktl.c) */
int pc_update(struct var_str *, struct phs_str *, real);
mat *pc_mktl(mat *, struct phs_str *, int, real);
mat *pc_mktl_nd(mat *, struct phs_str *, int, real);

    /* temperature dependent scattering factors */
mat pc_temtl(mat , mat , real , real , int , int );
mat pc_cumtl(mat , mat , real , real , real , real , int , int );
int pc_mk_ms(mat * , mat *, mat *, mat *, mat *, mat *, int );

/*********************************************************************
 Output
*********************************************************************/
int out_head(struct cryst_str *, FILE *);
int out_head_2(struct cryst_str *, const char *, const char *, FILE *);
int out_bmlist(struct beam_str **, struct beam_str *, struct eng_str *, FILE *);
int out_int(mat , struct beam_str *, struct beam_str *, struct var_str *, FILE * );
int out_intsym(mat , struct beam_str *, struct beam_str *, struct var_str *, FILE * );

    /* check cpu time */
double cpu_time(FILE *, const char *);

/*********************************************************************
 Layer doubling
*********************************************************************/
   /* LD for 2 layers */
int ld_2lay (mat *, mat *, mat *, mat *, 
             mat, mat, mat, mat, mat, mat, mat, mat, 
             struct beam_str *, real *);
mat ld_2lay_rpm (mat, mat, mat, mat, mat, mat,
             struct beam_str *, real *);
   /* LD for periodic layers */
mat ld_2n (mat, mat, mat, mat, mat, struct beam_str *, real *);
   /* LD for potential step */
mat ld_potstep ( mat , mat , struct beam_str *, real , real *);
mat ld_potstep0 ( mat , mat , struct beam_str *, real , real *);

/*********************************************************************
 Multiple scattering
*********************************************************************/

   /* Don't know yet */
int ms_bravl ( mat *, mat *,
               struct var_str *, struct layer_str *, struct beam_str *);
int ms_bravl_nd ( mat *, mat *, mat *, mat *,
               struct var_str *, struct layer_str *, struct beam_str *);
int ms_bravl_sym ( mat *, mat *,
               struct var_str *, struct layer_str *, struct beam_str *);
int ms_compl ( mat *, mat *, mat *, mat *, 
               struct var_str *, struct layer_str *, struct beam_str *);
int ms_compl_nd ( mat *, mat *, mat *, mat *, 
               struct var_str *, struct layer_str *, struct beam_str *);
int ms_complsym ( mat *, mat *, mat *, mat *,  
               struct var_str *, struct layer_str * ,struct beam_str *);

   /* lattice sum for one layer (lmslsumii.c) */
mat ms_lsum_ii (mat , real , real , real * , real * , int , real );

   /* lattice sum for two layers (lmslsumij(sym).c) */
int ms_lsum_ij (mat *, mat *, real , real , real * , real * , real *, int , real );
mat ms_lsum_ij_sym (mat, real , real , real * , real * , real *, int , real, int );

    /* partial inversion */
mat ms_partinv ( mat , mat , int , int );

   /* Green's function (lmstmatii/ij/ijsym.c, lmsgmatijsym.c) */
mat ms_tmat_ii (mat , mat, mat, int );
mat ms_tmat_nd_ii (mat , mat, mat, int );
mat ms_tmat_ij (mat , mat, mat, int );
mat ms_tmat_ij_sym (mat, mat, mat, int, int );

   /* Transformation L -> k (lmsymat.c/lmsymmat.c) */
mat ms_ymat  (mat , int , struct beam_str *, int );
mat ms_ymat_set  (mat , int , struct beam_str *, int );
mat ms_ymmat (mat , int , struct beam_str *, int );
mat ms_ymat_r (mat , int , struct beam_str *, int );

   /* Transformations of Ylm (lmsypy.c) */
mat ms_yp_ym  (mat , mat);
mat ms_yp_yxp (mat , mat);
mat ms_yp_yxm (mat , mat);

   /* Summation over all beams for composite layer in case of symmetry */
mat mscompksum(mat ,struct beam_str * ,struct atom_str * ,int, int, int );

/*********************************************************************
lower level functions
*********************************************************************/

