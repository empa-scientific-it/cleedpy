/*********************************************************************
GH/06.01.2012

include file for
 - additional data structures and type definitions
 - constant values
in the LEED programs

Changes:
GH/20.09.95
ST/16.07.97
GH/02.09.97 - change beam_str and cryst_str (add. elements),
              define S0, SX, SY, SXY.
WB/07.08.98 - add k_rsym (symmetry related beams) to beam_str.
WB/27.08.98 - include symmetry flags in crystal struct
GH/03.05.00 - include t_type in atom_str.
            - include t_type in phs_str.

version SYM 1.1 + TEMP 0.5
GH/27.09.00 - same include file for version SYM 1.1 + TEMP 0.5
GH/06.01.12 - set MIN_DIST to 1.9 (1.0 A)

*********************************************************************/

#ifndef LEED_DEF_H
#define LEED_DEF_H

#include "real.h"
#include "mat_def.h"

/*********************************************************************
 structures and types
*********************************************************************/

/*********************************************************************
  struct atom_str contains all properties of a single atom
*********************************************************************/

struct atom_str
{
 int  layer;      /* No of layer where the atom belongs to */
 int  type;       /* type of atom (i. e. set of phase shifts to be used) */
 int  t_type;     /* type of t matrix (T_DIAG or T_NOND) */
 real pos[4];     /* relative position inside the unit cell */
 real dwf;        /* Debye-Waller factor */
};

/*********************************************************************
  struct layer_str contains all properties of a single layer
*********************************************************************/

struct layer_str
{
 int  no_of_layer;       /* No of layer in array */
 int  bulk_over;         /* BULK (0) or OVER */
 int  periodic;          /* 1: layer is part of the periodically repeated bulk
                               unit cell
                            0: layer is only used once
                         */
 int  natoms;            /* number of atoms in the layer */
 real a_lat[5];          /* basis vectors of the real 2-dim unit cell stored as
                            standard matrix (a1,a2): a1x = a_lat[1], a2x = a_lat[2]
                                                     a1y = a_lat[3], a2y = a_lat[4]
                         */
 real rel_area;          /* area of the unit cell relative to 1x1 */
 real reg_shift[4];      /* 2-dim vector pointing to the axis of rot. symmetry/mirror plane */
 real vec_from_last[4];  /* vector from the origin of layer (n-1) */
 real vec_to_next[4];    /* vector to the origin of layer (n+1) */
 struct atom_str *atoms; /* properties of the atoms within the composite
                            (or Bravais) layers */
};

/*********************************************************************
  struct cryst_str contains all crystal specific program parameters
*********************************************************************/
struct cryst_str
{

/* general parameters */
 real vr;         /* real part of the optical potential */
 real vi;         /* imaginary part optical potential */
 real temp;       /* crystal temperature */

/* symmetries */
 int n_rot;       /* degree of rotational symmetry */
 real rot_axis[4];/* axis of rotational symmetry */
 int n_mir;       /* number of mirror plane */
 real *m_plane;   /* points define mirror plane */
 real *alpha;     /* angle in degree */
 int symmetry;    /* NOSYM(0) ROTSYM(1) MIRRORSYM(2) RMSYM(3)*/

/* 1x1 unit cell */
 real a[5];       /* basis vectors of the real 2-dim unit cell stored as
                     standard matrix (a1,a2): a1x = a[1], a2x = a[2]
                                             a1y = a[3], a2y = a[4] */
 real a_1[5];     /* inverse of a (= 2PI * (a^-1))
                     => the rows are the basis vectors of the reciprocal
                     2-dim unit cell: a*1x = a_1[1], a*1y = a_1[2]
                                      a*2x = a_1[3], a*2y = a_1[4] */

 real area;       /* Area of the real 2-dim (1x1) unit cell */

 real m_trans[5]; /* 2x2 transformation matrix: unit cell used in the programs
                     with respect to the input: a_inp = m_trans * a_prg */

/* superstructure */
 real m_super[5]; /* 2x2 superstructure matrix: b = m_super * a */
 real m_recip[5]; /* reciprocal 2x2 superstructure matrix: Mt^-1 */

 real b[5];       /* basis vectors of the real 2-dim superstructure unit cell stored as
                    standard matrix (b1,b2): b1x = b[1], a2x = b[2]
                                             b1y = b[3], a2y = b[4] */
 real b_1[5];     /* inverse of b (= 2PI * (b^-1))
                     => the rows are the basis vectors of the reciprocal 2-dim super
                     structure unit cell: b*1x = b_1[1], b*1y = b_1[2]
                                          b*2x = b_1[3], b*2y = b_1[4] */
 real rel_area_sup; /* Area of the real 2-dim superstructure unit cell
                     relative to 1x1 */

/* stacking of layers */
 int  nlayers;    /* number of layers in array layers */
 struct layer_str *layers;
                  /* information concerning the atomic layers */
 real dmin;       /* min. interlayer distance */

 int  natoms;     /* total number of atoms */
 int  ntypes;     /* total number atom types */

 char **comments; /* comments */
};

/*********************************************************************
  struct phs_str contains all parameters concerning the phase shifts
*********************************************************************/
struct phs_str
{
 int  lmax;
 int  neng;
 int  t_type;         /* type of scattering matrix: T_DIAG or T_NOND */
 real eng_max;
 real eng_min;
 real *energy;
 real *pshift;

 real dr[4];
 char *input_file;
};

/*********************************************************************
  struct beam_str contains all parameters of a specific beam in k-space.
*********************************************************************/
struct beam_str
{
 real ind_1;     /* beam indices in (1x1) basis A (real) */
 real ind_2;
 int  b_ind_1;   /* beam indices in super-structure basis B (integer) */
 int  b_ind_2;
 real k_par;     /* length of the parallel components */
 real k_r[4];    /* k_r[0] = |k_r|, k_r[i] = kx,ky,kz */
 real k_i[4];    /* imag. part */
 real cth_r;     /* cos (theta(k)), real and imag. part */
 real cth_i;
 real phi;       /* phi(k) */
 real Akz_r;     /* factor 1/(Akz), needed to calculate the scattering matrix */
 real Akz_i;

/* symmetry related phase factors */
 real k_p_sym[12]; /* angular difference (phi) between all symmetry-related beams and the representant beam (12 is max. number of beams) */
 real k_x_sym[12]; /* x component of all symmetry-related beams represented by this beam (0) (12 is max. number of beams) */
 real k_y_sym[12]; /* y component of all symmetry-related beams represented by this beam (0) (12 is max. number of beams) */
 int  n_eqb_b;   /* number of equivalent beams represented by this beam (bulk layers) */
 int  n_eqb_s;   /* number of equivalent beams represented by this beam (superstructure) */
 real *eout_b_r; /* phase factor for the outgoing beam sqrt(n_rot) S exp(+is*g')), */
 real *eout_b_i; /* needed to calculate the scattering matrix, length = n_layer*/
 real *eout_s_r; /* b: bulk layers, s: superstructure */
 real *eout_s_i;
 real *ein_b_r;  /* phase factor for the incoming beam 1/sqrt(n_rot) S exp(-i(m*phi + s*g)), */
 real *ein_b_i;  /* needed to calculate the scattering matrix, length = (2*l_max + 1) * n_layer*/
 real *ein_s_r;  /* b: bulk layers, s: superstructure */
 real *ein_s_i;

 int  set;       /* beam set, where the beam belongs to */
};

/*********************************************************************
  struct var_str contains all parameters that change during the energy
  loop and the parameters controlling them.
*********************************************************************/
struct var_str
{
 real eng_r;    /* current energy in crystal (real  part) */
 real eng_i;    /* current energy in crystal (imag. part) */
 real eng_v;    /* current vacuum energy */
 real vr;       /* real part of the optical potential */
 real vi_pre;   /* prefactor: imag. part of the optical potential */
 real vi_exp;   /* exponent: imag. part of the optical potential */

 real theta;    /* polar angle of incidence */
 real phi;      /* azimuth angle of incidence */
 real k_in[4];  /* incident beam:
                   k[0] = |k_par| k[1/2] = k_par_x/y k[3] = k_z */

 real epsilon;  /*   */
 int  l_max;    /* max. l quantum number used in the calculation */
 mat  *p_tl;    /* array of diagonal atomic scattering matrices (1st dim = lmax, 2nd dim = 1) */
};

/*********************************************************************
  struct eng_str contains the parameters that control the energy loop.
*********************************************************************/
struct eng_str  /* contains all parameters that change during the energy loop
                   and the parameters controlling them */
{
 real ini;      /* initial energy */
 real fin;      /* final energy */
 real stp;      /* energy step */
};

/*********************************************************************
 Fundamental constants/conversion factors
 (Source: CRC Handbook, 73rd Edition)
*********************************************************************/

#define HART  27.2113962       /* Hartree in eV */
#define BOHR  0.529177249      /* Bohr radius in Angstroms */
#define MEL_U 5.48579903e-4    /* electron mass in amu (atomic mass units) */
#define U_MEL 1822.88850636    /* 1 amu in multiples of the electron mass */
#define KB    3.16682941e-6    /* Boltzmann constant in Hartree/K */

/*********************************************************************
 Program parameters
  - tolerance
  - threshold values etc.
*********************************************************************/

/* Defaults */

#define PHASE_PATH "CLEED/PHASE"
                               /* default path name for phase shift input
                                  relative to HOME directory */
#define DEF_TEMP 300.          /* default temperature for calculating thermal
                                  vibrations */
#define R_FOR_LMAX 2.0         /* default muffin-tin-radius (in BOHR) used to
                                  calculate l_max (if not given) */

/* Threshold values */

#define MIN_DIST 1.9           /* Min distance for Layer Doubling in BOHR
                                  (1.0 A) */
#define VI_START 3.67493       /* start value for variable imag. opt. potential
                                  (100. eV/HART) */

/* Tolerances */

#define E_TOLERANCE    1.e-3   /* tolerance for energies in HARTREE (0.027eV)*/
#define GEO_TOLERANCE  1.e-3   /* tolerance for geometrical parameters in BOHR */
#define INT_TOLERANCE  1.e-10  /* min. intensity > 0. */
#define K_TOLERANCE    1.e-4   /* tolerance for k_par in (BOHR)^-1 */
#define LD_TOLERANCE   1.e-4   /* convergence criterion for layer doubling */
#define WAVE_TOLERANCE 1.e-4   /* tolerance for wave amplitudes */

/* Flags for mirror planes etc. */

#define BULK 0
#define OVER 1

#define T_DIAG  0
#define T_NOND  1

#define S0  0
#define SX  1
#define SY  2
#define SXY 3

#define NOSYM      101
#define MONO_2ROT  102
#define MONO_1MIR  111

#define REC_2ROT   202
#define REC_1MIR   211
#define REC_2MIR   221

#define HEX_1MIR   311
#define HEX_3ROT   303
#define HEX_3MIR   331

#define SQ_2ROT    402
#define SQ_4ROT    404
#define SQ_1MIR    411
#define SQ_2MIR    421
#define SQ_4MIR    441


/* current version No. now in file leed_ver.h */

#endif /* LEED_DEF_H */
