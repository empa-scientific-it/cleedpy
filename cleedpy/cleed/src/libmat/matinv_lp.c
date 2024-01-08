/*********************************************************************
  GH/08.06.94
  file contains function:

  matinv

Changes
GH/08.06.94 - Creation
GH/20.07.95 - Change call of function c_luinv

  mgjf 18.07.2014 - workaround
                replace NumRec inversion by LAPACK routines

*********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mat_blas.h"
#include "mat_lapack.h"
#include "mat.h"


/*
#define CONTROL
*/
#define ERROR

/********************************************************************/

mat matinv( mat A_1, mat A)

/*********************************************************************
  Invert a matrix by LU decomposition (real or complex).

  parameters:
  A_1- input: pointer to the inverse matrix.
  A  - input: pointer to the matrix to be inverted.

  return value:
     pointer to the inverse matrix.

*********************************************************************/

{
int i, d, n, nb, lwork;
const int p1 = 1, m1 = -1;
int *ipiv;
int info;

real *ptr, *ptr_end;

real *cblas_a, *work;

/*********************************************************************
  check input matrix
*********************************************************************/
/* check type of input matrix */
 if ( (A->mat_type != MAT_SQUARE) && (A->mat_type != MAT_DIAG) )
 {
#ifdef ERROR
  fprintf(STDERR," *** error (matinv): improper input type of matrix\n");
#endif
  return(NULL);
 }

/* check dimensions of input matrix A */
 else if (A->cols != A->rows)
 {
#ifdef ERROR
  fprintf(STDERR,
          " *** error (matinv): dimensions of input matrix do not match\n");
#endif
  return(NULL);
 }
 n = A->cols;

/*********************************************************************
  DIAGONAL MATRIX
*********************************************************************/

/*
 If the matrix is diagonal, the inversion is trivial
 (but currently not implemented)
*/

 if (A->mat_type == MAT_DIAG)
 {
#ifdef ERROR
   fprintf(STDERR,
   " *** error (matinv): inversion of diagonal matrices is not implemented \n");
#endif
   return(NULL);
 }

 else

/*********************************************************************
  NON-DIAGONAL MATRIX
*********************************************************************/

 {

/*********************************************************************
  Create lapack and result matrices
*********************************************************************/

 ipiv = (int *)calloc( (n+1), sizeof(int));

 switch(A->num_type)
 {
                    /*******************************
                         Real matrix (NUM_REAL)
                    *******************************/
   case (NUM_REAL):
   {

      cblas_a = calloc(n*n, sizeof(real)) ;
      mat2cblas( cblas_a, NUM_REAL, A);

      if ( sizeof(real) == sizeof(float) ) {
        nb = ilaenv_( &p1, "SGETRI", " ", &n, &m1, &m1, &m1);
        lwork = n*nb;
        work = calloc(lwork, sizeof(float));
        sgetrf_(&n, &n, cblas_a, &n, ipiv, &info);
	info_check("sgetrf", info);
        sgetri_(&n, (float*)cblas_a, &n, ipiv, work, &lwork, &info);
	info_check("sgetri", info);
	free(work);
      }
      else if ( sizeof(real) == sizeof(double) ) {
        nb = ilaenv_( &p1, "DGETRI", " ", &n, &m1, &m1, &m1);
        lwork = n*nb;
        work = calloc(lwork, sizeof(double));
        dgetrf_(&n, &n, cblas_a, &n, ipiv, &info);
	info_check("dgetrf", info);
        dgetri_(&n, (double*)cblas_a, &n, ipiv, work, &lwork, &info);
	info_check("dgetri", info);
	free(work);
     } else {
       fprintf(stderr, "matinv: unexpected sizeof(real)=%lu\n", sizeof(real));
       exit(1);
     }


#ifdef CONTROL
     fprintf(STDCTR, " (matinv) ipiv: ");
     for (i=0; i<n; i++) fprintf(STDCTR," %d", ipiv[i]);
     fprintf(STDCTR,"\n");
#endif

     /* Allocate A_1 (if it does not exist.) */
     A_1 = matalloc( A_1, n, n, A->num_type);
     cblas2mat(A_1, cblas_a);
     free(cblas_a);
     break;
   }  /* REAL */

                    /*******************************
                       Complex matrix (NUM_COMPLEX)
                    *******************************/

   case (NUM_COMPLEX):
   {

      cblas_a = calloc(n*n, 2*sizeof(real)) ;
      mat2cblas( cblas_a, NUM_COMPLEX, A);

      if ( sizeof(real) == sizeof(float) ) {
        nb = ilaenv_( &p1, "CGETRI", " ", &n, &m1, &m1, &m1);
        lwork = n*nb;
        work = calloc(lwork, 2*sizeof(float));
        cgetrf_(&n, &n, cblas_a, &n, ipiv, &info);
	info_check("cgetrf", info);
        cgetri_(&n, (float*)cblas_a, &n, ipiv, work, &lwork, &info);
	info_check("cgetri", info);
	free(work);
      }
      else if ( sizeof(real) == sizeof(double) ) {
        nb = ilaenv_( &p1, "ZGETRI", " ", &n, &m1, &m1, &m1);
        lwork = n*nb;
        work = calloc(lwork, 2*sizeof(double));
        zgetrf_(&n, &n, cblas_a, &n, ipiv, &info);
	info_check("zgetrf", info);
        zgetri_(&n, (double*)cblas_a, &n, ipiv, work, &lwork, &info);
	info_check("zgetri", info);
	free(work);
     } else {
       fprintf(stderr, "matinv: unexpected sizeof(real)=%lu\n", sizeof(real));
       exit(1);
     }

#ifdef CONTROL
     fprintf(STDCTR, " (matinv) ipiv: ");
     for (i=0; i<n; i++) fprintf(STDCTR," %d", ipiv[i]);
     fprintf(STDCTR,"\n");
#endif

     /* Allocate A_1 (if it does not exist.) */
     A_1 = matalloc( A_1, n, n, A->num_type);
     cblas2mat(A_1, cblas_a);
     free(cblas_a);
     break;
   }     /* COMPLEX */

 }   /* switch num_type */
 free(ipiv);
 return(A_1);
 } /* non-diagonal */
}
/********************************************************************/
