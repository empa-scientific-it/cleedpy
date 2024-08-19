/*********************************************************************
  GH/03.06.94
mat matmul( mat Mr, mat M1, mat M2 )
  Matrix multiplication.

Changes:
  GH/16.08.94 - Mr can be equal to  M1 or M2
                (Diagonal matrices are not included)
  GH/26.08.94 - Error in the multiplication for complex matrices
                corrected.

  mgjf 18.07.2014 - workaround
                replace naive matrix multiplication by cblas_Xgemm

*********************************************************************/
#include <math.h>
#include <stdio.h>

/*
#ifdef USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif
*/

#include "cblas.h"
#include "cblas_f77.h"
/*
#include "f77blas.h"
*/

#include "mat_blas.h"
#include "mat.h"

/*
#define CONTROL
*/

#define ERROR

#define EXIT_ON_ERROR

mat matmul( mat Mr, mat M1, mat M2 )

/*********************************************************************
  Multiply two matrices: Mr = M1*M2

  parameters:
  Mr - pointer to the matrix containing the result of the multiplication.
       If NULL, the pointer will be created and returned.

  M1, M2 - pointers to the matrices to multiply.

  return value: Mr

*********************************************************************/

{

long int i;
long int i_c2, i_r1;
long int i_cr1, i_cr2;

int result_num_type;

real *cblas_m1, *cblas_m2, *cblas_mr; /* passed to cblas_Xgemm */

/*********************************************************************
  check input matrices
*********************************************************************/

/* check validity of the input matrices */
 if ((matcheck(M1) < 1) || (matcheck(M2) < 1))
 {
#ifdef ERROR
  fprintf(STDERR,"*** error (matmul): ivalid input matrices\n");
#endif
#ifdef EXIT_ON_ERROR
  exit(1);
#else
  return(NULL);
#endif
 }

/* check dimensions of input matrices */
 if (M1->cols != M2->rows)
 {
#ifdef ERROR
  fprintf(STDERR,
  "*** error (matmul): dimensions of input matrices do not match\n");
#endif
#ifdef EXIT_ON_ERROR
  exit(1);
#else
  return(NULL);
#endif
 }

/*********************************************************************
  Create cblas matrices
*********************************************************************/

// printf("matmul: (%d,%d) x (%d,%d)\n", M1->rows, M1->cols, M2->rows, M2->cols);

 if((M1->num_type ==  NUM_REAL) && (M2->num_type ==  NUM_REAL) )
 {
   /* in this case we
   ** - need no intermediary storage for operands
   ** - need no conversion to complex
   */
   cblas_mr = calloc(M1->rows * M2->cols, sizeof(real));

   cblas_m1 = M1->rel + 1;   /* matrices are stored row major */
   cblas_m2 = M2->rel + 1;

   result_num_type = NUM_REAL;

 }
 else
 {
   /* at least one operand is complex */
   cblas_mr = calloc(M1->rows * M2->cols, 2*sizeof(real));

   cblas_m1 = calloc(M1->rows * M1->cols, 2*sizeof(real));
   cblas_m2 = calloc(M2->rows * M2->cols, 2*sizeof(real));

   mat2cblas ( cblas_m1, NUM_COMPLEX, M1 ) ;
   mat2cblas ( cblas_m2, NUM_COMPLEX, M2 ) ;

   result_num_type = NUM_COMPLEX;
 }

/*********************************************************************
  Perform the multiplication
*********************************************************************/
#ifdef CONTROL
  fprintf(STDCTR," (matmul) start multiplication \n");
#endif


  switch(result_num_type)
  {
   case (NUM_REAL):
   {
     if ( sizeof(real) == sizeof(float) ) {
       float alpha = 1.0 ;
       float beta = 0.0 ;
       cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                   M1->rows, M2->cols, M1->cols,
		   alpha, (float*)cblas_m1, M1->cols, /* lda=cols row major */
		   (float*)cblas_m2, M2->cols,
		   beta, (float*)cblas_mr, M2->cols);
     }
     else if ( sizeof(real) == sizeof(double) ) {
       double alpha = 1.0 ;
       double beta = 0.0 ;
       cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                   M1->rows, M2->cols, M1->cols,
		   alpha, (double*)cblas_m1, M1->cols, /* lda=cols row major */
		   (double*)cblas_m2, M2->cols,
		   beta, (double*)cblas_mr, M2->cols);
     } else {
       fprintf(stderr, "matmul: unexpected sizeof(real)=%lu\n", sizeof(real));
       exit(1);
     }
   }  /* endcase NUM_REAL */
   break;
   case (NUM_COMPLEX):
   {
     if ( sizeof(real) == sizeof(float) ) {
       float alpha[2] = { 1.0, 0.0 } ;
       float beta[2] =  { 0.0, 0.0 } ;
       cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                   M1->rows, M2->cols, M1->cols,
		   alpha, (float*)cblas_m1, M1->cols, /* lda=cols row major */
		   (float*)cblas_m2, M2->cols,
		   beta, (float*)cblas_mr, M2->cols);
     }
     else if ( sizeof(real) == sizeof(double) ) {
       double alpha[2] = { 1.0, 0.0 } ;
       double beta[2] =  { 0.0, 0.0 } ;
       cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                   M1->rows, M2->cols, M1->cols,
		   alpha, (double*)cblas_m1, M1->cols, /* lda=cols row major */
		   (double*)cblas_m2, M2->cols,
		   beta, (double*)cblas_mr, M2->cols);
     } else {
       fprintf(stderr, "matmul: unexpected sizeof(real)=%lu\n", sizeof(real));
       exit(1);
     }
   }  /* endcase NUM_COMPLEX */
   break ;
  }   /* endswitch */

/*
  Copy cblas_mr into Mr
*/

  Mr = matalloc(Mr, M1->rows, M2->cols, result_num_type);
  cblas2mat ( Mr, cblas_mr );

/* WAS
  Mr = matcop(Mr, Maux);
  matfree(Maux);
*/
  if ( result_num_type == NUM_COMPLEX ) {
    free(cblas_m1);
    free(cblas_m2);
  }
  free(cblas_mr);
  return(Mr);

}  /* end of function matmul */
