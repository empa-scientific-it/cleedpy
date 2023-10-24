/* mgjf 18.7.2014
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mat.h"

void mat2cblas ( real *cblas_mx, int cblas_num, mat Mx ) {

   int i, incre;
   real *cblas_px;
   real *ptrx, *ptix;

   /* we permit real<-real, complex<-real, complex<-complex */
   if ( (cblas_num == NUM_REAL) && (Mx->num_type == NUM_REAL) ) {
     incre = 1;
   } else if ( cblas_num == NUM_COMPLEX ) {
     incre = 2;
   } else {
     fprintf(stderr, "mat2cblas: invalid operand type: %d %d\n",
       cblas_num, Mx->num_type);
     exit(1);
   }

   if (Mx->num_type ==  NUM_REAL) {
     for ( i = 1 , cblas_px = cblas_mx , ptrx = Mx->rel+1 ;
           i <= Mx->rows * Mx->cols ;
           i++ , cblas_px += incre , ptrx++ ) {
       *cblas_px = *ptrx;
       if ( cblas_num == NUM_COMPLEX ) {
         *(cblas_px+1) = 0;
       }
     }

   } else if (Mx->num_type ==  NUM_COMPLEX) {
     for ( i = 1 , cblas_px = cblas_mx , ptrx = Mx->rel+1 , ptix = Mx->iel+1 ;
           i <= Mx->rows * Mx->cols ;
           i++ , cblas_px += incre , ptrx++ , ptix++ ) {
       *cblas_px = *ptrx;
       *(cblas_px+1) = *ptix;
     }

   } else {
     fprintf(stderr, "mat2cblas: invalid operand type: %d\n", Mx->num_type);
     exit(1);
   }
}

void cblas2mat ( mat Mx, real *cblas_mx ) { /* silently assume cblas same type as Mx */
   int i;
   real *cblas_px;
   real *ptrx, *ptix;

   if (Mx->num_type ==  NUM_REAL) {
     for ( i = 1 , cblas_px = cblas_mx , ptrx = Mx->rel+1 ;
           i <= Mx->rows * Mx->cols ;
           i++ , cblas_px++ , ptrx++ )
     {
       *ptrx = *cblas_px;
     }

   } else if (Mx->num_type ==  NUM_COMPLEX) {
     for ( i = 1 , cblas_px = cblas_mx , ptrx = Mx->rel+1 , ptix = Mx->iel+1 ;
           i <= Mx->rows * Mx->cols ;
           i++ , cblas_px += 2 , ptrx++ , ptix++ )
     {
       *ptrx = *cblas_px;
       *ptix = *(cblas_px+1);
     }

   } else {
     fprintf(stderr, "mat2cblas: invalid operand type: %d\n", Mx->num_type);
     exit(1);
   }
}

void info_check(const char *routine, const int info) {
  if ( info != 0 ) {
    fprintf(stderr, "%s failed: info = %d\n", routine, info);
    exit(1);
  }
}
