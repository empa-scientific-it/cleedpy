/* mgjf 18.7.2014
*/

#ifndef MAT_H
#include "mat.h"
#endif

#include <cblas.h>
/*
#include <cblas.h>
#include <mkl.h>
*/


void mat2cblas ( real *cblas_mx, int cblas_num, mat Mx ) ;
void cblas2mat ( mat Mx, real *cblas_mx ) ;

void info_check(const char *routine, const int info) ;

