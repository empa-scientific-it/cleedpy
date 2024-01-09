/*********************************************************************
  GH/08.06.94

 int matfree( mat M )
     Free the memory allocated for a matrix

 Changes:

 GH/26.08.94 - Remove MAT_ERROR
*********************************************************************/
#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif
#include "mat.h"

/*
#define CONTROL
*/
#define ERROR

int matfree( mat M )

/*********************************************************************
  Free the memory allocated for a matrix

  parameters:
  M - pointer to the matrix.

  return value: 1 if successful, 0 if failed.

*********************************************************************/
{
/*
  check input matrix
*/
 if ( matcheck(M) < 1)
 {
#ifdef ERROR
   fprintf(STDERR," *** error (matfree): improper input \n");
#endif
   return(0);
 }

 if (M->rel != NULL) free(M->rel);
 if (M->iel != NULL) free(M->iel);

 free(M);
 return(1);

}  /* end of function matfree */
