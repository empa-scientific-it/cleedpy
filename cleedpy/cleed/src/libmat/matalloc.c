/*********************************************************************
  GH/02.03.95

  matalloc
     Allocate memory for a cols x rows matrix.
     all matrix elements are set to zero.

  Changes:

  GH/11.08.94 - Diagonal matrices
  GH/15.08.94 - set all matrix elements to zero.
  GH/26.08.94 - num_type has a different meaning: num_type + mat_type.
  GH/20.01.95 - default blk_type = BLK_SINGLE

*********************************************************************/

#include <stdio.h>
#if defined (__MACH__)
  #include <stdlib.h>
#else
  #include <malloc.h>
#endif
#include "mat.h"

/*
#define ERROR
#define CONTROL
#define CONTROL_X
*/
#define ERROR
#define EXIT_ON_ERROR

/*======================================================================*/
/*======================================================================*/

mat matalloc(mat M, int rows, int cols, int num_type)
/*********************************************************************

  Allocate memory for a cols x rows matrix of type num_type. All matrix
  elements are set to zero.

  INPUT
    mat M
    int rows     - max first index (i.e. number of rows)
    int cols     - max second index  (i.e. number of columns)
    int num_type - type of matrix elements (NUM_REAL or NUM_COMPLEX)
                   and type of matrix (MAT_NORMAL, MAT_DIAG, etc.).
                   input as:
                   mat_type | num_type.

  RETURN VALUE:
    pointer to the matrix (mat) (if successful)
    NULL else.

*********************************************************************/
{
    int mat_type;
    size_t no_of_elts;
    real *ptr, *ptr_end;

    /*********************************************************************
     Check the validity of the pointer M and cols/rows first
    *********************************************************************/

    if(matcheck(M) < 0)
    {
        fprintf(STDERR," *** error (matalloc): Invalid pointer \n");
        exit(1);
    }

    if( (cols < 1) || (rows < 1) )
    {
        fprintf(STDERR," *** error (matalloc): Invalid number of cols/rows \n");
        exit(1);
    }

    /*********************************************************************
     Find the correct matrix type (high byte of num_type) and number
    types (low byte of num_type).
    *********************************************************************/

    mat_type = num_type & MAT_MASK;

    if(mat_type == 0)
    {
        if(cols != rows) mat_type = MAT_NORMAL;
        else
            if(cols != 1) mat_type = MAT_SQUARE;
            else mat_type = MAT_SCALAR;
    }
    else if( (mat_type == MAT_DIAG) && (cols != rows) )
    {
        fprintf(STDERR," *** error (matalloc): different numbers of cols and rows are incompatible with diagonal matrix type\n");
        exit(1);
    }

    num_type = num_type & NUM_MASK;
    /*******************************************************************************
     If M points to the right matrix type already, only reset all matrix elements.
    *******************************************************************************/

    if( (matcheck(M) > 0) && (M->cols == cols) && (M->rows == rows) && (M->num_type == num_type) && (M->mat_type == mat_type) )
    {
        if(M->num_type == NUM_COMPLEX)
        {
            // M->mat_type < MAT_DIAG means square, normal or scalar
            if(M->mat_type < MAT_DIAG)
            {
                ptr_end = M->rel + M->rows*M->cols;
                for(ptr = M->rel; ptr <= ptr_end; ptr ++) *ptr = 0.;

                ptr_end = M->iel + M->rows*M->cols;
                for(ptr = M->iel; ptr <= ptr_end; ptr ++) *ptr = 0.;
            }
            if(M->mat_type == MAT_DIAG)
            {
                ptr_end = M->rel + M->cols;
                for(ptr = M->rel; ptr <= ptr_end; ptr ++) *ptr = 0.;

                ptr_end = M->iel + M->cols;
                for(ptr = M->iel; ptr <= ptr_end; ptr ++) *ptr = 0.;
            }
        } /* NUM_COMPLEX */

        if(M->num_type == NUM_REAL)
        {
            if(M->mat_type < MAT_DIAG)
            {
            ptr_end = M->rel + M->rows*M->cols;
            for(ptr = M->rel; ptr <= ptr_end; ptr ++) *ptr = 0.;
            }
            if(M->mat_type == MAT_DIAG)
            {
            ptr_end = M->rel + M->cols;
            for(ptr = M->rel; ptr <= ptr_end; ptr ++) *ptr = 0.;
            }
        } /* NUM_REAL */
        // Set magic number and return
        M->mag_no = MATRIX;
        return(M);
    }  /* reset existing matrix */

    /*********************************************************************
    If not: allocate new memory

    - if M == NULL: allocate a new matrix header, blk_type = BLK_SINGLE
    - if M != NULL: update the old header, use former blk_type.
    *********************************************************************/

    if( M == NULL)
    {
        M = ( mat )malloc(sizeof(struct mat_str));
        M->blk_type = BLK_SINGLE;
    }
    else  // M != NULL
    {
        if (M->iel != NULL) free(M->iel);
        if (M->rel != NULL) free(M->rel);
    }

    M->cols = cols;
    M->rows = rows;
    M->num_type = num_type;
    M->mat_type = mat_type;

    if(M->mat_type == MAT_DIAG)
        no_of_elts = M->cols+1;
    else
        no_of_elts = M->rows * M->cols + 1;

    switch(num_type)
    {
        case(NUM_REAL):
        {
            M->iel = NULL;
            M->rel = (real*)calloc( no_of_elts, sizeof(real));
            if (M->rel == NULL)
            {
                free(M);
                fprintf(STDERR,"*** error (matalloc) allocation error\n");
                exit(1);
            }
            break;
        }  /* NUM_REAL */

        case(NUM_COMPLEX):
        {
            M->rel = (real*)calloc( no_of_elts, sizeof(real));
            M->iel = (real*)calloc( no_of_elts, sizeof(real));

            if( (M->rel == NULL) || (M->iel == NULL) )
            {
                if (M->rel != NULL) free(M->rel);
                if (M->iel != NULL) free(M->iel);
                free(M);
                fprintf(STDERR,"*** error (matalloc) allocation error\n");
                exit(1);
            }
            break;
        }  /* NUM_COMPLEX */

        default: // Unknown number type.
        {
            fprintf(STDERR, " *** error (matalloc) wrong number type: %d\n", num_type);
            exit(1);

        }
    } /* switch */

    // Set magic number and return
    M->mag_no = MATRIX;
    return(M);

}  /* end of function matalloc */
