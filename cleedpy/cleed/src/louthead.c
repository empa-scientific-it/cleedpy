/*********************************************************************
GH/11.08.95
  file contains function:

  int out_head( FILE * outfile)

 Write header information to output file.

 Changes:

 GH/11.08.95 - Creation

*********************************************************************/

#include <stdio.h>
#include <time.h>

#include "leed.h"

/*
#define WARNING
*/

#define CONTROL
#define ERROR

#define EXIT_ON_ERROR

int out_head(struct cryst_str *cryst_par,
             FILE * outfile)

/************************************************************************

 write header information to output file.

 INPUT:

  struct cryst_str *bulk_par - (input)
  FILE * outfile - (input) pointer to the output file were the output
            is written to.

 DESIGN:

  "#vn" No of version.
  "#ts" start time and date.

 RETURN VALUES:

   1   if o.k.
  -1   if failed

*************************************************************************/
{

struct tm *l_time;
time_t t_time;

 fprintf(outfile, "# ####################################### #\n");
 fprintf(outfile, "#            output from CLEED            #\n");
 fprintf(outfile, "# ####################################### #\n");

/************************************************************************
  Write version number and start time to output file
************************************************************************/

 fprintf(outfile, "#vn 0.0.1\n");

 t_time = time(NULL);
 l_time = localtime(&t_time);

 fprintf(outfile, "#ts %s", asctime(l_time) );
 fprintf(outfile, "#\n");

#ifdef CONTROL
 fprintf(STDCTR,"(out_head): Start date: %s", asctime(l_time) );
#endif


 return(1);
}  /* end of function out_head */
/************************************************************************/
