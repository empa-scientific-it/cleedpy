#include <stdio.h>
#include <string.h>

#include "gh_stddef.h"
#include "leed.h"

main(){
  // File names
  char *bul_file="Ni111_Cu.inp"; // bulk and non-geometrical parameters.
  char *par_file="Ni111_Cu.inp"; // overlayer parameters of all parameters (if bul_file does not exist)
  char *res_file="leed.res"; // output file containing the LEED intensities

  // Calculate LEED intensities
  leed(bul_file, par_file, res_file);

  exit(0);
}
