/*
 *    sca_grid.c
 *
 *    David Garen  1/94
 *
 *    Calculates snow covered area (in percent) by determining number of
 *    grid cells above snow line
 */

#include <stdio.h>

#include "dk_x.h"

float sca_grid(sl)
float sl;                        /* snow line elevation (thousands) */
{
   int i;                        /* loop index */
   int n = 0;                    /* counter */

   for (i = 0; i < ngrid; i++) {
      if (grid[i].elev > sl) {
         if (imask == 0 || (imask == 1 && grid[i].mask == 1))
            n++;
      }
   }
/* Debug
fprintf(fpout, "\nn = %d   nmask = %d   ngrid = %d", n, nmask, ngrid);
   End debug */
   if (imask == 0)
      return (float) ((100.0 * ((float) n) / ((float) ngrid)));
   else
      return (float) ((100.0 * ((float) n) / ((float) nmask)));
}
