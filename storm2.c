/*
 *    storm2.c
 *
 *    David Garen  9/92
 *
 *    Compute MAP based on storms.
 */

#include <stdio.h>

#include "dk_x.h"

void storm2()
{
   int i, j, k, l, m, n;         /* loop indexes */
   int ns;                       /* number of stations with data */
   int *staflg;				 	 /* station use flags*/

   // set station use flags
   staflg = ivector(nsta);
   for (m = 0; m < nsta; m++)
	   staflg[m] = 1;


   /* Storm loop */


   for (m = 0; m < nstorm; m++) {
      if (b0[m][0] <= 99998 && b1[m][0] <= 99998) {

         /* Day loop */

         j = storm[m].dstart;
         k = storm[m].ystart;
         dstop = lastday[k] - 1;
         for (n = 0; n < storm[m].slen; n++) {
            if (map[j][k] > missing) {

               /* Check to see if any stations have missing data */

               imiss = ns = 0;
               for (i = 0; i < nsta; i++)
                  if (sta[i].data[j][k] < accum)
                     ns++;
               if (ns <= 1)
                  continue;
               if (ns < nsta) {
                  imiss = 1;
                  if (iwt == 2) {
                     for (i = 0; i < nsta; i++) {
                        if (sta[i].data[j][k] < accum)
                           w[i] = 1.0 / ns;
                        else
                           w[i] = 0.0;
                     }
                  }
               }

               /* Grid loop */

               dum = 0;
               for (l = 0; l < ngrid; l++) {
                  if (grid[l].use == 1) {

                     /* If one or more stations have missing data,
                        calculate kriging weights excluding those stations;
                        otherwise, use weights for all stations
                        that have already been calculated */

                     if (iwt == 1 && imiss == 1) {
                        for (i = 0; i < nsta; i++) {
                           if (sta[i].data[j][k] < accum)
                              staflg[i] = 1;
                           else
                              staflg[i] = 0;
                        }
                        krige(l, ns);
                     }

                     /* Compute detrended precipitation at grid cell */

                     gprec[l] = 0;
                     if (imiss == 1) {
                        for (i = 0; i < nsta; i++)
                           gprec[l] += (float) ((w[i] * sta[i].data[j][k]));
                     }
                     else {
                        for (i = 0; i < nsta; i++)
                           gprec[l] += (wall[l][i] * sta[i].data[j][k]);
                     }

                     /* Re-trend grid prec/temp */

                     gprec[l] += (b0[m][0] + b1[m][0] * grid[l].elev);

                     /* Set grid prec values to zero if estimate is
                        less than zero */

                     if (gprec[l] < 0)
                        gprec[l] = 0;

                     /* Add grid prec to basin sum */

                     if (imask == 0 || (imask == 1 && grid[l].mask == 1))
                        dum += gprec[l];
                  }
               }

               /* If requested, write out grid in GRASS format */

               if (iout == 2 && j >= igridout1 && j <= igridout2)
                  grassout(year[k], j);

               /* If requested, write out grid in ARC/INFO format */

               if (iout == 3 && j >= igridout1 && j <= igridout2)
                  arcout(year[k], j);

               /* If requested, write out grid in IPW format */

               if (iout == 4 && j >= igridout1 && j <= igridout2)
                  ipwout(year[k], j);

               /* If requested, compute and write out zonal means for day */

               if (izone == 1)
                  zoneout(year[k], j);

               /* Compute MAP for day */

               if (imask == 0)
                  map[j][k] = dum / ngriduse;
               else
                  map[j][k] = dum / nmask;
            }
            j++;
            if (j > dstop) {
               j = 0;
               k++;
            }
         }
      }
   }
}
