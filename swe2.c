/*
 *    swe2.c
 *
 *    David Garen  1/94, 3/94
 *
 *    Compute MASWE based on dpp-day periods.
 */

#include <stdio.h>

#include "dk_x.h"

void swe2()
{

   int i, j, jj, k, l, m, n;     /* loop indexes */
   int ns;                       /* number of stations with data */

   /* Year loop */

   for (k = 0; k < nyear; k++) {
      n = lastday[k] - firstday[k] + 1;
      nper = n / dpp;
      nperm1 = nper - 1;
      dppl = n - dpp * nperm1;
      nstop = dpp;

      /* Period loop */

      for (m = 0; m < nper; m++) {
         if (m == nperm1)
            nstop = dppl;
         jj = dpp * m + firstday[k] - 1;
         if (b0[m][k] <= 99998 && b1[m][k] <= 99998) {

            /* Day loop */

            for (n = 0; n < nstop; n++) {
               j = jj + n;
               /* @@@ progress "..." still needed?
               printf("   %s%d%s%d%s%d%s\n",
                      "Processing year ", year[k], ", period ", m+1,
                      ", day ", j+1, " ...");
               */
               if (map[j][k] > missing) {

                  /* Check to see if any stations have missing data */

                  imiss = ns = 0;
                  for (i = 0; i < nsta; i++)
                     if (sta[i].data[j][k] < missing)
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

                        /* Compute detrended swe at grid cell */

                        gprec[l] = 0;
                        if (grid[l].elev > snolin[m][k]) {
                           if (imiss == 1) {
                              for (i = 0; i < nsta; i++)
                                 gprec[l] += (float) ((w[i] * sta[i].data[j][k]));
                           }
                           else {
                              for (i = 0; i < nsta; i++)
                                 gprec[l] += (wall[l][i] * sta[i].data[j][k]);
                           }

                           /* Re-trend grid swe */

                           if (iswehz[m][k] >= 0 &&
                               grid[l].elev < sta[isweln[m][k]].elev) {
                              gprec[l] += (b02[m][k] + b12[m][k] * grid[l].elev);
/* Debug
printf("\nswe2: Period %d -- hz/ln retrending ...", m+1);
   End debug */
                           }
                           else if (b1[m][k] > 0.0000001)
                              gprec[l] += (b0[m][k] + b1[m][k] * grid[l].elev);

                           /* Set grid swe values to zero if estimate is
                              less than zero */

                           if (gprec[l] < 0)
                              gprec[l] = 0;
                        }

                        /* Add grid swe to basin sum */

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

                  /* Compute MASWE for day */

                  if (imask == 0)
                     map[j][k] = dum / ngriduse;
                  else
                     map[j][k] = dum / nmask;
               }
            }
         }
      }
   }
}
