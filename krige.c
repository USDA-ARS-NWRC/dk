/*
 *    krige.c
 *
 *    David Garen   8/91, 3/94
 *
 *    Calculate kriging weights
 *
 *    26 May 2000:
 *    Change solution method to LU decomposition
 */

#include <stdio.h>
#include <stdlib.h>
#include "malloc.h"

#include "dk_x.h"

void krige(l, ns)
int l;                           /* grid index */
int ns;                          /* number of stations used */
{
   float elevsave;               /* stored value of station elevation */
   int m, mm, n, nn;             /* loop indexes */
   int msave;                    /* stored value of m index */
   int nsp1;                     /* ns plus 1 */
   double *wcalc;                /* calculation vector for weights */

   wcalc = dvector(nsta+1);
   while (1) {
      nsp1 = ns + 1;

      /* Load matrix for calculating kriging weights using only
         the desired stations (staflg = 1) */

      mm = -1;
      for (m = 0; m < nsta; m++) {
         if (staflg[m] == 1) {
            mm++;
            nn = -1;
            for (n = 0; n < nsta; n++) {
               if (staflg[n] == 1) {
                  nn++;
                  a[mm][nn] = ad[m][n];
               }
            }
            a[mm][ns] = a[ns][mm] = 1;
            a[mm][nsp1] = dgrid[l][m];
         }
      }
      a[ns][ns] = 0;
      a[ns][nsp1] = 1;
      n = nsp1;
/* Debug
fprintf(fpout, "\n\n%s %d, %s %d, %s %d:\n",
        "Coefficient matrix for grid cell",
        l+1, "day", j+1, "year", year[k]);
for (mm = 0; mm < n; mm++) {
   fprintf(fpout, "\n");
   for (nn = 0; nn <=n; nn++)
      fprintf(fpout, "%9.4f", a[mm][nn]);
}
   End debug */

      /* Solve linear system for kriging weights */
   
      if ((luret = lusolv(n, a, wcalc)) != 0) {
         if (icoord == 1)
            fprintf(fpout, "\n\n%s\n%s%d%s%5.2f%s%6.2f%s%6.0f\n\n%s\n",
                    "Indeterminate linear system ... ",
                     "   Grid cell ", l+1, ":  lat ", grid[l].north,
                     "   long ", grid[l].east, "   elev ", grid[l].elev*1000,
                     "Program terminating ...");
         else
            fprintf(fpout, "\n\n%s\n%s%d%s%10.2f%s%10.2f%s%6.0f\n\n%s\n",
                    "Indeterminate linear system ... ",
                     "   Grid cell ", l+1, ":  northing ", grid[l].north,
                     "   easting ", grid[l].east, "   elev ", grid[l].elev*1000,
                     "Program terminating ...");
         exit(0);
      }
/* Debug
fprintf(fpout, "\n\n%s %d, %s %d, %s %d:\n",
        "Kriging weights for grid cell",
        l+1, "day", j+1, "year", year[k]);
for (nn = 0; nn < ns; nn++)
   fprintf(fpout, "%9.4f", wcalc[nn]);
fprintf(fpout, "\n\n");
   End debug */

      /* Check for negative weights, throw out the most distant station with
         a negative weight, and recalculate weights until all are positive */

      elevsave = 0.0;
      mm = msave = -1;
      for (m = 0; m < nsta; m++) {
         if (staflg[m] == 1) {
            mm++;
            if (wcalc[mm] < 0.0) {
               if (sta[m].elev > elevsave) {
                  msave = m;
                  elevsave = sta[m].elev;
               }
            }
         }
      }
      if (msave >= 0) {
         staflg[msave] = 0;
         ns--;
      }
      else {
         mm = -1;
         for (m = 0; m < nsta; m++) {
            if (staflg[m] == 1) {
               mm++;
               w[m] = wcalc[mm];
            }
            else
               w[m] = 0.0;
         }
         break;
      }
   }
   free(wcalc);
}

