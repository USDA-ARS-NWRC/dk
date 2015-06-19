/*
 *    storm1.c
 *
 *    David Garen  9/92, 8/95, 5/97
 *
 *    Compute precipitation vs. elevation trend lines and residuals for
 *    storms.  The following rules are applied:
 *
 *    1)  In general, a storm is defined as a consecutive sequence of days
 *    during which at least three stations have at least 5 mm of precipitation.
 *
 *    2)  If only one or two stations have precipitation on any day, no
 *    detrending is done, and the day is not combined with any others.
 *
 *    3)  If three or more stations have precipitation, but all stations have
 *    5 mm or less, no detrending is done, and the day is not combined with
 *    any others.
 *
 *    Modification 18 December 2012:
 *       Small changes in wording of output header lines (first lines
 *       written to fpout in code below)
 */

#include <stdio.h>

#include "dk_x.h"

void storm1()
{
   int i, j, jj, k, kk, n;       /* loop indexes */
   int k5;                       /* flag to indicate presence of precipitation
                                    > 5 mm on a given day */
   int ksta;                     /* number of stations with precipitation on
                                    a given day */

   if (ireg == 1) {
      fprintf(fpout, "Precipitation-elevation regressions ");
      if (irmeth == 1)
         fprintf(fpout, "(least squares regression)");
      if (irmeth == 2)
         fprintf(fpout, "(least absolute deviations regression)");
      fprintf(fpout, "\n%s\n%s", "for individual storms:",
              "(Note:  Slopes based on elevation / 1000)");
      fprintf(fpout, "\n\n\n%s",
              "STORM  YEAR  DAY  LENGTH  INTERCEPT    SLOPE");
      if (irmeth == 1)
         fprintf(fpout, "       R      SE       T   N\n");
      if (irmeth == 2)
         fprintf(fpout, "     MAE   N\n");
   }
   nstorm = 0;
   storm[0].dstart = -1;
   storm[0].slen = 0;
   for (k = 0; k < nyear; k++) {
      for (j = (firstday[k]-1); j < lastday[k]; j++) {
         izero = 1;
         ksta = k5 = 0;
         for (i = 0; i < nsta; i++) {
            if (sta[i].data[j][k] > 0.01 && sta[i].data[j][k] < missing)
               ksta++;
            if (sta[i].data[j][k] > 5.0 && sta[i].data[j][k] < missing)
               k5 = 1;
         }
         if (ksta >= 3 && k5 == 1) {
            if (storm[nstorm].dstart == -1) {
               storm[nstorm].dstart = j;
               storm[nstorm].ystart = k;
            }
            izero = 0;
            storm[nstorm].slen++;
         }
         else if (storm[nstorm].dstart > -1) {
            izero = 1;
            j--;
         }
         else if (ksta >= 1) {
            if (storm[nstorm].dstart == -1) {
               storm[nstorm].dstart = j;
               storm[nstorm].ystart = k;
            }
            storm[nstorm].slen++;
            b0[nstorm][0] = b1[nstorm][0] = 0.0;
/*          printf("\nProcessing storm %d:  year = %d, day = %d\n",
                   nstorm+1, year[storm[nstorm].ystart],
                   storm[nstorm].dstart+1); */
            if (ireg == 1)
               fprintf(fpout,
                       "\n%5d%6d%5d%8d%11.4f%9.4f",
                       nstorm+1, year[storm[nstorm].ystart],
                       storm[nstorm].dstart+1, storm[nstorm].slen,
                       b0[nstorm][0], b1[nstorm][0]);
            nstorm++;
            storm[nstorm].dstart = -1;
            storm[nstorm].slen = 0;
         }
         if (j == (lastday[k]-1) && k == (nyear-1))
            izero = 1;
         if (izero == 1 && storm[nstorm].dstart > -1) {
/*          printf("\nProcessing storm %d:  year = %d, day = %d\n",
                   nstorm+1, year[storm[nstorm].ystart],
                   storm[nstorm].dstart+1); */

            /* Compute period prec totals for each station that has
               no missing data during the storm */

            for (i = 0; i < nsta; i++) {
               adata[i] = 0;
               jj = storm[nstorm].dstart;
               kk = storm[nstorm].ystart;
               dstop = lastday[kk] - 1;
               for (n = 0; n < storm[nstorm].slen; n++) {
                  if (sta[i].data[jj][kk] < accum)
                     adata[i] += sta[i].data[jj][kk];
                  else if (sta[i].data[jj][kk] > missing) {
                     adata[i] = 99999;
                     break;
                  }
                  jj++;
                  if (jj > dstop) {
                     jj = 0;
                     kk++;
                  }
               }
            }
/* Debug
fprintf(fpout, "\n\nadata :  ");
for (i = 0; i < nsta; i++)
   fprintf(fpout, "%10.6f", adata[i]);
   End debug */

            /* Compute average daily prec by dividing storm total
               by the length of the storm (in days) */

            for (i = 0; i < nsta; i++)
               if (adata[i] <= 99998)
                  adata[i] /= storm[nstorm].slen;

            /* Load data arrays */

            n = -1;
            for (i = 0; i < nsta; i++) {
               if (adata[i] <= 99998) {
                  n++;
                  x[n] = sta[i].elev;
                  y[n] = adata[i];
               }
            }
/* Debug
fprintf(fpout, "\nY data:  ");
for (i = 0; i <= n; i++)
   fprintf(fpout, "%10.6f", y[i]);
   End debug */

            /* Calculate trend lines */

            if (n > 0) {
               if (irmeth == 1)
                  ret = sreg(x, y, &b0dum, &b1dum, &r, &se, &t, n+1);
               if (irmeth == 2) {
                  ret = medfit(x, y, &b0dum, &b1dum, &mae, n+1);
/* Debug
fprintf(fpout, "\nReturning from medfit() ...\nintercept = %f   slope = %f   mae = %f\n",
        b0dum, b1dum, mae);
   End debug */
               }
               if (ret == 0) {
                  if (b1dum < 0.0) {
                     b0[nstorm][0] = b1[nstorm][0] = 0;
                     if (irmeth == 1)
                        r = se = t = 0;
                     if (irmeth == 2)
                        mae = 0;
                  }
                  else {
                     b0[nstorm][0] = (float) b0dum;
                     b1[nstorm][0] = (float) b1dum;
                  }
                  if (ireg == 1) {
                     if (irmeth == 1)
                        fprintf(fpout,
                                "\n%5d%6d%5d%8d%11.4f%9.4f%8.3f%8.3f%8.3f%4d",
                                nstorm+1, year[storm[nstorm].ystart],
                                storm[nstorm].dstart+1, storm[nstorm].slen,
                                b0[nstorm][0], b1[nstorm][0], r, se, t, n+1);
                     if (irmeth == 2)
                        fprintf(fpout,
                                "\n%5d%6d%5d%8d%11.4f%9.4f%8.3f%4d",
                                nstorm+1, year[storm[nstorm].ystart],
                                storm[nstorm].dstart+1, storm[nstorm].slen,
                                b0[nstorm][0], b1[nstorm][0], mae, n+1);
                  }

                  /* Compute residuals */

                  for (i = 0; i < nsta; i++) {
                     dum = b0[nstorm][0] + b1[nstorm][0] * sta[i].elev;
                     jj = storm[nstorm].dstart;
                     kk = storm[nstorm].ystart;
                     dstop = lastday[kk] - 1;
                     for (n = 0; n < storm[nstorm].slen; n++) {
                        if (sta[i].data[jj][kk] < accum)
                           sta[i].data[jj][kk] -= dum;
                        jj++;
                        if (jj > dstop) {
                           jj = 0;
                           kk++;
                        }
                     }
                  }
               }
               else if (ret == 1 && ireg == 1)
                  fprintf(fpout, "\n%5d%6d%5d%8d  %s",
                          nstorm+1, year[storm[nstorm].ystart],
                          storm[nstorm].dstart+1, storm[nstorm].slen,
                          "No regression possible -- all x data are equal.");
               else if (ret == 3 && ireg == 1)
                  fprintf(fpout, "\n%5d%6d%5d%8d  %s%s",
                          nstorm+1, year[storm[nstorm].ystart],
                          storm[nstorm].dstart+1, storm[nstorm].slen,
                          "No regression possible -- ",
                          "all x and y data are equal.");
            }
            else {
               if (ireg == 1)
                  fprintf(fpout, "\n%5d%6d%5d%8d  %s",
                          nstorm+1, year[storm[nstorm].ystart],
                          storm[nstorm].dstart+1, storm[nstorm].slen,
                          "No regression possible -- insufficient data pairs.");
            }
            nstorm++;
            storm[nstorm].dstart = -1;
            storm[nstorm].slen = 0;
         }
      }
   }
   if (ireg == 1)
      fprintf(fpout, "\n\n\n");
}

