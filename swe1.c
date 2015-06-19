/*
 *    swe1.c
 *
 *    David Garen  10/93, 12/93, 10/97
 *
 *    Compute snow water equivalent vs. elevation trend lines for
 *    dpp-day periods.
 *
 *    Modified handling of snow line and trend between first zero and
 *    first nonzero value.   DCG  28 October 1997
 *
 *    Modification 18 December 2012:
 *       Small changes in wording of output header lines (first lines
 *       written to fpout in code below)
 */

#include <stdio.h>

#include "dk_x.h"

void swe1()
{
   int i, j, jj, k, m, n, nn;    /* loop indexes */
   double *elev;                 /* station elevation vector for sorting */
   void indexx();                /* sorting function */
   int *inx;                     /* index vector */
   float sca;                    /* snow covered area as defined by snolin
                                    and grid */
   float sca_grid();             /* function to calculate sca by determining
                                    number of grid cells above snolin */

   /* Sort swe station elevations */

   elev = dvector(nsta);
   inx = ivector(nsta);
   for (i = 0; i < nsta; i++)
      elev[i] = sta[i].elev;
   indexx(elev, inx, nsta);

   /* Begin regression calculations */

   if (ireg == 1) {
      fprintf(fpout, "SWE-elevation regressions ");
      if (irmeth == 1)
         fprintf(fpout, "(least squares regression)");
      if (irmeth == 2)
         fprintf(fpout, "(least absolute deviations regression)");
      fprintf(fpout, "\n%s\n%s", "for individual years and periods:",
             "(Note:  Slopes based on elevation / 1000)");
   }
   for (k = 0; k < nyear; k++) {
      if (ireg == 1) {
         fprintf(fpout, "\n\n\n%s%d\n\n%s",
                 "YEAR ", year[k], "PERIOD  INTERCEPT      SLOPE");
         if (irmeth == 1)
            fprintf(fpout, "%s\n",
                    "        R       SE        T    N  SNOW LINE    SCA (%)");
         if (irmeth == 2)
            fprintf(fpout, "%s\n",
                    "      MAE    N  SNOW LINE    SCA (%)");
      }
      n = lastday[k] - firstday[k] + 1;
      nper = n / dpp;
      nperm1 = nper - 1;
      dppl = n - dpp * nperm1;
      nstop = dpp;
      for (m = 0; m < nper; m++) {
         if (m == nperm1)
            nstop = dppl;
         jj = dpp * m + firstday[k] - 1;

         /* Compute period swe totals for each station that has
            no missing data during the period */

         for (i = 0; i < nsta; i++) {
            adata[i] = 0;
            for (n = 0; n < nstop; n++) {
               j = jj + n;
               if (sta[i].data[j][k] < missing)
                  adata[i] += sta[i].data[j][k];
               else {
                  adata[i] = 99999;
                  break;
               }
            }
         }

         /* Of the stations that have no missing data during the period,
            determine the number of days where at least one station
            has nonzero swe */

         nn = nstop;
         for (n = 0; n < nstop; n++) {
            j = jj + n;
            izero = 1;
            for (i = 0; i < nsta; i++) {
               if (adata[i] <= 99998) {
                  if (sta[i].data[j][k] > 0.001) {
                     izero = 0;
                     break;
                  }
               }
            }
            if (izero == 1)
               nn--;
         }

         /* Compute average daily swe by dividing period swe total
            by number of days where at least one station had nonzero swe */

         for (i = 0; i < nsta; i++)
            if (adata[i] <= 99998 && nn > 0)
               adata[i] /= nn;

         /* In order of elevation, find first nonzero swe and see if
            there is at least one other nonzero value above it;
            if there is, perform the regression */

         iswehz[m][k] = isweln[m][k] = -1;
         for (i = 0; i < nsta; i++) {
            if (adata[inx[i]] <= 99998 && adata[inx[i]] > 0.0001)
               break;
            iswehz[m][k] = inx[i];
         }
         isweln[m][k] = i;
         n = 0;
         for (i = isweln[m][k]; i < nsta; i++)
            if (adata[inx[i]] <= 99998 && adata[inx[i]] > 0.0001)
               n++;
         if (n > 0)
            isweln[m][k] = inx[isweln[m][k]];
/* Debug
printf("\nPeriod %d:  Highest zero = %d     Lowest nonzero = %d", m+1,
       iswehz[m][k], isweln[m][k]);
   End debug */
         if (n >= 2) {

            /* Calculate line between highest zero value and lowest
               nonzero value */

            if (iswehz[m][k] >= 0) {
               b12[m][k] = adata[isweln[m][k]] /
                           (sta[isweln[m][k]].elev - sta[iswehz[m][k]].elev);
               b02[m][k] = -b12[m][k] * sta[iswehz[m][k]].elev;
            }
            else
               b12[m][k] = b02[m][k] = 0;

            /* Load data arrays */

            n = -1;
            for (i = isweln[m][k]; i < nsta; i++) {
               if (adata[inx[i]] <= 99998) {
                  n++;
                  x[n] = sta[inx[i]].elev;
                  y[n] = adata[inx[i]];
               }
            }
/* Debug
if (m == 8 && year[k] == 69) {
   fprintf(fpout, "\n\n\nRegression data:\n\n");
   if (type == 1)
      fprintf(fpout, "nn = %d\n", nn);
   for (i = 0; i <= n; i++)
      fprintf(fpout, "%15.8f%15.8f\n", x[i], y[i]);
}
   End debug */

            /* Calculate trend lines */

            if (irmeth == 1)
               ret = sreg(x, y, &b0dum, &b1dum, &r, &se, &t, n+1);
            if (irmeth == 2)
               ret = medfit(x, y, &b0dum, &b1dum, &mae, n+1);
            if (ret == 0) {
               if (b1dum < 0.0) {
/* Debug
printf("\nPeriod %d: medfit slope = %f", m+1, b1dum);
   End debug */
                  b0[m][k] = b1[m][k] = 0;
                  if (isweln[m][k] > 0)
                     snolin[m][k] = sta[inx[iswehz[m][k]]].elev;
                  else
                     snolin[m][k] = 0;
                  if (irmeth == 1)
                     r = se = t = 0;
                  if (irmeth == 2)
                     mae = 0;
               }
               else {
                  b0[m][k] = (float) b0dum;
                  b1[m][k] = (float) b1dum;
                  snolin[m][k] = -b0[m][k] / b1[m][k];
                  if (isweln[m][k] > 0)
                     if (snolin[m][k] < sta[inx[iswehz[m][k]]].elev)
                        snolin[m][k] = sta[inx[iswehz[m][k]]].elev;
               }
               sca = sca_grid(snolin[m][k]);
               if (ireg == 1) {
                  if (irmeth == 1)
                     fprintf(fpout,
                             "\n%6d%11.4f%11.4f%9.3f%9.4f%9.3f%5d%11.0f%11.1f",
                             m+1, b0[m][k], b1[m][k], r, se, t, n+1,
                             snolin[m][k]*1000, sca);
                  if (irmeth == 2)
                     fprintf(fpout,
                             "\n%6d%11.4f%11.4f%9.4f%5d%11.0f%11.1f",
                             m+1, b0[m][k], b1[m][k], mae, n+1,
                             snolin[m][k]*1000, sca);
               }

               /* Compute residuals */

               for (i = 0; i < nsta; i++) {
                  if (sta[i].elev <= snolin[m][k])
                     dum = 0;
                  else
                     dum = b0[m][k] + b1[m][k] * sta[i].elev;
                  for (n = 0; n < nstop; n++) {
                     j = jj + n;
                     if (map[j][k] > missing && sta[i].data[j][k] < missing)
                        sta[i].data[j][k] -= dum;
                  }
               }
            }
            else if (ret == 1 && ireg == 1)
               fprintf(fpout, "\n%6d  %s", m+1,
                       "No regression possible -- all x data are equal.");
            else if (ret == 3 && ireg == 1)
               fprintf(fpout, "\n%6d  %s%s", m+1,
                       "No regression possible -- ",
                        "all x and y data are equal.");
         }
         else {

            /* Can't perform the regression; if there are two or more
               valid data values (including zero), consider the trend line
               to be zero and go ahead and allow spatial interpolation */

            n = 0;
            for (i = 0; i < nsta; i++)
               if (adata[i] <= 99998)
                  n++;
            if (n >= 2) {
               b1[m][k] = b0[m][k] = sca = 0;
               snolin[m][k] = sta[iswehz[m][k]].elev;
               sca = sca_grid(snolin[m][k]);
/* Debug
printf("\nPeriod %d: Highest zero = %d     Snow line = %f", m+1, iswehz[m][k],
       snolin[m][k]);
   End debug */
               if (irmeth == 1)
                  r = se = t = 0;
               if (irmeth == 2)
                  mae = 0;
               if (ireg == 1) {
                  if (irmeth == 1)
                     fprintf(fpout,
                             "\n%6d%11.4f%11.4f%9.3f%9.4f%9.3f%5d%11.0f%11.1f",
                             m+1, b0[m][k], b1[m][k], r, se, t, n,
                             snolin[m][k]*1000, sca);
                  if (irmeth == 2)
                     fprintf(fpout,
                             "\n%6d%11.4f%11.4f%9.4f%5d%11.0f%11.1f",
                             m+1, b0[m][k], b1[m][k], mae, n, snolin[m][k]*1000,
                             sca);
               }
            }
            else {
               if (ireg == 1)
                  fprintf(fpout, "\n%6d  %s", m+1,
                     "No estimation possible -- fewer than two data values.");
            }
         }
      }
   }
   if (ireg == 1)
      fprintf(fpout, "\n\n\n");
}
