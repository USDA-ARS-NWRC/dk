/*
 *    period1.c
 *
 *    David Garen  9/92
 *
 *    Compute precipitation, temperature, or other parameter vs. elevation
 *    trend lines for dpp-day periods.
 *
 *    Modification, 8 November 2006:
 *       Created "other" parameter type category (type = 4) to allow detrending
 *       without screening or restrictions
 *
 *    Modification 18 December 2012:
 *       Small changes in wording of output header lines (first lines
 *       written to fpout in code below)
 */

#include <stdio.h>

#include "dk_m.h"
#include "dk_x.h"

void period1()
{
	int i, j, jj, k, m, n, nn;    /* loop indexes */
	float interp();               /* accumulated precip interpolation function */
	float intval[MSTA];           /* interpolated accumulated precip value */

	for (i = 0; i < nsta; i++)
		intval[i] = -1;

	/* Compute regressions and residuals */

	if (ireg == 1) {
		if (type == 1)
			fprintf(fpout, "\n\n\nPrecipitation");
		else if (type == 2)
			fprintf(fpout, "\n\n\nTemperature");
		else
			fprintf(fpout, "\n\n\nParameter");
		fprintf(fpout, "-elevation regressions ");
		if (irmeth == 1)
			fprintf(fpout, "(least squares regression)");
		if (irmeth == 2)
			fprintf(fpout, "(least absolute deviations regression)");
		fprintf(fpout, "\n%s\n%s", "for individual years and aggregation periods:",
				"(Note:  Slopes based on elevation / 1000)");
	}

	for (k = 0; k < nyear; k++) {
		if (ireg == 1) {
			fprintf(fpout, "\n\n\n%s%d\n\n%s\n%s",
					"YEAR ", year[k], "AGG.", "PERIOD  INTERCEPT      SLOPE");
			if (irmeth == 1)
				fprintf(fpout, "        R       SE        T    N\n");
			if (irmeth == 2)
				fprintf(fpout, "      MAE    N\n");
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

			/* Compute period totals for each station that has
            no missing data during the period */

			for (i = 0; i < nsta; i++) {
				adata[i] = 0;
				for (n = 0; n < nstop; n++) {
					j = jj + n;
					if (type == 1 && (sta[i].data[j][k] > accum
							&&  sta[i].data[j][k] < replace)) {
						if (intval[i] <= 0)
							intval[i] = interp(i, j, k);
						adata[i] += intval[i];
					}
					else if (type == 1 && (sta[i].data[j][k] > replace
							&&  sta[i].data[j][k] < missing)) {
						adata[i] += intval[i];
						intval[i] = -1;
					}
					else if (sta[i].data[j][k] < missing)
						adata[i] += sta[i].data[j][k];
					else {
						adata[i] = 99999;
						break;
					}
				}
			}

			/* Of the stations that have no missing data during the period,
            determine the number of days where at least one station
            has prec (type = 1) */

			if (type == 1) {
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
			}

			/* Compute average daily prec by dividing period prec total
            by number of days where at least one station had prec;
            compute average daily temp or other parameter by dividing
            period total by number of days in period */

			for (i = 0; i < nsta; i++) {
				if (type == 1) {
					if (adata[i] <= 99998 && nn > 0)
						adata[i] /= nn;
				}
				else {
					if (adata[i] <= 99998)
						adata[i] /= nstop;
				}
			}
			/* Debug
fprintf(fpout, "\n\nadata:  ");
for (i = 0; i < nsta; i++)
   fprintf(fpout, "%10.6f", adata[i]);
fprintf(fpout, "\nnn = %d", nn);
   End debug */

			/* Load data arrays */

			n = -1;
			for (i = 0; i < nsta; i++) {
				if (adata[i] <= 99998) {
					if (type == 1){
						if (adata[i] > 0){ // if the precip is greater than zero keep
							n++;
							x[n] = sta[i].elev;
							y[n] = adata[i];
						}
					}
					else {
						n++;
						x[n] = sta[i].elev;
						y[n] = adata[i];
					}
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

			/* If only one point, make slope and int zero */
			if (n == 0) {
				b0[m][k] = b1[m][k] = 0;
			}

			/* if 2 or more points, fit line, constrain */
			if (n > 0) {
				if (irmeth == 1)
					ret = sreg(x, y, &b0dum, &b1dum, &r, &se, &t, n+1);
				if (irmeth == 2)
					ret = medfit(x, y, &b0dum, &b1dum, &mae, n+1);
				if (ret == 0) {
					if ((type == 1 && b1dum < 0.0)
							|| (type == 2 && b1dum > 0.0)) {
						b0[m][k] = b1[m][k] = 0;
						if (irmeth == 1)
							r = se = t = 0;
						if (irmeth == 2)
							mae = 0;
					}
					else {
						b0[m][k] = (float) b0dum;
						b1[m][k] = (float) b1dum;
					}
					if (ireg == 1) {
						if (irmeth == 1)
							fprintf(fpout, "\n%6d%11.4f%11.4f%9.3f%9.4f%9.3f%5d",
									m+1, b0[m][k], b1[m][k], r, se, t, n+1);
						if (irmeth == 2)
							fprintf(fpout, "\n%6d%11.4f%11.4f%9.4f%5d",
									m+1, b0[m][k], b1[m][k], mae, n+1);
					}

					/* Compute residuals */

					if (type != 1) {
						for (i = 0; i < nsta; i++) {
							dum = b0[m][k] + b1[m][k] * sta[i].elev;
							for (n = 0; n < nstop; n++) {
								j = jj + n;
								if (map[j][k] > missing && sta[i].data[j][k] < accum)
									sta[i].data[j][k] -= dum;
							}
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
			else
				if (ireg == 1)
					fprintf(fpout, "\n%6d  %s", m+1,
							"No regression possible -- insufficient data pairs.");
		}
	}
	if (ireg == 1)
		fprintf(fpout, "\n\n\n");
}

