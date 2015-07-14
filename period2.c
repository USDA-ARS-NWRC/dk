/*
 *    period2.c
 *
 *    David Garen  9/92
 *
 *    Compute MAP/MAT based on dpp-day periods.
 */

#include <stdio.h>

#include "dk_x.h"

int netcdf_create();
int netcdf_write();
int netcdf_close();

void period2()
{
	int i, j, jj, k, l, m, n;  /* loop indexes */
	int ns;                    /* number of stations with data */
	int *staflg;				 	 /* station use flags*/
	int *ncid;		/* file id for netcdf file */

	// set station use flags
	staflg = ivector(nsta);
	for (m = 0; m < nsta; m++)
		staflg[m] = 1;

	/* Year loop */

	for (k = 0; k < nyear; k++) {
		n = lastday[k] - firstday[k] + 1;
		nper = n / dpp;
		nperm1 = nper - 1;
		dppl = n - dpp * nperm1;
		nstop = dpp;

		/* Create the netcdf file if wanted */
		if (iout == 5)
			netcdf_create(year[k], &grid, arc.cols, arc.rows, &ncid);

		/* Period loop */

		for (m = 0; m < nper; m++) {
			if (m == nperm1)
				nstop = dppl;
			jj = dpp * m + firstday[k] - 1;

			/* Process all days that have valid detrending coefficients */

			if (b0[m][k] <= 99998 && b1[m][k] <= 99998) {

				/* Day loop */

				for (n = 0; n < nstop; n++) {
					j = jj + n;
					/* @@@ Progress "..." still needed?
               printf("   %s%d%s%d%s%d%s\n",
                      "Processing year ", year[k], ", period ", m+1,
                      ", day ", j+1, " ...");
					 */

					/* Process days that have not been set to zero because
                  all stations have zero values (precipitation, SWE) */

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
									/* Debug
fprintf(fpout, "\nRevised weights:  Grid point %d, Year %d, Period %d\n",
        l+1, year[k], m+1);
for (i = 0; i < nsta; i++)
   fprintf(fpout, "%8.4f", w[i]);
fprintf(fpout, "\n");
   End debug */
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

								gprec[l] += (b0[m][k] + b1[m][k] * grid[l].elev);
								/* Experiment -- scale factor
                        gprec[l] /= 100.;
   End experiment */

								/* Set grid prec values to zero if estimate is
                           less than zero */

								if (gprec[l] < 0 && type == 1)
									gprec[l] = 0;

								/* Add grid prec/temp to basin sum */

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
							zoneout(year[k], j, 1);

						/* Compute MAP/MAT for day */

						if (imask == 0)
							map[j][k] = dum / ngriduse;
						else
							map[j][k] = dum / nmask;
					}
					else if (iout == 5) {
						/* no data to calculate, create a zero image mainly to fill precip */
						for (l = 0; l < ngrid; l++)
							if (grid[l].use == 1)
								gprec[l] = 0;
					}
					else {

						/* For zone output, write a line of output anyway
                     if all stations have zero values */

						if (izone == 1)
							zoneout(year[k], j, 0);
					}
				}
				/* If requested, write out grid in NETCDF format */

				if (iout == 5 && j >= igridout1 && j <= igridout2)
					netcdf_write(&ncid, j, gprec, arc.cols, arc.rows);
			}
		}

		if (iout == 5)
			netcdf_close(&ncid);

	}
}
