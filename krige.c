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
#include <malloc/malloc.h>

#include "dk_x.h"

double *krige(l, nsta, ad, dgrid, elevations, w)
int l;                           /* grid index */
int nsta;                          /* number of stations used */
float **ad;                      /* matrix of distances between prec/temp
                                    stations for computing kriging weights */
float **dgrid;                   /* matrix of distances between grid cells
                                    and prec/temp stations */
int *elevations;				 /* vector of station elevations */
double *w;                    /* kriging weights */
{
	float elevsave;               /* stored value of station elevation */
	int m, mm, n, nn, i, j;             /* loop indexes */
	int msave;                    /* stored value of m index */
	int nsp1;                     /* ns plus 1 */
	double *wcalc;                /* calculation vector for weights */
	int ns;					 	 /* number of stations */
	int luret;                    /* return value from lusolv() */
	int *staflg;				 	 /* station use flags*/
	float temp;						 /* temporary variable */
	int itemp;						 /* temporary variable */
	float *dist;					 /* sorted distance */
	int *idx;					 /* index sorted distance */
	double **a;                   /* data matrix for solving for kriging
                                       weights (input to m_inv()) */
//	double *w;                    /* kriging weights */
//	double w[nsta+1];

	//   nsta = ns;

	// find the N closest stations
	dist = vector(nsta);
	idx = ivector(nsta);
	for (i = 0; i < nsta; ++i){
		dist[i] = dgrid[l][i];
		idx[i] = i;
	}
	for (i = 0; i < nsta; ++i){
		for (j = i + 1; j < nsta; ++j){
			if (dist[i] > dist[j]) {
				// sort the distance
				temp = dist[i];
				dist[i] = dist[j];
				dist[j] = temp;

				// re-sort the index
				itemp = idx[i];
				idx[i] = idx[j];
				idx[j] = itemp;
			}
		}
	}


	// set station use flags
	ns = 0;
	staflg = ivector(nsta);
	for (m = 0; m < nsta; m++) {
		staflg[idx[m]] = 1;
		ns++;
	}
	//   for (i = 0; i < nsta; ++i){
	//	   printf("%f - %i - %i\n",dist[i],idx[i],staflg[i]);
	//   }
	//   exit(0);

	a = dmatrix(nsta+1, nsta+2);
//	w = dvector(nsta+1);

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

		/* Check for negative weights, throw out the most distant station by elevation with
         a negative weight, and recalculate weights until all are positive */

		elevsave = 0.0;
		mm = msave = -1;
		for (m = 0; m < nsta; m++) {
			if (staflg[m] == 1) {
				mm++;
				if (wcalc[mm] < 0.0) {
					if (elevations[m] > elevsave) {
						msave = m;
						elevsave = elevations[m];
					}
				}
			}
		}
		if (msave >= 0) {
			staflg[msave] = 0; // set station use flag to zero for furthest station
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

	/* clean up 1D arrays*/
	free(dist);
	free(idx);
	free(staflg);
	free(wcalc);

	/* clean up 2D arrays*/
	for (m = 0; m < nsta+1; m++) {
		free(a[m]);
	}
	free(a);

//	return w;
}

