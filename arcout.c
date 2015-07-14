/*
 *    arcout.c
 *
 *    David Garen  1/2000
 *
 *
 *    Write daily grids to output files in ARC/INFO format
 *
 *    File names are as follows:  parameter (prc, tmp, swe), year, day . arc
 *       Examples:  prc89058.arc, tmp91255.arc, swe93121.arc
 *
 *    Modification 28 March 2005:
 *       Changed file naming convention.  Example:  prc_2004_004.asc
 *
 *    Modification 6 November 2006:
 *       Added day fraction to file naming convention.
 *       Example:  prc_2004_004_125.asc (125 corresponds to hour 3)
 *
 *    Modification for Version 4.7:
 *       Changed file naming convention to reflect generic time periods
 *       instead of days, and removed day fraction.
 *       Example:  prc_2004_6358.asc
 */

#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void arcout(iy, ip)
int iy;                          /* year */
int ip;                          /* period (sequential number beginning Oct 1) */
{
   char buf[6];                  /* buffer for file name building */
   FILE *fparc;                  /* output file pointer */
   int i, j;                     /* loop indexes */
   int k;                        /* grid value counter */
   char outfile[21];             /* output file name */

   /* Build output file name and open file */

   if (type == 1)
      strcpy(outfile, "prc_");
   else if (type == 2)
      strcpy(outfile, "tmp_");
   else if (type == 3)
      strcpy(outfile, "swe_");
   else
      strcpy(outfile, "dat_");
   sprintf(buf, "%04d_", iy);
   strcat(outfile, buf);
   sprintf(buf, "%04d", (ip+1));
   strcat(outfile, buf);
/* sprintf(buf, "%03d", dayfrac);
   strcat(outfile, buf); */
   strcat(outfile, ".asc");
   if ((fparc = fopen(outfile, "w")) == NULL) {
      printf("\n\nError opening file %s.\n", outfile);
      return;
   }

   /* Write header information */

   fprintf(fparc, "ncols %d\n", arc.cols);
   fprintf(fparc, "nrows %d\n", arc.rows);
   fprintf(fparc, "xllcorner %.2f\n", arc.xll);
   fprintf(fparc, "yllcorner %.2f\n", arc.yll);
   fprintf(fparc, "cellsize %.2f\n", arc.cell);
   fprintf(fparc, "nodata_value %.0f\n", (arc.nodata-0.1));

   /* Write grid values */

   k = -1;
   for (i = 0; i < arc.rows; i++) {
      for (j = 0; j < arc.cols; j++) {
         k++;
         if (grid[k].use == 1) {
            if (igridpr == 1)
               fprintf(fparc, "%.1f ", gprec[k]);
            else if (igridpr == 2)
               fprintf(fparc, "%.0f ", gprec[k]);
            else if (igridpr == 3)
               fprintf(fparc, "%.0f ", (gprec[k]*10));
         }
         else
            fprintf(fparc, "%.0f ", (arc.nodata-0.1));
      }
      fprintf(fparc, "\n");
   }
   fclose(fparc);
}
