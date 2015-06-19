/*
 *    grassout.c
 *
 *    David Garen  8/95
 *
 *    Write daily grids to output files in GRASS format
 *
 *    File names are as follows:  parameter (prc, tmp, swe), year, day . grs
 *       Examples:  prc89058.grs, tmp91255.grs, swe93121.grs
 *
 *    Modification 21 September 2006 (to match changes in arcout.c
 *    and ipwout.c made 28 March 2005):
 *       Changed file naming convention.  Example:  prc_2004_004.grs
 *
 *    Modification 6 November 2006:
 *       Added day fraction to file naming convention.
 *       Example:  prc_2004_004_125.grs (125 corresponds to hour 3)
 *
 *    Modification for Version 4.7:
 *       Changed file naming convention to reflect generic time periods
 *       instead of days, and removed day fraction.
 *       Example:  prc_2004_6358.grs
 */

#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void grassout(iy, ip)
int iy;                          /* year */
int ip;                          /* period (sequential number beginning Oct 1) */
{
   char buf[6];                  /* buffer for file name building */
   FILE *fpgrs;                  /* output file pointer */
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
   strcat(outfile, ".grs");
   if ((fpgrs = fopen(outfile, "w")) == NULL) {
      printf("\n\nError opening file %s.\n", outfile);
      return;
   }

   /* Write header information */

   fprintf(fpgrs, "north: %.2f\n", grass.north);
   fprintf(fpgrs, "south: %.2f\n", grass.south);
   fprintf(fpgrs, "east: %.2f\n", grass.east);
   fprintf(fpgrs, "west: %.2f\n", grass.west);
   fprintf(fpgrs, "rows: %d\n", grass.rows);
   fprintf(fpgrs, "cols: %d\n", grass.cols);

   /* Write grid values */

   k = -1;
   for (i = 0; i < grass.rows; i++) {
      for (j = 0; j < grass.cols; j++) {
         if (igridpr == 1)
            fprintf(fpgrs, "%.1f ", gprec[++k]);
         else if (igridpr == 2)
            fprintf(fpgrs, "%.0f ", gprec[++k]);
         else if (igridpr == 3)
            fprintf(fpgrs, "%.0f ", (gprec[++k]*10));
      }
      fprintf(fpgrs, "\n");
   }
   fclose(fpgrs);
}
