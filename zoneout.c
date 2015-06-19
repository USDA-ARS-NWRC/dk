/*
 *    zoneout.c
 *
 *    David Garen  11/2005
 *
 *    Compute and write out zonal mean values
 *
 *    Modification, 29 May 2007:
 *       Added flag to indicate all input values are zero.  This corresponds
 *       to modification in period2.c to write a line in zone output for every
 *       day, even if all values are zero (pertains to precipitation and SWE).
 *
 *    Modification, 18 July 2011:
 *       Changed output format from tabular to OMS-compatible csv format.
 *
 *    Modification, 28 November 2012:
 *       Changed order of zone output to be in numerical order, using
 *       zoneseq array to indicate array index for ordering
 */

#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void zoneout(iy, id, iz)
int iy;                          /* year */
int id;                          /* day (sequential number beginning Oct 1) */
int iz;                          /* zero flag (0 = all values are zero */
{
   void caldate();               /* julian day to calendar day conversion function */
   int day;                      /* day of month */
   int i, j;                     /* loop indexes */
   int month;                    /* calendar month number */

   /* Set all zonal mean values to zero */

   for (j = 0; j < nzone; j++)
      zone[j].mean = 0.0;

   /* Compute zonal means if input not all zero */

   if (iz != 0) {
      for (i = 0; i < ngrid; i++) {
         for (j = 0; j < nzone; j++) {
            if (grid[i].zone == zone[j].number)
               zone[j].mean += gprec[i];
         }
      }

      for (j = 0; j < nzone; j++)
         zone[j].mean /= zone[j].ncells;
   }

   /* Determine month and day for given water year julian day */

   caldate(iy, (id+1), &month, &day);

   /* Write year, month, and day */

   if (month >= 10)
      fprintf(fpzone, ",%d %d %d 0 0 0", (iy-1), month, day);
   else
      fprintf(fpzone, ",%d %d %d 0 0 0", iy, month, day);

   /* Write zonal means for day */

   for (j = 0; j < nzone; j++)
      fprintf(fpzone, ",%.2f", zone[zoneseq[j]].mean);
   fprintf(fpzone, "\n");

}
