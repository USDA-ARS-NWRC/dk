/*
 *    readdata.c
 *
 *    David Garen   6/92
 *
 *    Read input data for DK program
 *
 *    Modified for Version 4.7 by adding variable mtper (to replace 366)
 *    and removing dayfrac
 */

#include "malloc.h"
#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void readdata()
{
   double atof();                /* ascii-to-float function */
   char buf[4];                  /* buffer for parsing lat. and lng. */
   char clat[7];                 /* character version of latitude */
   char clng[8];                 /* character version of longitude */
   float decmin;                 /* decimal minutes */
   float decsec;                 /* decimal seconds */
   int i, j, k;                  /* loop indexes */
   int iyear;                    /* year -- temporary variable for reading */
   int len;                      /* string length */


   /* Read number of stations, number of years, first and last days for each year */
   /* (Delete day fraction) */

   fscanf(fpin1, "%d%d", &nsta, &nyear);
   year = ivector(nyear);
   firstday = ivector(nyear);
   lastday = ivector(nyear);

   for (k = 0; k < nyear; k++)
      fscanf(fpin1, "%d%d", &firstday[k], &lastday[k]);

/* fscanf(fpin1, "%d", &dayfrac); */

   /* Read station i.d., elevation, northing (or latitude),
      and easting (or longitude), and allocate array space for data */

   for (i = 0; i < nsta; i++) {
      if (icoord == 1)
         fscanf(fpin1, "%s%f%s%s", (char*) &sta[i].id, &sta[i].elev, clat, clng);
      else
         fscanf(fpin1, "%s%f%f%f", (char*) &sta[i].id, &sta[i].elev, &sta[i].north,
                &sta[i].east);
      sta[i].elev /= 1000;
      sta[i].data = matrix(mtper, nyear);
      for (j = 0; j < mtper; j++)
         for (k = 0; k < nyear; k++)
            sta[i].data[j][k] = missing;

      if (icoord == 1) {

         /* Convert latitude into decimal degrees */

         len = strlen(clat);
         strncpy(buf, &clat[len-2], 2);
         buf[2] = '\0';
         decsec = (float) (atof(buf) / 3600);
         strncpy(buf, &clat[len-4], 2);
         buf[2] = '\0';
         decmin = (float) (atof(buf) / 60);
         strncpy(buf, clat, len-4);
         buf[len-4] = '\0';
         sta[i].north = (float) (atof(buf) + decmin + decsec);

         /* Convert longitude into decimal degrees */

         len = strlen(clng);
         strncpy(buf, &clng[len-2], 2);
         buf[2] = '\0';
         decsec = (float) (atof(buf) / 3600);
         strncpy(buf, &clng[len-4], 2);
         buf[2] = '\0';
         decmin = (float) (atof(buf) / 60);
         strncpy(buf, clng, len-4);
         buf[len-4] = '\0';
         sta[i].east = (float) (atof(buf) + decmin + decsec);
      }
   }

   /* Read column format data */

   k = 0;
   fscanf(fpin1, "%d%d", &year[k], &j);
   j--;
   for (i = 0; i < nsta; i++)
      fscanf(fpin1, "%f", &sta[i].data[j][k]);
   while (1) {
      if (fscanf(fpin1, "%d%d", &iyear, &j) == EOF)
         break;
      if (iyear != year[k])
         year[++k] = iyear;
      j--;
      for (i = 0; i < nsta; i++)
         fscanf(fpin1, "%f", &sta[i].data[j][k]);
   }

   fclose(fpin1);


/* Experiment -- scaling factor
for (i = 0; i < nsta; i++)
   for (j = 0; j < mtper; j++)
      for (k = 0; k < nyear; k++)
         sta[i].data[j][k] *= 100.;
   End experiment */

}
