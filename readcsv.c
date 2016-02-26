/*
 *    readcsv.c
 *
 *    David Garen   8/2011, 10/2012
 *
 *    Read input data for DK program in OMS-compatible csv format
 *
 *    Modification, 28 November 2012:
 *       Added reading of dataname (precip, tmax, or tmin) from
 *       data file header line and changed keywords for easting
 *       and northing
 *
 *    Modification, 29 November 2012:
 *       Added variables yr_start, mo_start, dy_start, yr_end,
 *       mo_end, dy_end to save starting and ending dates of data
 *       and pass back to dk.c to put in zone output file header
 *
 *    Modification, 5 December 2012:
 *       Added reading missing data value from input (val_miss)
 *       and setting value in data matrix only if non-missing
 *       (otherwise leave at DK's internal missing value)
 *
 *    Modification, 10 July 2013:
 *       Removed leading space in reading "date_start", "date_end",
 *       and "missing_value".  Apparently this space was in early
 *       versions of the data csv file, but it is no longer there.
 */

#include <malloc/malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dk_x.h"

void readcsv()
{
   double atof();                /* ascii-to-float function */
   char buf[51];                 /* buffer for reading data fields */
   int fwy;                      /* final water year */
   int fwyjd;                    /* final water year julian day */
   int i, j, k;                  /* loop indexes */
   int ieast = 0;                /* flag indicating eastings found */
   int ielev = 0;                /* flag indicating elevations found */
   int iend = 0;                 /* flag indicating ending date found */
   int imiss = 0;                /* flag indicating missing data value found */
   int inorth = 0;               /* flag indicating northings found */
   int istart = 0;               /* flag indicating starting date found */
   int iday;                     /* day -- temporary variable for reading */
   int imonth;                   /* month -- temporary variable for reading */
   int iyear;                    /* year -- temporary variable for reading */
   int iwy;                      /* initial water year */
   int iwyjd;                    /* initial water year julian day */
   int iyr;                      /* index for year array */
   float val_miss;               /* missing value code in input file */
   float value;                  /* data value read from input file */
   void wyjdate();               /* function to determine water year julian
                                    date from calendar date */

   /* Find start and end dates of data and convert to water year format;
      also read missing data value */

   while (getln(line, fpin1) != EOF) {
      /* printf("Reading line of csv file -- first section:\n%s\n", line); */
      if (strncmp(line, "@H", 2) == 0)
         break;
      else if (strncmp(line, "date_start", 10) == 0 ||
               strncmp(line, "DATE_START", 10) == 0 ||
               strncmp(line, "Date_Start", 10) == 0 ||
               strncmp(line, "Date_start", 10) == 0) {
         sscanf(&line[11], "%d%d%d", &yr_start, &mo_start, &dy_start);
         wyjdate(yr_start, mo_start, dy_start, &iwy, &iwyjd);
         istart = 1;
         /* printf("Start wy = %d   Start wyjd = %d\n", iwy, iwyjd); */
     }
      else if (strncmp(line, "date_end", 8) == 0 ||
               strncmp(line, "DATE_END", 8) == 0 ||
               strncmp(line, "Date_End", 8) == 0 ||
               strncmp(line, "Date_end", 8) == 0) {
         sscanf(&line[9], "%d%d%d", &yr_end, &mo_end, &dy_end);
         wyjdate(yr_end, mo_end, dy_end, &fwy, &fwyjd);
         iend = 1;
         /* printf("End wy = %d   End wyjd = %d\n", fwy, fwyjd); */
      }
      else if (strncmp(line, "missing_value", 13) == 0 ||
               strncmp(line, "MISSING_VALUE", 13) == 0 ||
               strncmp(line, "Missing_Value", 13) == 0 ||
               strncmp(line, "Missing_value", 13) == 0) {
         sscanf(&line[14], "%f", &val_miss);
         imiss = 1;
      }
   }
   if (istart != 1) {
      printf("\nStarting date not found in input file ... terminating ...\n");
      exit(0);
   }
   if (iend != 1) {
      printf("\nEnding date not found in input file ... terminating ...\n");
      exit(0);
   }
   if (imiss != 1)
      printf("\nMissing data value not found in input file ... continuing ...\n");

   /* Calculate number of years and allocate some array space */

   nyear = fwy - iwy + 1;
   firstday = ivector(nyear);
   lastday = ivector(nyear);
   year = ivector(nyear);
   year[0] = iwy;

   /* Set first and last days for each year */

   for (k = 0; k < nyear; k++) {
      if (k == 0)
         firstday[k] = iwyjd;
      else
         firstday[k] = 1;
      if (k == (nyear - 1))
         lastday[k] = fwyjd;
      else
         lastday[k] = 365 + isleap(iwy + k);
      /* printf("Year = %d   firstday = %d   lastday = %d\n", year[k], firstday[k], lastday[k]); */
   }

   /* Set dayfrac to zero for now -- only daily data handled */

/* dayfrac = 0; variable deleted in Version 4.7*/

   /* Read name of data type (precip, tmax, or tmin) from line
      starting with "@H" (already read in previous while loop) */

   j = 0;
   k = 8;
   while (line[k] != '[') {
      dataname[j]= line[k];
      j++;
      k++;
   }
   dataname[j] = '\0';

   /* Read station names and count the number of stations */

   while (getln(line, fpin1) != EOF) {
      /* printf("Reading line of csv file -- second section:\n%s\n", line); */
      if (strncmp(line, "name", 4) == 0 ||
          strncmp(line, "NAME", 4) == 0 ||
          strncmp(line, "Name", 4) == 0) {
         /* printf("Parsing line\n%s\n", line); */
         i = 0;
         j = 0;
         k = 6;
         while (1) {
            while (line[k] != ',' && line[k] != '\0') {
               buf[j] = line[k];
               j++;
               k++;
            }
            buf[j] = '\0';
            strcpy(sta[i].id, buf);
            if (line[k] == '\0')
               break;
            else {
               i++;
               j = 0;
               k++;
            }
         }
         break;
      }
   }
   if (i > 0) {
      nsta = i + 1;
	  /* printf("Number of stations = %d\n", nsta); */
   }
   else {
      printf("\nStation names not found and station count not made ...\n");
      exit(0);
   }

   /* Read station elevations, eastings, northings */

   while (getln(line, fpin1) != EOF) {
      /* printf("Reading line of csv file -- third section:\n%s\n", line); */

      /* Elevations */

      if (strncmp(line, "elevation", 9) == 0 ||
          strncmp(line, "ELEVATION", 9) == 0 ||
          strncmp(line, "Elevation", 9) == 0) {
         i = 0;
         j = 0;
         k = 11;
         while (1) {
            while (line[k] != ',' && line[k] != '\0') {
               buf[j] = line[k];
               j++;
               k++;
            }
            buf[j] = '\0';
            sta[i].elev = (float) (atof(buf) / 1000);
            if (line[k] == '\0')
               break;
            else {
               i++;
               j = 0;
               k++;
            }
         }
         ielev = 1;
      }

      /* Eastings */

      else if (strncmp(line, "easting", 7) == 0 ||
               strncmp(line, "EASTING", 7) == 0 ||
               strncmp(line, "Easting", 7) == 0) {
         i = 0;
         j = 0;
         k = 9;
         while (1) {
            while (line[k] != ',' && line[k] != '\0') {
               buf[j] = line[k];
               j++;
               k++;
            }
            buf[j] = '\0';
            sta[i].east = (float) atof(buf);
            if (line[k] == '\0')
               break;
            else {
               i++;
               j = 0;
               k++;
            }
         }
         ieast = 1;
      }

      /* Northings */

      else if (strncmp(line, "northing", 8) == 0 ||
               strncmp(line, "NORTHING", 8) == 0 ||
               strncmp(line, "Northing", 8) == 0) {
         i = 0;
         j = 0;
         k = 10;
         while (1) {
            while (line[k] != ',' && line[k] != '\0') {
               buf[j] = line[k];
               j++;
               k++;
            }
            buf[j] = '\0';
            sta[i].north = (float) atof(buf);
            if (line[k] == '\0')
               break;
            else {
               i++;
               j = 0;
               k++;
            }
         }
         inorth = 1;
      }

      if (ielev == 1 && ieast == 1 && inorth == 1)
         break;
   }

   if (ielev == 0) {
      printf("\nStation elevations not found ... terminating ...\n");
      exit(0);
   }
   if (ieast == 0) {
      printf("\nStation eastings not found ... terminating ...\n");
      exit(0);
   }
   if (inorth == 0) {
      printf("\nStation northings not found ... terminating ...\n");
      exit(0);
   }

   /* Allocate space for data matrix and initialize to missing */

   for (i = 0; i < nsta; i++) {
      sta[i].data = matrix(366, nyear);
      for (j = 0; j < 366; j++)
         for (k = 0; k < nyear; k++)
            sta[i].data[j][k] = missing;
   }

   /* Read data */

   iyr = 0;
   while (getln(line, fpin1) != EOF) {
      if (line[0] == ',') {

         /* First get the date */

         sscanf(&line[1], "%d%d%d", &iyear, &imonth, &iday);
         wyjdate(iyear, imonth, iday, &iwy, &iwyjd);
         if (iwy != year[iyr])
            year[++iyr] = iwy;

         /* Then advance to the comma following the date
            and proceed parsing out the station values */

         k = 1;
         while (line[k] != ',')
            k++;

         i = 0;
         j = 0;
         k++;
         while (1) {
            while (line[k] != ',' && line[k] != '\0') {
               buf[j] = line[k];
               j++;
               k++;
            }
            buf[j] = '\0';
            value = (float) atof(buf);
            /* Set data value if not missing */
            if ((val_miss < 0.0 && value > (val_miss + 0.1)) ||
                (val_miss > 0.0 && value < (val_miss - 0.1)))
               sta[i].data[iwyjd-1][iyr] = value;
            /* Check for end of line, otherwise advance a character */
            if (line[k] == '\0')
               break;
            else {
               i++;
               j = 0;
               k++;
            }
         }
      }
   }

   fclose(fpin1);
}
