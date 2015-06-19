/*
 *    caldate.c
 *
 *    David Garen  7/2011
 *
 *    Determine calendar month and day from water year julian day for a given year
 *
 */

#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void caldate(iy, id, month, day)
int iy;                          /* water year */
int id;                          /* water year julian day (sequential number beginning Oct 1) */
int *day;                        /* day of month (return value) */
int *month;                      /* calendar year month number (return value) */
{
   int leapday;                  /* 1 = leap year, 0 = not leap year */

   /* Determine if the year is a leap year */

   leapday = isleap(iy);

   /* Set month and day number based on year and water year julian day */

   if (id >= 1 && id <= 31) {
      /* October */
      *month = 10;
      *day = id;
   }
   else if (id >= 32 && id <= 61) {
      /* November */
      *month = 11;
      *day = id - 31;
   }
   else if (id >= 62 && id <= 92) {
      /* December */
      *month = 12;
      *day = id - 61;
   }
   else if (id >= 93 && id <= 123) {
      /* January */
      *month = 1;
      *day = id - 92;
   }
   else if (id >= 124 && id <= (151 + leapday)) {
      /* February */
      *month = 2;
      *day = id - 123;
   }
   else if (id >= (152 + leapday) && id <= (182 + leapday)) {
      /* March */
      *month = 3;
      *day = id - 151 - leapday;
   }
   else if (id >= (183 + leapday) && id <= (212 + leapday)) {
      /* April */
      *month = 4;
      *day = id - 182 - leapday;
   }
   else if (id >= (213 + leapday) && id <= (243 + leapday)) {
      /* May */
      *month = 5;
      *day = id - 212 - leapday;
   }
   else if (id >= (244 + leapday) && id <= (273 + leapday)) {
      /* June */
      *month = 6;
      *day = id - 243 - leapday;
   }
   else if (id >= (274 + leapday) && id <= (304 + leapday)) {
      /* July */
      *month = 7;
      *day = id - 273 - leapday;
   }
   else if (id >= (305 + leapday) && id <= (335 + leapday)) {
      /* August */
      *month = 8;
      *day = id - 304 - leapday;
   }
   else if (id >= (336 + leapday) && id <= (365 + leapday)) {
      /* September */
      *month = 9;
      *day = id - 335 - leapday;
   }
}