/*
 *    wyjdate.c
 *
 *    David Garen  10/2012
 *
 *    Determine water year julian date from calendar date
 *
 */

#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void wyjdate(iy, im, id, iwy, iwyjd)
int iy;                          /* calendar year */
int im;                          /* month */
int id;                          /* calendar day */
int *iwy;                        /* water year (return value) */
int *iwyjd;                      /* water year julian day (return value) --
                                    sequential day number beginning 1 October */
{
   int leapday;                  /* 1 = leap year, 0 = not leap year */

   /* Determine if the year is a leap year */

   leapday = isleap(iy);

   /* Set water year and water year julian day based on calendar date */

   if (im == 10) {
      /* October */
      *iwy = iy + 1;
      *iwyjd = id;
   }
   else if (im == 11) {
      /* November */
      *iwy = iy + 1;
      *iwyjd = id + 31;
   }
   else if (im == 12) {
      /* December */
      *iwy = iy + 1;
      *iwyjd = id + 61;
   }
   else if (im == 1) {
      /* January */
      *iwy = iy;
      *iwyjd = id + 92;
   }
   else if (im == 2) {
      /* February */
      *iwy = iy;
      *iwyjd = id + 123;
   }
   else if (im == 3) {
      /* March */
      *iwy = iy;
      *iwyjd = id + 151 + leapday;
   }
   else if (im == 4) {
      /* April */
      *iwy = iy;
      *iwyjd = id + 182 + leapday;
   }
   else if (im == 5) {
      /* May */
      *iwy = iy;
      *iwyjd = id + 212 + leapday;
   }
   else if (im == 6) {
      /* June */
      *iwy = iy;
      *iwyjd = id + 243 + leapday;
   }
   else if (im == 7) {
      /* July */
      *iwy = iy;
      *iwyjd = id + 273 + leapday;
   }
   else if (im == 8) {
      /* August */
      *iwy = iy;
      *iwyjd = id + 304 + leapday;
   }
   else if (im == 9) {
      /* September */
      *iwy = iy;
      *iwyjd = id + 335 + leapday;
   }
}