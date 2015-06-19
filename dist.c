/*
 *    dist_ll.c
 *
 *    David Garen  11/89
 *
 *    Compute distances (kilometers) between stations
 *    based on latitude and longitude
 */

#include <math.h>
#include <stdio.h>

#define ABS(a) ((a) >= 0 ? (a) : -(a))      /* absolute value operator */

/* Debug
extern FILE *fpout;
   End debug */

double sqrt();                   /* square root function */

float dist_ll(lat1, long1, lat2, long2, ewdst, nsdst)
float lat1, long1, lat2, long2;  /* decimal latitude and longitude for
                                    two stations */
float *ewdst;                    /* east-west distance */
float *nsdst;                    /* north-south distance */
{
   int i, ip1;                   /* array indexes */
   float len1, len2;             /* lengths */
   float lenavg;                 /* average length */
   struct lldata {
      float deg;                 /* degrees */
      float len;                 /* length */
   };
   struct lldata latd[28] = {
      {25.5f, 8.833f}, {26.5f, 8.842f}, {27.5f, 8.852f}, {28.5f, 8.862f},
      {29.5f, 8.873f}, {30.5f, 8.883f}, {31.5f, 8.894f}, {32.5f, 8.905f},
      {33.5f, 8.916f}, {34.5f, 8.928f}, {35.5f, 8.939f}, {36.5f, 8.951f},
      {37.5f, 8.962f}, {38.5f, 8.974f}, {39.5f, 8.986f}, {40.5f, 8.998f},
      {41.5f, 9.011f}, {42.5f, 9.023f}, {43.5f, 9.035f}, {44.5f, 9.047f},
      {45.5f, 9.060f}, {46.5f, 9.072f}, {47.5f, 9.084f}, {48.5f, 9.096f},
      {49.5f, 9.108f}, {50.5f, 9.121f}, {51.5f, 9.133f}, {52.5f, 9.145f},
   };                            /* latitude and length of one degree of
                                    latitude (statute miles) */
   struct lldata longd[28] = {
      {25.0f, 2.729f}, {26.0f, 2.212f}, {27.0f, 1.676f}, {28.0f, 1.122f},
      {29.0f, 0.548f}, {30.0f, 9.956f}, {31.0f, 9.345f}, {32.0f, 8.716f},
      {33.0f, 8.071f}, {34.0f, 7.407f}, {35.0f, 6.725f}, {36.0f, 6.027f},
      {37.0f, 5.311f}, {38.0f, 4.579f}, {39.0f, 3.829f}, {40.0f, 3.063f},
      {41.0f, 2.281f}, {42.0f, 1.483f}, {43.0f, 0.669f}, {44.0f, 9.840f},
      {45.0f, 8.995f}, {46.0f, 8.136f}, {47.0f, 7.261f}, {48.0f, 6.372f},
      {49.0f, 5.469f}, {50.0f, 4.552f}, {51.0f, 3.621f}, {52.0f, 2.676f},
   };                            /* latitude and length of one degree of
                                    longitude (statute miles ) */
/* Debug
fprintf(fpout, "\n\ndist:  lat1=%5.2f  long1=%6.2f  lat2=%5.2f  long2=%6.2f",
        lat1, long1, lat2, long2);
   End debug */

   /* Compute north-south distance between stations (km) */

   for (i = 0; i < 27; i++) {
      ip1 = i + 1;
      if (lat1 >= latd[i].deg && lat1 <= latd[ip1].deg)
         len1 = (float) ((lat1 - latd[i].deg) / (latd[ip1].deg - latd[i].deg))
                * (latd[ip1].len - latd[i].len) + latd[i].len;
      if (lat2 >= latd[i].deg && lat2 <= latd[ip1].deg)
         len2 = (float) ((lat2 - latd[i].deg) / (latd[ip1].deg - latd[i].deg))
                * (latd[ip1].len - latd[i].len) + latd[i].len;
   }
   lenavg = (float) (1.609 * ((len1 + len2) / 2));
   *nsdst = (float) (ABS(lat1 - lat2) * lenavg);
/* Debug
fprintf(fpout, "\n       nsdst=%5.2f  len1=%5.2f  len2=%5.2f  lenavg=%5.2f",
        *nsdst, len1, len2, lenavg);
   End debug */

   /* Compute east-west distance between stations */

   for (i = 0; i < 27; i++) {
      ip1 = i + 1;
      if (lat1 >= longd[i].deg && lat1 <= longd[ip1].deg)
         len1 = (float) (((lat1 - longd[i].deg) / (longd[ip1].deg - longd[i].deg))
                * (longd[ip1].len - longd[i].len) + longd[i].len);
      if (lat2 >= longd[i].deg && lat2 <= longd[ip1].deg)
         len2 = (float) (((lat2 - longd[i].deg) / (longd[ip1].deg - longd[i].deg))
                * (longd[ip1].len - longd[i].len) + longd[i].len);
   }
   lenavg = (float) (1.609 * ((len1 + len2) / 2));
   *ewdst = (float) (ABS(long1 - long2) * lenavg);
/* Debug
fprintf(fpout, "\n       ewdst=%5.2f  len1=%5.2f  len2=%5.2f  lenavg=%5.2f",
        *ewdst, len1, len2, lenavg);
   End debug */

   /* Compute straight-line distance between stations using
      Pythagorean theorem */

   return((float) sqrt((double) (*nsdst * *nsdst + *ewdst * *ewdst)));
}

/*
 *    dist_en.c
 *
 *    David Garen   6/92
 *
 *    Compute distances between stations based on easting and northing.
 *    Eastings and northings are in meters, computed distances are in
 *    kilometers.
 */

float dist_en(north1, east1, north2, east2)
float north1, east1, north2, east2;    /* northings and eastings for two
                                          stations */
{
   float ewdst;                  /* east-west distance */
   float nsdst;                  /* north-south distance */

   ewdst = (float) (ABS(east1 - east2));
   nsdst = (float) (ABS(north1 - north2));
   return (float) ((sqrt((double) (nsdst * nsdst + ewdst * ewdst)) / 1000));
}

