/*
 *    ipwout.c
 *
 *    David Garen  11/1997, 1/2000, 5/2006
 *
 *    Write daily grids to output files in IPW format
 *
 *    All images use 8-bit quantization.  The maximum and minimum values used
 *    are:
 *       Precipitation:          MAX = 127.5 mm    MIN = 0 mm
 *       Temperature:            MAX = 63.5 C      MIN = -64.0 C
 *       Snow water equivalent:  MAX = 2550 mm     MIN = 0 mm
 *    The quantization interval (precision) is calculated as
 *    (max - min) / (2^nbits - 1).  Using 8 bits, the denominator is 255.
 *    With these values, the resolution of the data in the images is 0.5 mm for
 *    precipitation, 0.5 degrees C for temperature, and 10 mm for swe.
 *
 *    File names are as follows:  parameter (prc, tmp, swe), year, day . IPW
 *       Examples:  prc89058.ipw, tmp91255.ipw, swe93121.ipw
 *
 *    Modification 24 January 2000:
 *       Changed precipitation resolution to 0.5 mm from 1.0 mm.
 *
 *    Modification 28 March 2005:
 *       Changed file naming convention.  Example:  prc_2004_004.ipw
 *
 *    Modification 15 May 2006:
 *       Changed quantization from previous fixed values for each variable
 *       (see table above) to varying value depending on the maximum and
 *       minimum of the individual image.  In general, this will increase
 *       the precision by tailoring the quanitzation interval to each
 *       image rather than using generic values that had to cover the
 *       entire possible range for each variable type.  This still uses
 *       8 bit quantization, which allows for 256 possible values.
 *
 *    Modification 6 November 2006:
 *       Added day fraction to file naming convention.
 *       Example:  prc_2004_004_125.ipw (125 corresponds to hour 3)
 *
 *    Modification 22 November 2006:
 *       Fixed bug in calculation of delta -- fraction was upside down
 *
 *    Modification for Version 4.7:
 *       Changed file naming convention to reflect generic time periods
 *       instead of days, and removed day fraction.
 *       Example:  prc_2004_6358.ipw
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void ipwout(iy, ip)
int iy;                          /* year */
int ip;                          /* period (sequential number beginning Oct 1) */
{
   char buf[6];                  /* buffer for file name building */
   float delta;                  /* number of image units per data unit
                                    (reciprocal of precision) */
   FILE *fpipw;                  /* output file pointer */
   int i;                        /* loop index */
/* Debug
   int j, k;
   End debug */
   int ival;                     /* int quantized data value */
   float max;                    /* maximum value of variable */
   float min;                    /* minimum value of variable */
   char outfile[21];             /* output file name */

   /* Find minimum and maximum values in image and compute delta */

   min = gprec[0];
   max = gprec[0];
   for (i = 1; i < ngrid; i++) {
      if (gprec[i] < min)
         min = gprec[i];
      if (gprec[i] > max)
         max = gprec[i];
   }
   delta = (float) (255.0 / (max - min));

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
   strcat(outfile, ".ipw");
   if ((fpipw = fopen(outfile, "w")) == NULL) {
      printf("\n\nError opening file %s.\n", outfile);
      return;
   }

   /* Write header information: */

   /* Basic image header */

   if (icoord == 3)
      fprintf(fpipw, "%s\n%s\n%s%d\n%s%d\n%s\n%s",
                     "!<header> basic_image_i -1", "byteorder = 3210",
                     "nlines = ", grass.rows, "nsamps = ", grass.cols,
                     "nbands = 1", "annot = ");
   else if (icoord == 4)
      fprintf(fpipw, "%s\n%s\n%s%d\n%s%d\n%s\n%s",
                     "!<header> basic_image_i -1", "byteorder = 3210",
                     "nlines = ", arc.rows, "nsamps = ", arc.cols,
                     "nbands = 1", "annot = ");
   if (type == 1)
      fprintf(fpipw, "precipitation");
   else if (type == 2)
      fprintf(fpipw, "temperature");
   else if (type == 3)
      fprintf(fpipw, "snow water equivalent");
   else
      fprintf(fpipw, "general data type");
   fprintf(fpipw, " for year %04d, period %04d\n", iy, ip+1);
/* fprintf(fpipw, " for year %04d, day %03d.%03d\n", iy, id+1, dayfrac); */

   fprintf(fpipw, "%s\n%s\n%s\n%s\n",
                  "!<header> basic_image 0", "bytes = 1", "bits = 8",
                  "history = created by DK program");

   /* Geo header */

   if (icoord == 3)
      fprintf(fpipw, "%s\n%s%f\n%s%f\n%s%f\n%s%f\n%s\n%s\n",
                     "!<header> geo 0", "bline = ",
                     (grass.north - 0.5 * grass.nsres),
                     "bsamp = ", (grass.west + 0.5 * grass.ewres), "dline = ",
                     -grass.nsres, "dsamp = ", grass.ewres, "units = meters",
                     "coord_sys_ID = UTM");
   else if (icoord == 4)
      fprintf(fpipw, "%s\n%s%f\n%s%f\n%s%f\n%s%f\n%s\n%s\n",
                     "!<header> geo 0", "bline = ",
                     (arc.yll + ((double) arc.rows - 0.5) * arc.cell),
                     "bsamp = ", (arc.xll + 0.5 * arc.cell), "dline = ",
                     -arc.cell, "dsamp = ", arc.cell, "units = meters",
                     "coord_sys_ID = UTM");

   /* Linear quantization header */

   fprintf(fpipw, "!<header> lq 0\nunits = ");
   if (type == 1 || type == 3)
      fprintf(fpipw, "mm\nmap = 0 %6.1f\nmap = 255 %6.1f\n", min, max);
   else if (type == 2)
      fprintf(fpipw, "C\nmap = 0 %5.1f\nmap = 255 %5.1f\n", min, max);
   else
      fprintf(fpipw, "unknown\nmap = 0 %5.1f\nmap = 255 %5.1f\n", min, max);

   /* Image header and CNTL-L to separate header from data */

   fprintf(fpipw, "!<header> image -1\n");

   /* Write grid values */

   for (i = 0; i < ngrid; i++) {
      ival = (int) ((gprec[i] - min) * delta + 0.5);
      fprintf(fpipw, "%c", (char) (ival & 0xFF));
   }

/* Debug
   k = -1;
   for (i = 0; i < arc.rows; i++) {
      for (j = 0; j < arc.cols; j++) {
         ival = (int) ((gprec[++k] - min) * delta + 0.5);
         fprintf(fpipw, " %d", ival);
      }
      fprintf(fpipw, "\n");
   }
   End debug */

   fclose(fpipw);
}
