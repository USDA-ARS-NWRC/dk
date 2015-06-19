/*
 *    readgrid.c
 *
 *    David Garen   7/91, 3/94, 7/95, 1/00
 *
 *    Read grid data -- latitude, longitude, and elevation or GRASS and
 *    ARC/INFO grid formats
 *
 *    Changed strncmp() to strnicmp() function use to do better job
 *    of case-insensitive read of grid header information
 *    D. Garen   9 May 2003
 *
 *    Added zones
 *    D. Garen   8 November 2005
 *
 *    Added zoneseq array to allow zone output in numerical order
 *    D. Garen   28 November 2012
 *
 *    Added "NODATA_value" as possibility in ARC/INFO grid format read
 *    D. Garen   30 November 2012
 *
 *    Added recognition of NODATA in ARC/INFO-format elevation grid and
 *    assigning a use flag
 *    D. Garen   14 February 2013
 *
 *    Modified for Version 4.8:
 *    Use mask and zone grids to exclude grid cells that have valid
 *    elevation but lie outside mask and/or zone grid
 */

#include <stdio.h>
#include <string.h>

#include "dk_x.h"

void readgrid()
{
   double atof();                /* ascii-to-float function */
   char buf[4];                  /* buffer for parsing lat and long */
   char clat[7];                 /* latitude (character string) */
   char clng[8];                 /* longitude (character string) */
   float decmin;                 /* decimal minutes for lat and long */
   float decsec;                 /* decimal seconds for lat and long */
   int i, j, k;                  /* loop indexes and counters */
   int len;                      /* string length */
   double rnorth;                /* northing for row in grid */

   i = -1;
   ngriduse = 0;
   nmask = 0;

   if (icoord == 1 || icoord == 2) {

      /* Column grid format */

      while (getln(line, fpin2) != EOF) {
         i++;
         if (icoord == 1)
            sscanf(line, "%s%s%f", clat, clng, &grid[i].elev);
         else if (icoord == 2)
            sscanf(line, "%f%f%f", &grid[i].north, &grid[i].east,
                   &grid[i].elev);
         grid[i].elev /= 1000;
   
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
            grid[i].north = (float) (atof(buf) + decmin + decsec);
   
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
            grid[i].east = (float) (atof(buf) + decmin + decsec);
         }
      }
      ngrid = i + 1;
   }
   else if (icoord == 3) {

      /* GRASS grid format */

      for (j = 0; j < 6; j++) {
         getln(line, fpin2);
         if (strncmp(line, "north", 5) == 0 ||
             strncmp(line, "NORTH", 5) == 0)
            sscanf(line, "%*s%lf", &grass.north);
         else if (strncmp(line, "south", 5) == 0 ||
                  strncmp(line, "SOUTH", 5) == 0)
            sscanf(line, "%*s%lf", &grass.south);
         else if (strncmp(line, "east", 4) == 0 ||
                  strncmp(line, "EAST", 4) == 0)
            sscanf(line, "%*s%lf", &grass.east);
         else if (strncmp(line, "west", 4) == 0 ||
                  strncmp(line, "WEST", 4) == 0)
            sscanf(line, "%*s%lf", &grass.west);
         else if (strncmp(line, "rows", 4) == 0 ||
                  strncmp(line, "ROWS", 4) == 0)
            sscanf(line, "%*s%d", &grass.rows);
         else if (strncmp(line, "cols", 4) == 0 ||
                  strncmp(line, "COLS", 4) == 0)
            sscanf(line, "%*s%d", &grass.cols);
         if (imask == 1)
            getln(line, fpin3);
         if (izone == 1)
            getln(line, fpin4);
      }
///* Debug
fprintf(fpout, "\nreadgrid -- grass.north = %.8f\n", grass.north);
fprintf(fpout, "readgrid -- grass.south = %.8f\n", grass.south);
fprintf(fpout, "readgrid -- grass.east = %.8f\n", grass.east);
fprintf(fpout, "readgrid -- grass.west = %.8f\n", grass.west);
fprintf(fpout, "readgrid -- grass.rows = %d\n", grass.rows);
fprintf(fpout, "readgrid -- grass.cols = %d\n", grass.cols);
   //End debug */
      grass.nsres = (grass.north - grass.south) / (double) grass.rows;
      grass.ewres = (grass.east - grass.west) / (double) grass.cols;
/* Debug
fprintf(fpout, "\nreadgrid:  grass.nsres = %.8f     grass.ewres = %.8f\n",
        grass.nsres, grass.ewres);
   End debug */
      for (j = 0; j < grass.rows; j++) {
         rnorth = grass.north - ((((double) j) + 0.5) * grass.nsres);
/* Debug
fprintf(fpout, "\nrow %d -- rnorth = %f\n", j, rnorth);
   End debug */
         for (k = 0; k < grass.cols; k++) {
            i++;
            fscanf(fpin2, "%f", &grid[i].elev);
            grid[i].elev /= 1000;
            grid[i].north = (float) rnorth;
            grid[i].east = (float) (grass.west + ((((double) k) + 0.5) * grass.ewres));
            if (imask == 1) {
               fscanf(fpin3, "%d", &grid[i].mask);
               if (grid[i].mask == 1)
                  nmask++;
            }
            if (izone == 1)
               fscanf(fpin4, "%d", &grid[i].zone);
/* Debug
fprintf(fpout, "grid cell %d:  north = %.8f   east = %.8f   elev = %.3f\n",
        i, grid[i].north, grid[i].east, grid[i].elev);
   End debug */
         }
      }
      ngrid = i + 1;
   }
   else if (icoord == 4) {

      /* ARC/INFO grid format */

      for (j = 0; j < 6; j++) {
         getln(line, fpin2);
         if (strncmp(line, "ncols", 5) == 0 ||
             strncmp(line, "NCOLS", 5) == 0)
            sscanf(line, "%*s%d", &arc.cols);
         else if (strncmp(line, "nrows", 5) == 0 ||
                  strncmp(line, "NROWS", 5) == 0)
            sscanf(line, "%*s%d", &arc.rows);
         else if (strncmp(line, "xllcorner", 9) == 0 ||
                  strncmp(line, "XLLCORNER", 9) == 0)
            sscanf(line, "%*s%lf", &arc.xll);
         else if (strncmp(line, "yllcorner", 9) == 0 ||
                  strncmp(line, "YLLCORNER", 9) == 0)
            sscanf(line, "%*s%lf", &arc.yll);
         else if (strncmp(line, "cellsize", 8) == 0 ||
                  strncmp(line, "CELLSIZE", 8) == 0)
            sscanf(line, "%*s%f", &arc.cell);
         else if (strncmp(line, "nodata_value", 12) == 0 ||
                  strncmp(line, "NODATA_VALUE", 12) == 0 ||
                  strncmp(line, "NODATA_value", 12) == 0) {
            sscanf(line, "%*s%f", &arc.nodata);
            arc.nodata += 0.1;
         }
         if (imask == 1)
            getln(line, fpin3);
         if (izone == 1)
            getln(line, fpin4);
      }
      for (j = 0; j < arc.rows; j++) {
         rnorth = arc.yll + ((double) arc.rows - (double) j - 0.5) * arc.cell;
         for (k = 0; k < arc.cols; k++) {
            i++;
            fscanf(fpin2, "%f", &grid[i].elev);
            if (grid[i].elev > arc.nodata) {
               ngriduse++;
               grid[i].use = 1;
               grid[i].elev /= 1000;
            }
            else
               grid[i].use = 0;
            grid[i].north = (float) rnorth;
            grid[i].east = (float) (arc.xll + ((double) k + 0.5) * arc.cell);
            if (imask == 1) {
               fscanf(fpin3, "%d", &grid[i].mask);
               if (grid[i].mask == 1)
                  nmask++;
            }
            if (izone == 1)
               fscanf(fpin4, "%d", &grid[i].zone);
         }
      }
      ngrid = i + 1;
   }
   fclose(fpin2);
   if (imask == 1)
      fclose(fpin3);
   if (izone == 1)
      fclose(fpin4);

   /* Finished with reading grids */

   /* Use mask and zone grids to exclude grid cells that have a valid
      elevation but lie outside the mask and/or zone grid.  This assumes
      that grid cells outside one or both of these grids are not of
      interest, therefore processing can be speeded up by not including
      these cells in the calculations. The grid use flag is turned off
      for these cells. */

   if ((icoord == 3 || icoord == 4) && (imask == 1 || izone == 1)) {
      for (i = 0; i < ngrid; i++) {
         if (grid[i].use == 1) {
            /* Only mask grid specified */
            if (imask == 1 && izone != 1) {
               if (grid[i].mask != 1) {
                  grid[i].use = 0;
                  ngriduse--;
               }
            }
            /* Only zone grid specified */
            else if (imask != 1 && izone == 1) {
               if (grid[i].zone <= 0) {
                  grid[i].use = 0;
                  ngriduse--;
               }
            }
            /* Both mask and zone grid specified --
               exclude cell only if outside both */
            else {
               if (grid[i].mask != 1 && grid[i].zone <= 0) {
                  grid[i].use = 0;
                  ngriduse--;
               }
            }
         }
      }
   }

   /* Some final zone processing: */

   /* If zones are used, determine the zone numbers
      and the cell counts for each zone */

   if (izone == 1 && (icoord == 3 || icoord == 4)) {
      nzone = 0;
      for (i = 0; i < ngrid; i++) {
         if (grid[i].zone > 0) {
            for (j = 0; j < nzone; j++) {
               if (grid[i].zone == zone[j].number) {
                  zone[j].ncells++;
                  break;
               }
            }
            if (j >= nzone) {
               zone[nzone].number = grid[i].zone;
               zone[nzone].ncells++;
               nzone++;
            }
         }
      }

      /* Set values in zoneseq array to allow zone output in
         numerical order */

      for (j = 0; j < nzone; j++)
         zoneseq[zone[j].number-1] = j;

/* Debug
fprintf(fpzone, "Zone summary:\n\nNumber of zones = %d\n\nZone      Cells\n",
        nzone);
for (j = 0; j < nzone; j++)
   fprintf(fpzone, "%4d   %8d\n", zone[j].number, zone[j].ncells);
fprintf(fpzone, "\n");
fprintf(fpzone, "Year WYD");
for (j = 0; j < nzone; j++)
   fprintf(fpzone, "%7d", zone[j].number);
fprintf(fpzone, "\n\n");
   End debug */

   }
}
