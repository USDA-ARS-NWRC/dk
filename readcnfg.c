/*
 *    readcnfg.c
 *
 *    Chance Lerro   11/2007
 *    John Langland  04/2008;  cleaned up, clarified, simplified
 *
 *    Read input data from a file for batch run settings/config
 *    This is an alternative to the interactive input.
 *
 *    This was originally implemented in support of a Java user
 *    interface.  This is no longer supported, but the configuration
 *    file method of specifying program run information has proven
 *    to be useful on its own.
 *
 *    D. Garen, 17 October 2012:
 *    Added parameter "input-format-csv" to configuration file,
 *    which corresponds to variable "iomscsv" in program.  This
 *    enables input to be in OMS-compatible csv format.
 *
 *    Modified for Version 4.7:
 *    Added time step, changed keywords from referring to "day"
 *    to referring to "period"
 *
 *    Modified for Version 4.8:
 *    Removed replacement of double backslashes in file path
 *    names.  This was originally put in to accommodate use
 *    of a Java user interface.  This interface has not been
 *    maintained and is no longer supported.  Source file
 *    replace.c has also been removed from the system.
 *    
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "dk_x.h"
#include "dk_m.h"

extern char elevfile[];
extern char maskfile[];
extern char zonefile[];
extern char infile[];
extern char outputdir[];
extern char kw_filename[];
extern int ioutputdir;
extern int iprintdistances;
extern int iprintelevation;
extern int iprintinput;
extern int iprintresiduals;
extern int iprintweights;
extern int iomscsv;
extern int ikwfile;
extern int mtper;

/*
 *  Read through the configuation file and set parameters accordingly.
 *  Parameters are set the same way as the with the interactive input.
 *
 *  Lines in the configuration file are name-value pairs.
 *  They are separated by an equal sign.
 *  Lines beginning with a '#' are comments
 *
 */

void get_file_configuration(char configuration_filename[101])
{
   FILE *fp_configuration_file;   // settings file
   char *value;                   // value of name-value pair
   char *name;                    // name of name-value pair
   char buffer[200];              // line buffer for a line read from input
   int line_length = 0;           // length of C-string in buffer

   int buffer_size = sizeof(buffer)/sizeof(char);

   /* Open the configuration file */

   if ((fp_configuration_file = fopen(configuration_filename, "r")) == NULL) {
      printf("\n\nError opening file %s\nProgram terminated ...\n",
         configuration_filename);
      exit(0);
   }

   /*  Loop through every line in the configuration file.  
    *  Parse each line into a name-value pair.
    *  Check if the name matches one of the expected names,
    *  and if it does, take action according to the value.
    *  NOTE: if a line in the configuration file has a name but no value,
    *  it is skipped entirely -- the program goes on to the next line.
    */

   while (fgets(buffer, buffer_size, fp_configuration_file) != NULL)
   {
      if ((line_length = strlen(buffer)) == 0)
         continue;

      if (buffer[0] == '#')
         continue;

      buffer[line_length - 1] = '\0';  // overwrite newline with null character

      name = strtok (buffer, "=");

      if (name ==  NULL || strlen(name) < 3) // e.g. x=1
         continue;

      value = strtok (NULL, "=");

      if (value == NULL)
         continue;
          
      if (strcmp(name, "input-data-file-name") == 0) {
         if (strlen(value) == 0) {
            printf("\n\nNo station-file-name\nProgram terminated ...\n");
            exit(0);
         }
         if ((fpin1 = fopen(value, "r")) == NULL) {
            printf("\n\nError opening file %s\nProgram terminated ...\n", value);
            exit(0);
         }
         strcpy(infile, value);
      }
      else if (strcmp(name, "type-of-data") == 0) {
         switch (value[0]) {
            case '1':
               // Precipitation
               type = 1;
               missing = 9999.8f;
               accum = 8888.7f;
               replace = 8999.9f;
               break;
            case '2':
               // Temperature
               type = 2;
               missing = accum = 9999.8f;
               break;
            case '3':
               // Snow-water equivalent
               type = 3;
               missing = accum = 9999.8f;
               break;
            case '4':
               // Other
               type = 4;
               missing = accum = 9999.8f;
               break;
            default:
               printf("\n\nError, input-data-type = %c not allowed\n"
                  "Program terminated ...\n", value[0]);
               exit(0);
         }
      }
      else if (strcmp(name, "time-step") == 0) {
         switch (value[0]) {
            case '1':
               // Hourly
               mtper = 8784;
               break;
            case '2':
               // Daily
               mtper = 366;
               break;
            case '3':
               // Monthly
               mtper = 12;
               break;
            case '4':
               // Yearly
               mtper = 1;
               break;
            default:
               printf("\n\nError, time-step = %c not allowed\n"
                  "Program terminated ...\n", value[0]);
               exit(0);
         }
      }
      else if (strcmp(name, "coord-system") == 0) {
         switch (value[0]) {
            case '1':
               // Latitude and longitude, column grid format
               icoord = 1;
               break;
            case '2':
               // Easting and northing, column grid format
               icoord = 2;
               break;
            case '3':
               // Easting and northing, GRASS grid format
               icoord = 3;
               break;
            case '4':
               // Easting and northing, ARC/INFO grid format
               icoord = 4;
               break;
            default:
               printf("\n\nError, coord-system = %c not allowed\n"
                  "Program terminated ...\n", value[0]);
               exit(0);
         }
      }
      else if (strcmp(name, "elevation-grid-file-name") == 0) {
         if (strlen(value) == 0) {
            printf("\n\nNo elevation-grid-file-name\nProgram terminated ...\n");
            exit(0);
         }
         if ((fpin2 = fopen(value, "r")) == NULL) {
            printf("\n\nError opening file %s\nProgram terminated ...\n", value);
            exit(0);
         }             
         strcpy(elevfile, value);
      }
      else if (strcmp(name, "watershed-mask-file-name") == 0) {
         if (strlen(value) == 0) {
            imask = 0;
         }
         else {
            if ((fpin3 = fopen(value, "r")) == NULL) {
               printf("\n\nError opening file %s\nProgram terminated ...\n", value);
               exit(0);
            }
            strcpy(maskfile, value);
            imask = 1;
            nmask = 0;
         }
      }
      else if (strcmp(name, "zone-grid-file-name") == 0) {
         if (strlen(value) == 0) {
            izone = 0;
         }
         else {
            if ((fpin4 = fopen(value, "r")) == NULL) {
               printf("\n\nError opening file %s\nProgram terminated ...\n", value);
               exit(0);
            }
            strcpy(zonefile, value);
            izone = 1;
         }
      }
      else if (strcmp(name, "output-format") == 0) {
         switch (value[0]) {
            case '1':
               // Mean areal values -- tabular (one column per year)
               iout = 1;
               break;
            case '2':
               // Timestep grids -- GRASS format, plus (1) above
               iout = 2;
               break;
            case '3':
               // Timestep grids -- ARC/INFO format, plus (1) above
               iout = 3;
               break;
            case '4':
               // Timestep grids -- IPW format, plus (1) above
               iout = 4;
               break;
            default:
               printf("\n\nError, output-format = %c not allowed\n"
                  "Program terminated ...\n", value[0]);
               exit(0);
         }
      }
      else if (strcmp(name, "beginning-period-number") == 0) {
         int result = sscanf(value, "%d", &igridout1);
         if (result != 1) {
            printf("\n\nError, no beginning-period-number\nProgram terminated ...\n");
            exit(0);
         }
         igridout1--;
      }
      else if (strcmp(name, "ending-period-number") == 0) {
         int result = sscanf(value, "%d", &igridout2);
         if (result != 1) {
            printf("\n\nError, no ending-period-number\nProgram terminated ...\n");
            exit(0);
         }
         igridout2--;
      }
      else if (strcmp(name, "output-precision") == 0) {
         switch (value[0]) {
            case '2':
               // Integer
               igridpr = 2;
               break;
            case '3':
               // Integer multiplied by 10
               igridpr = 3;
               break;
            case '1':
            default:
                // Floating point, one decimal place (default)
                igridpr = 1;
         }
      }
      else if (strcmp(name, "output-file-name") == 0) {
         if ((fpout = fopen(value, "w")) == NULL) {
            printf("\n\nError opening file %s\nProgram terminated ...\n", name);
            exit(0);
         }
      }
      else if (strcmp(name, "zone-output-file-name") == 0) {
         if ((fpzone = fopen(value, "w")) == NULL) {
            printf("\n\nError opening file %s\nProgram terminated ...\n", value);
            exit(0);
         }
      }
      else if (strcmp(name, "regression-method") == 0) {
         switch (value[0]) {
            case '2':
               // least absolute deviations
               irmeth = 2;
               break;
            default:
               // least squares
               irmeth = 1;
         }
      }
      else if (strcmp(name, "station-weighting-method") == 0) {
         switch (value[0]) {
            case '2':
               // Equal
               iwt = 2;
               break;
            default:
               // Distance
               iwt = 1;
         }
      }
      else if (strcmp(name, "timesteps-per-period") == 0) {
         dpp = atoi(value);
         nper = mtper / dpp;
      }
      else if (strcmp(name, "grid-output-file-dir") == 0) {
          strcpy(outputdir, value);
          ioutputdir = 1;
      }
      else if (strcmp(name, "print-distances") == 0) {
         if (strcmp(value, "true") == 0)
            iprintdistances = 1;
      }
      else if (strcmp(name, "print-regressions") == 0) {
         if (strcmp(value, "true") == 0)
            iprintelevation = 1;
      }
      else if (strcmp(name, "print-input") == 0) {
         if (strcmp(value, "true") == 0)
            iprintinput = 1;
      }
      else if (strcmp(name, "print-residuals") == 0) {
         if (strcmp(value, "true") == 0)
            iprintresiduals = 1;
      }
      else if (strcmp(name, "print-weights") == 0) {
         if (strcmp(value, "true") == 0)
            iprintweights = 1;
      }
      else if (strcmp(name, "input-format-csv") == 0) {
         if (strcmp(value, "true") == 0)
            iomscsv = 1;
      }
      else if (strcmp(name, "kriging-weights-file-name") == 0) {
         if (strlen(value) == 0) {
            ikwfile = 0;
         }
         else {
            strcpy(kw_filename, value);
            ikwfile = 1;
         }
      }
   }

   fclose(fp_configuration_file);
}
