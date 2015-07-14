/*
 *    dk.c
 *
 *    Detrended kriging program for computing spatially distributed
 *    precipitation, temperature, and snow water equivalent fields
 *
 *    David Garen   6/92, 8/92, 9/92, 10/93, 8/95, 5/97, 11/97
 *
 *    This program is a modified version of a previous program called
 *    map_d (Mean Areal Precipitation -- Daily data).  This modification
 *    involves a name change to be more descriptive of its current usage
 *    and capabilities, some enhancements to be more convenient to use,
 *    and a removal of some experimental and unnecessary code in the
 *    older program.  Improvements in coding and documentation have also
 *    been made to improve the program's usability and maintainability.
 *
 *    For a few years, I gave the program the name SPAM (for SPAtially
 *    distributed hydroMeteorological variables).  This was an attempt to
 *    be humorous, and the name had some personal significance as a
 *    nickname for my older son.  Since this word has taken on a negative
 *    connotation in relation to e-mail, I decided to rename the program.
 *    It is now called DK, for Detrended Kriging, which simply indicates
 *    the mathematical algorithm used.
 *
 *    Spatial fields on a grid cell basis for daily values of precipitation,
 *    temperature, and snow water equivalent are calculated by interpolating
 *    point measurements at hydrometeorological stations.  Spatial
 *    interpolation is done by detrended kriging.  Spatial fields can be
 *    output in the ASCII column-and-row format used by the GRASS GIS, and
 *    an arithmetic average of the values for all grid cells in the region
 *    is also given.
 *
 *    Elevation is accounted for by detrending, whereby a linear regression of
 *    elevation vs. precipitation, temperature, or snow water equivalent is
 *    fit to each individual period's data, and the residuals from this
 *    regression are used in the kriging calculations.                       
 *
 *    Kriging weights are calculated from the distances among stations and
 *    distances between stations and grid cells.  This is equivalent to
 *    using a linear semivariogram.  Alternatively, equal weighting can be
 *    used, which is equivalent to a flat semivariogram.
 *
 *    Grid cells are represented by latitude or northing, longitude or
 *    easting, and elevation.  These are most easily defined by a digital
 *    elevation model (DEM).
 *
 *    Input data are daily precipitation, temperature, or snow water equivalent
 *    values, and the elevation, latitude or northing, and longitude or easting
 *    of each station.
 *
 *    For calculating the relationships of the data with elevation,
 *    the daily data can either be aggregated into periods of a specified
 *    number of days in length or (for precipitation only) into storms
 *    defined as consecutive sequences of days where one or more stations
 *    has precipitation.
 *
 *    Command line switches:
 *       -c or /c    writes input data to output file in column format then
 *                   quits
 *       -d or /d    gives printout of distances among stations
 *       -e or /e    gives printout of prec, temp, or swe vs. elevation
 *                   regressions
 *       -f or /f    reads kriging weights from a file instead of
 *                   calculating them (file name follows the -f)
 *       -i or /i    gives printout of input data
 *       -o or /o    reads data in OMS-compatible csv format
 *       -r or /r    gives printout of detrended residuals
 *       -w or /w    gives printout of kriging weights
 *       -k or /k    reads all program info from a configuration file
 *                   (the name of which follows the -k)
 *
 *    Version 2.0, 8 May 1997:
 *       First version called SDHV.  Version number of 2.0 was used because
 *       this is not a new program, being derived from MAP_D, which would
 *       be considered to be version 1.x.
 *
 *    25 September 1997:  Name of program changed to SPAM.
 *
 *    Version 2.1, 28 October 1997:
 *       Small bug fix to modules swe1.c and swe2.c
 *       (added -1 to index jj calculation)
 *       Modified handling of snow line and trend between highest zero
 *       and lowest nonzero swe
 *
 *    Version 2.2, 3 November 1997:
 *       Added IPW output option.
 *
 *    Version 2.3, 7 November 1997:
 *       Added watershed mask option.
 *
 *    Version 2.4, 16 December 1997:
 *       Changed variable type for grid[i].mask to long instead of short int
 *       so that huge model for PC could be used (each element of array has
 *       to have a size in bytes equal to a power of 2, 16 bytes in this case).
 *
 *    Version 2.5, 23 December 1997:
 *       Added option for different precision of values in GRASS-format
 *       raster files.
 *
 *    Version 2.51, 9 March 1999:
 *       Fixed small bug in line 465 -- from nper=365/dpp to
 *       nper=366/dpp.  This is used to dimension some arrays, but crashed
 *       in leap years.  Problem reported by Joachim Geyer.
 *
 *    Version 2.6, 24 January 2000:
 *       Added ARC/INFO grid format and changed IPW precipitation resolution
 *       to 0.5 mm from 1.0 mm.
 *
 *    Version 2.7, 26 May 2000:
 *       Changed method of solving for kriging weights from Gauss-Jordan to
 *       LU decomposition (much faster)
 *
 *    Version 2.8, 14 April 2003
 *       Name of program changed to DK
 *
 *    Version 2.81, 9 May 2003
 *       Function readgrid.c modified to use strnicmp() instead of strncmp()
 *       to do a better job of non-case-sensitive reads of grid headers
 *
 *    Version 2.82, 28 March 2005
 *       Changed output grid file naming convention for IPW and ARC formats.
 *
 *    Version 2.9, 8 November 2005
 *       Added computation and output of zone means
 *
 *    Version 2.91, 15 May 2006
 *       Changed precision for IPW output by computing quantization interval
 *       based on maximum and minimum in each image rather than using generic
 *       values
 *
 *    Version 2.92, 6 November 2006
 *       1) Added day fraction field to input file and added day fraction to
 *       grid output file names.
 *       2) Added run information to top of main output file (input file names,
 *       day fraction).
 *       3) Changed strnicmp() back to strncmp() in readgrid.c because
 *       strnicmp() function for case-insensitive reads exists only in PC
 *       Visual C++ compiler and not in UNIX gnu C compiler.  Thus un-does
 *       the modification made for Version 2.81.  The read now assumes that
 *       all the letters are lower case or all are upper case; it will fail
 *       for mixed case keywords.
 *
 *    Version 3.0, 8 November 2006
 *       Created "other" data category to allow detrending without screening
 *       or restrictions on the slope
 *
 *    Version 3.01, 22 November 2006
 *       Fixed bug in calculation of delta in ipwout.c
 *
 *    Version 3.1, 29 May 2007
 *       Modified zoneout.c and period2.c to write a line of zone output
 *       for every day, even if all values are zero
 *
 *    Version 4.0, 20 April 2009
 *       Modified to read parameters from a "configuration" file
 *       (in addition to the interactive, scroll-up prompting).
 *       This enables, among other things, for the use of a GUI input
 *       form that writes the user's choices to the configuration file.
 *       The program can then be invoked with the -k option to read a given
 *       configuration file.  Reading the configuration file is typically
 *       created by a GUI front end (that can also read the file).
 *       This front end can also be used to invoke this program
 *
 *    Version 4.1, 1 August 2011
 *       Modified zone output to produce OMS-compatible csv format
 *
 *    Version 4.2, 17 October 2012
 *       Added option to read input data in OMS-compatible csv format
 *
 *    Version 4.3, 28 November 2012
 *       1) Added variable dataname to contain name of data type in
 *       OMS-csv data input file (precip, tmax, or tmin)
 *       2) Added zoneseq array to contain index for output of zones
 *       in numerical order
 *       3) Changed keyword in OMS-csv input read from "x" to "easting"
 *       and "y" to "northing"
 *
 *    Version 4.4, 29 November 2012
 *       1) Added variables yr_start, mo_start, dy_start, yr_end,
 *       mo_end, dy_end, read in readcsv.c, to write out in OMS-csv
 *       formatted zone output
 *       2) Renamed function getline.c to getln.c to avoid conflict
 *       with built-in function "getline" in gcc compiler
 *
 *    Version 4.5, 5 December 2012:
 *       Added reading of missing data code in OMS-csv input
 *
 *    Version 4.6, 14 February 2013:
 *       Added recognition of NODATA in ARC/INFO elevation grid and
 *       assignment of use flag to those that are non-NODATA
 *
 *    Version 4.6.1, 10 July 2013
 *       Slight changes to OMS-csv reading and output formats
 *
 *    Version 4.7, 13 March 2014
 *       Added capability for hourly, monthly, and yearly data, in
 *       addition to existing daily data, by handling time periods
 *       more generically; added variable mtper to support this
 *
 *    Version 4.8, 15 May 2015
 *       Various clean-up and minor modifications, including:
 *       - Remove -g command line switch (GeoEAS output -- obsolete)
 *       - Minor output formatting and wording clarifications
 *       More significant changes:
 *       - Added -f command line switch to allow reading of kriging
 *         weights from file instead of calculating them
 *       - Changed readgrid.c to exclude from use any grid cells that
 *         have a valid elevation but are outside of the watershed
 *         and zone grids
 *       Regarding Java user interface:
 *          This has not been maintained and at this point is no
 *          longer supported.  The previous double backslash requirement
 *          for filenames in config file (a Java user interface
 *          idiosyncracy) is therefore no longer necessary.
 *          
 */

#include "malloc.h" 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "dk_m.h"                /* program header with array dimensions */

double **a;                      /* data matrix for solving for kriging
                                    weights (input to m_inv()) */
float accum;                     /* accumulated precip code */
float **ad;                      /* matrix of distances between prec/temp
                                    stations for computing kriging weights */
float *adata;                    /* vector of aggregated data */
struct {
	int cols;                     /* number of columns in ARC/INFO raster */
	int rows;                     /* number of rows in ARC/INFO raster */
	double xll;                   /* x-coordinate of lower left corner of
                                    ARC/INFO raster */
	double yll;                   /* y-coordinate of lower left corner of
                                    ARC/INFO raster */
	float cell;                   /* cellsize of ARC/INFO raster */
	float nodata;                 /* nodata value of ARC/INFO raster */
} arc;
void arcout();                   /* function to write out daily grids in
                                    ARC/INFO format */
float **b0, **b1;                /* matrices of regression intercepts
                                    and slopes */
float **b02, **b12;              /* matrices of intercepts and slopes for line
                                    between highest zero and lowest nonzero
                                    swe values */
double b0dum, b1dum;             /* temporary intercept and slope variables */
char dataname[21];               /* name of data type in csv input file
                                    (precip, tmax, or tmin) */
/* int dayfrac;                     day fraction of data
                                    (beginning of time period) */
float **dgrid;                   /* matrix of distances between grid cells
                                    and prec/temp stations */
float dist_en();                 /* function to calculate distances between
                                    stations based on eastings and northings */
float dist_ll();                 /* function to calculate distances between
                                    stations based on latitude and longitude */
double **dmatrix();              /* double matrix space allocation function */
int dpp;                         /* days (time steps) per period */
int dppl;                        /* days (time steps) in last period */
int dstop;                       /* stopping day (time step) for storm index */
float dum;                       /* intermediate calculation variable */
double *dvector();               /* double vector space allocation function */
int dy_end;                      /* ending day of OMS-csv input file */
int dy_start;                    /* starting day of OMS-csv input file */
float *elevations;				 /* vector of elevation for each station */
double exp();                    /* exponential function */
int *firstday;                   /* vector of first day (period) of data for each year */
//FILE *fopen();                   /* file open function */
FILE *fpin1;                     /* pointer to data input file */
FILE *fpin2;                     /* pointer to elevation grid file */
FILE *fpin3;                     /* pointer to watershed mask grid file */
FILE *fpin4;                     /* pointer to zone grid file */
FILE *fpout;                     /* pointer to main output file */
FILE *fpzone;                    /* pointer to zone output file */
FILE *fpkw;                      /* pointer to kriging weight file */
int getln();                     /* function to read line from file */
float *gprec;                    /* vector of precip at grid cells for
                                    one day */
struct {
	double north;                 /* northernmost extent of GRASS raster */
	double south;                 /* southernmost extent of GRASS raster */
	double east;                  /* easternmost extent of GRASS raster */
	double west;                  /* westernmost extent of GRASS raster */
	int rows;                     /* number of rows in GRASS raster */
	int cols;                     /* number of columns in GRASS raster */
	double nsres;                 /* north-south resolution */
	double ewres;                 /* east-west resolution */
} grass;
void grassout();                 /* function to write out daily grids in
                                    GRASS format */
struct {
	float east;                   /* easting (or longitude) of grid cell */
	float north;                  /* northing (or latitude) of grid cell */
	float elev;                   /* elevation of grid cell (thousands) */
	int mask;                     /* 1 = cell is in watershed, 0 = outside */
	int use;                      /* 1 = cell is used (valid elev.), 0 = not used */
	int zone;                     /* zone number */
} grid[MGRID];
int icoord = 0;                  /* coordinate and grid flag
                                    (1 = lat. and long., column format;
                                     2 = easting and northing, column format;
                                     3 = easting and northing, GRASS format) */
int igridout1, igridout2;        /* beginning and ending day (period) numbers for
                                    grid output */
int igridpr;                     /* flag for precision of values in GRASS or
                                    ARC/INFO grid output files (1=float, one
                                    decimal place; 2=integer, 3=integer*10 */
int imask;                       /* flag indicating if watershed mask is to be
                                    used for calculating spatial averages */
int **imatrix();                 /* int matrix space allocation function */ 
int imiss;                       /* flag to indicate if one or more stations
                                    have missing data */
int iout = 0;                    /* output format (1 = tabular,
                                    2 = GRASS+tabular, 3 = ARC/INFO+tabular,
                                    4 = IPW+tabular) */
void ipwout();                   /* function to write out daily grids in
                                    IPW format */
int ireg = 0;                    /* flag to request printout of regressions */
int irmeth;                      /* regression method flag:
                                    1 = least squares regression
                                    2 = least absolute deviations */
int **iswehz;                    /* index of station with highest zero swe */
int **isweln;                    /* index of station with lowest nonzero swe */
int *ivector();                  /* int vector space allocation function */
int iwt;                         /* station weighting flag (1 = distance
                                    weighting; 2 = equal weighting) */
int izero;                       /* flag indicating a day where all stations
                                    have zero precipitation (izero = 1) */
int izone;                       /* flag indicating if zones (such as hydrologic 
                                    response units) are to be defined */
double *krige();                    /* kriging function */
int *lastday;                    /* vector of last day (period) of data for each year */
int len;                         /* string length */
char line[501];                  /* input line buffer */
int luret;                       /* return value from lusolv() */
int lusolv();                    /* linear equation solver - LU decomposition */
double mae;                      /* mean absolute error */
float **map;                     /* mean areal prec/temp matrix */
float **matrix();                /* float matrix space allocation function */
int medfit();                    /* least absolute deviations regression
                                    function */
float missing;                   /* missing data value (9999.8 internally) */
int mo_end;                      /* ending month of OMS-csv input file */
int mo_start;                    /* starting month of OMS-csv input file */
int mtper;                       /* maximum number of time periods in a year
                                    (8784 for hourly data; 366 for daily data;
                                    12 for monthly data; 1 for yearly data) */
int netcdfout();				 /* NETCDF output function */
int ngrid;                       /* number of grid cells */
int ngriduse;                    /* number of grid cells used (non-missing) */
int nmask;                       /* number of grid cells within mask */
int nper;                        /* number of periods */
int nperm1;                      /* nper minus 1 */
int nsta;                        /* number of stations */
int nstop;                       /* stopping value for loop index n */
int nstorm = 0;                  /* number of storms */
int nyear;                       /* number of years of data */
int nzone;                       /* number of zones */
double pow();                    /* power function */
double r;                        /* correlation coefficient */
void readcsv();                  /* function to read input data in OMS-csv
                                    format */
void readdata();                 /* function to read input data */
void readgrid();                 /* function to read grid data */
float replace;                   /* code to replace an accumulated precip
                                    value */
int ret;                         /* function return code */
double se;                       /* standard error */
float **snolin;                  /* snowline */
int sreg();                      /* simple linear regression function */
struct stations {
	char id[26];                  /* station identifier */
	float elev;                   /* elevation (thousands) */
	float east;                   /* easting (or longitude) of station */
	float north;                  /* northing (or latitude) of station */
	float **data;                 /* data matrix */
} sta[MSTA];
//int *staflg;                     /* station use flags */
struct {
	int dstart;                   /* index of starting day of storm */
	int ystart;                   /* index of starting year of storm */
	int slen;                     /* storm length (days) */
} storm[MSTORM];
double t;                        /* t-statistic */
int type;                        /* data type (1 = prec, 2 = temp, 3 = swe,
                                    4 = other) */
float *vector();                 /* float vector space allocation function */
double *w;                       /* kriging weights */
float **wall;                    /* kriging weight matrix for all stations */
double *x, *y;                   /* regression data vectors */
int *year;                       /* years of data */
int yr_end;                      /* ending year of OMS-csv input file */
int yr_start;                    /* starting year of OMS-csv input file */
struct {
	int number;                   /* zone number */
	int ncells;                   /* number of grid cells in zone */
	float mean;                   /* mean value of prec/temp/swe within zone */
} zone[MZONE];
void zoneout();                  /* function to compute and write out
                                    zonal means */
int zoneseq[MZONE];              /* array index number in zone structure
                                    used to produce zone output in numerical
                                    zone order */

int iomscsv = 0;
int iprintdistances = 0;
int iprintelevation = 0;
int iprintinput = 0;
int iprintresiduals = 0;
int iprintweights = 0;
int i_input_to_output = 0;
int ioutputdir = 0;
int nthreads = 1;
int use_config_file = 0;
int ikwfile = 0;
char config_filename[150];
char kw_filename[150];
char outputdir[150];

char infile[101];             /* input file name */
char elevfile[101];           /* elevation file name */
char maskfile[101];           /* mask file name */
char zonefile[101];           /* zone file name */
int istorm = 0;               /* flag for storm option */


void get_interactive_configuration();
void get_file_configuration(); 


int main(argc, argv)
int argc;
char *argv[];
{
	double atof();                /* ascii-to-float function */
	int atoi();                   /* ascii-to-int function */
	float ewdist;                 /* east-west distance -- argument to
                                    dist_ll() (not used here) */
	int i, j, k, m;               /* loop indexes and counters */
	float nsdist;                 /* north-south distance -- argument to
                                    dist_ll() (not used here) */
	int nstap1;                   /* nsta plus 1 */
	void period1();               /* prec/temp vs. elev calculation function
                                    for periods */
	void period2();               /* MAP/MAT calculation function for periods */
	void storm1();                /* prec vs. elev calculation function
                                    for storms */
	void storm2();                /* MAP calculation function for storms */
	void swe1();                  /* swe vs. elevation calculation function */
	void swe2();                  /* MASWE calculation function */


	/* First, evaluate command-line options and set flags accordingly */

	if (argc > 1) {
		for (i = 1; i < argc; i++) {
			if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "/o") == 0)
				iomscsv = 1;
			else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "/r") == 0)
				iprintresiduals = 1;
			else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "/d") == 0)
				iprintdistances = 1;
			else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "/e") == 0)
				iprintelevation = 1;
			else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "/i") == 0)
				iprintinput = 1;
			else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "/w") == 0)
				iprintweights = 1;
			else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "/c") == 0)
				i_input_to_output = 1;
			else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "/t") == 0) {
				if (sscanf (argv[i+1], "%i", &nthreads) !=1 ) {
					printf ("ERROR - t option not an integer\n");
					exit(0);
				}
				if (nthreads > omp_get_max_threads()){
					nthreads = omp_get_max_threads();
					printf("WARNING - maximum number of threads is %i, using %i\n", omp_get_max_threads(), nthreads);
				}
			}
			else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "/k") == 0) {
				use_config_file = 1;
				strcpy(config_filename, argv[i+1]);
			}
			else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "/f") == 0) {
				ikwfile = 1;
				strcpy(kw_filename, argv[i+1]);
			}
		}
	}

	if (use_config_file)
		get_file_configuration(config_filename);
	else
		get_interactive_configuration();

	/* Read input data */

	if (iomscsv == 1) {
		printf("\n\nNow reading input data in csv format ...\n");
		readcsv();
	}
	else {
		printf("\n\nNow reading input data in column format ...\n");
		readdata();
	}

	/* Write input data to output and exit, if requested (-c switch) */

	if (i_input_to_output == 1) {
		fprintf(fpout, "\n\nStation, Elevation, Easting, Northing:\n");
		for (i = 0; i < nsta; i++)
			fprintf(fpout, "\n%s,   %5.0f,   %8.2f,   %8.2f", sta[i].id,
					sta[i].elev*1000, sta[i].east, sta[i].north);
		fprintf(fpout, "\n\n\nData:\n");
		for (k = 0; k < nyear; k++) {
			j = 0;
			while (j < mtper) {
				for (i = 0; i < nsta; i++)
					if (sta[i].data[j][k] < missing)
						break;
				if (i < nsta) {
					fprintf(fpout, "%4d %4d", year[k], j+1);
					for (i = 0; i < nsta; i++)
						fprintf(fpout, "%8.2f", sta[i].data[j][k]);
					fprintf(fpout, "\n");
				}
				j++;
			}
		}
		exit(0);
	}

	/* Print out input data in main output file, if requested (-i switch) */

	if (iprintinput == 1) {
		fprintf(fpout, "\n\nStation, Elevation, Easting, Northing:\n");
		for (i = 0; i < nsta; i++)
			fprintf(fpout, "\n%s,   %5.0f,   %8.2f,   %8.2f", sta[i].id,
					sta[i].elev*1000, sta[i].east, sta[i].north);
		fprintf(fpout, "\n\n\nData:\n");
		for (k = 0; k < nyear; k++) {
			if (type == 1)
				fprintf(fpout, "\n\n\nPrecipitation");
			else if (type == 2)
				fprintf(fpout, "\n\n\nTemperature");
			else if (type == 3)
				fprintf(fpout, "\n\n\nSnow water equivalent");
			fprintf(fpout, " data for year %d:\n\nPeriod", year[k]);
			for (i = 0; i < nsta; i++)
				fprintf(fpout, "%8d", i+1);
			j = 0;
			while (j < mtper) {
				for (i = 0; i < nsta; i++)
					if (sta[i].data[j][k] < missing)
						break;
				if (i < nsta) {
					fprintf(fpout, "\n%4d", j+1);
					for (i = 0; i < nsta; i++)
						fprintf(fpout, "%8.2f", sta[i].data[j][k]);
				}
				j++;
			}
		}
		fprintf(fpout, "\n");
	}

	/* Read grid data */

	printf("\nNow reading grid data ...\n"); fflush(stdout);
	readgrid();
	/* Debug
fprintf(fpout, "\n\nGrid data:\n");
for (i = 0; i < ngrid; i++) {
   fprintf(fpout, "\n%6.2f   %6.2f   %5.0f", grid[i].north, grid[i].east,
           grid[i].elev*1000);
}
   fflush(fpout);
   End debug */

	/* Allocate array space */

	nstap1 = nsta + 1;
//	a = dmatrix(nstap1, nsta+2);
	ad = matrix(nsta, nsta);
	adata = vector(nsta);
	if (istorm == 1) {
		b0 = matrix(MSTORM, 1);
		b1 = matrix(MSTORM, 1);
	}
	else {
		b0 = matrix(nper, nyear);
		b1 = matrix(nper, nyear);
	}
	dgrid = matrix(ngrid, nsta);
	elevations = vector(nsta);
	gprec = vector(ngrid);
	map = matrix(mtper, nyear);
	//	staflg = ivector(nsta);
//	w = dvector(nstap1);
	wall = matrix(ngrid, nsta);
	x = dvector(nsta);
	y = dvector(nsta);
	if (type == 3) {
		b02 = matrix(nper, nyear);
		b12 = matrix(nper, nyear);
		iswehz = imatrix(nper, nyear);
		isweln = imatrix(nper, nyear);
		snolin = matrix(nper, nyear);
	}

	/* Initialize b0, b1, and map matrices to missing code */

	if (istorm == 1) {
		for (i = 0; i < MSTORM; i++)
			b0[i][0] = b1[i][0] = 99999;
	}
	else {
		for (i = 0; i < nper; i++)
			for (j = 0; j < nyear; j++)
				b0[i][j] = b1[i][j] = 99999;
	}
	for (i = 0; i < mtper; i++)
		for (j = 0; j < nyear; j++)
			map[i][j] = (float) (missing + 0.1);

	/* Set MAP and MASWE to zero for days where prec or swe for all stations
      is zero */

	if (type == 1 || type == 3) {
		for (k = 0; k < nyear; k++) {
			for (j = 0; j < mtper; j++) {
				izero = 0;
				for (i = 0; i < nsta; i++) {
					if (sta[i].data[j][k] < missing) {
						izero = 1;
						if (sta[i].data[j][k] > 0.001) {
							izero = 0;
							break;
						}
					}
				}
				if (izero == 1) {
					map[j][k] = 0;
					/* Debug
printf("\nMAP for period %d year %d = %8.4f\n", j+1, year[k], map[j][k]);
   fflush(stdout);
   End debug */
				}
			}
		}
	}

	/* Read or calculate kriging weights */

	if (iwt == 1) {

		if (ikwfile == 1) {

			/* If specified by command line switch, read kriging weights from
            file instead of calculating them */

			printf("\nNow reading kriging weights ...\n");

			fpkw = fopen(kw_filename, "r");

			/* Find spot in file where kriging weights start */

			m = 0;
			while (getln(line, fpkw) != EOF) {
				if (strstr(line, "Kriging weights") != NULL) {
					m = 1;
					break;
				}
			}
			if (m == 0) {
				printf("\nKriging weights not found in file %s.\nProgram terminated ...\n",
						kw_filename);
				exit(0);
			}

			/* Read grid cell number and station weights */

			while (fscanf(fpkw, "%d", &i) > 0) {
				for (j = 0; j < nsta; j++)
					fscanf(fpkw, "%f", &wall[i-1][j]);
			}
		}

		else {

			/* Calculate kriging weights */

			printf("\nNow calculating kriging weights ...\n");

			/* Compute distances between stations and load distances into
            ad matrix for later use in solving linear system for kriging weights */

			for (i = 0; i < nsta; i++) {
				ad[i][i] = 0;
				elevations[i] = sta[i].elev;
				for (j = i+1; j < nsta; j++) {
					if (icoord == 1)
						ad[i][j] = ad[j][i] = dist_ll(sta[i].north, sta[i].east,
								sta[j].north, sta[j].east,
								&ewdist, &nsdist);
					else
						ad[i][j] = ad[j][i] = dist_en(sta[i].north, sta[i].east,
								sta[j].north, sta[j].east);
				}
			}

			/* Compute distances between grid cells and prec/temp/swe stations */

			for (i = 0; i < ngrid; i++) {
				if (grid[i].use == 1) {
					for (j = 0; j < nsta; j++) {
						if (icoord == 1)
							dgrid[i][j] = dist_ll(grid[i].north, grid[i].east,
									sta[j].north, sta[j].east, &ewdist,
									&nsdist);
						else
							dgrid[i][j] = dist_en(grid[i].north, grid[i].east,
									sta[j].north, sta[j].east);
					}
				}
			}

			if (iprintdistances == 1) {
				/* print out distances among stations and grid cells */
				fprintf(fpout, "\n\n\nMatrix of distances between stations (km):\n");
				for (i = 0; i < nsta; i++) {
					fprintf(fpout, "\n");
					for (j = 0; j < nsta; j++)
						fprintf(fpout, "%9.2f", ad[i][j]);
				}
				fprintf(fpout, "\n\n\n%s\n",
						"Distances between grid cells and prec/temp/swe stations (km):");
				for (i = 0; i < ngrid; i++) {
					if (grid[i].use == 1) {
						fprintf(fpout, "\n%d", i+1);
						for (j = 0; j < nsta; j++)
							fprintf(fpout, "%9.2f", dgrid[i][j]);
					}
				}
			}

			/* Calculate kriging weights using all stations */
			omp_set_dynamic(0);     // Explicitly disable dynamic teams
			omp_set_num_threads(nthreads); // Use N threads for all consecutive parallel regions

			#pragma omp parallel shared(nsta, ad, dgrid, elevations, grid) private(i, j, w)
			{
				#pragma omp for
				for (i = 0; i < ngrid; i++) {

					if (grid[i].use == 1) {

						w = krige(i, nsta, ad, dgrid, elevations);

//						#pragma omp critical
//						{
						for (j = 0; j < nsta; j++){
							wall[i][j] = (float) w[j];
						}
//						}
					}
				}
			}
		}
	}

	/* For equal weighting, set weights equal to 1/nsta */

	else if (iwt == 2) {
		dum = (float) (1. / nsta);
		for (i = 0; i < ngrid; i++)
			if (grid[i].use == 1) {
				for (j = 0; j < nsta; j++)
					wall[i][j] = dum;
			}
	}

	if (iprintweights == 1) {
		/* Print out weights */
		fprintf(fpout, "\n\n\nGrid\nPt.:   Kriging weights:\n");
		for (i = 0; i < ngrid; i++) {
			if (grid[i].use == 1) {
				fprintf(fpout, "\n%d", i+1);
				for (j = 0; j < nsta; j++)
					fprintf(fpout, "%8.4f", wall[i][j]);
			}
		}
		fprintf(fpout, "\n\n");
	}

	/* If requested with command line switch, set flag for printing out
      elevation regression results */

	if (iprintelevation == 1)
		ireg = 1;


	/* Before starting computations, write some run information to
      main output file */

	fprintf(fpout, "Detrended Kriging (DK) Program\n\n");
	fprintf(fpout, "Input file:      %s\n", infile);
	fprintf(fpout, "Elevation file:  %s\n", elevfile);
	if (imask == 1)
		fprintf(fpout, "Mask file:       %s\n", maskfile);
	if (izone == 1)
		fprintf(fpout, "Zone file:       %s", zonefile);

	/* Also write some header information to zone output file */

	if (izone == 1) {
		fprintf(fpzone, "@T,obs\ndate_start, %d %d %d 0 0 0",
				yr_start, mo_start, dy_start);
		fprintf(fpzone, "\ndate_end, %d %d %d 0 0 0",
				yr_end, mo_end, dy_end);
		fprintf(fpzone, "\ndate_format, yyyy MM dd H m s\n@H,date");
		for (j = 0; j < nzone; j++)
			fprintf(fpzone, ",%s[%d]", dataname, j);
		fprintf(fpzone, "\ntype,Date");
		for (j = 0; j < nzone; j++)
			fprintf(fpzone, ",Real");
		fprintf(fpzone, "\n");
	}

	/* For detrending, compute regressions for each period and year
      or for each storm then compute residuals */

	if (istorm == 1) {
		printf("\nNow calculating prec-elevation regressions");
		printf(" by analyzing storms ...\n");
		storm1();
	}
	else if (type == 3) {
		printf("\nNow calculating swe-elevation regressions ...\n");
		swe1();
	}
	else {
		printf("\nNow calculating parameter-elevation regressions");
		printf(" by analyzing %d-time-step periods ...\n", dpp);
		period1();
	}


	if (iprintresiduals == 1) {
		/* print out detrended residuals */
		for (i = 0; i < nsta; i++) {
			fprintf(fpout, "\n\n%s%s, %5.0f, %6.2f, %6.2f:\n\n%s",
					"Detrended data for ", sta[i].id, (sta[i].elev*1000),
					sta[i].north, sta[i].east, "Period");
			for (k = 0; k < nyear; k++)
				fprintf(fpout, "%10d", year[k]);
			fprintf(fpout, "\n");
			for (j = 0; j < mtper; j++) {
				fprintf(fpout, "\n%4d", j+1);
				for (k = 0; k < nyear; k++) {
					if (sta[i].data[j][k] < accum)
						fprintf(fpout, "%10.6f", sta[i].data[j][k]);
					else
						fprintf(fpout, "          ");
				}
			}
		}
		fprintf(fpout, "\n\n");
	}


	/* For each day, compute kriging weights (if there are stations with
      missing data), estimate grid cell values, and compute areal averages */

	if (istorm == 1) {
		printf("\nNow calculating grid cell precipitation ...\n");
		storm2();
	}
	else if (type == 1) {
		printf("\nNow calculating grid cell precipitation ...\n");
		period2();
	}
	else if (type == 2) {
		printf("\nNow calculating grid cell temperature ...\n");
		period2();
	}
	else if (type == 3) {
		printf("\nNow calculating grid cell snow water equivalent ...\n");
		swe2();
	}
	else {
		printf("\nNow calculating grid cell values ...\n");
		period2();
	}

	/* Write out results in tabular format */

	if (type == 1)
		fprintf(fpout, "\nMean areal precipitation:\n\nTime\nPeriod");
	else if (type == 2)
		fprintf(fpout, "\nMean areal temperature:\n\nTime\nPeriod");
	else if (type == 3)
		fprintf(fpout, "\nMean areal snow water equivalent:\n\nTime\nPeriod");
	else
		fprintf(fpout, "\nMean areal values:\n\nTime\nPeriod");
	for (k = 0; k < nyear; k++)
		fprintf(fpout, "%8d", year[k]);
	fprintf(fpout, "\n");
	for (j = 0; j < mtper; j++) {
		fprintf(fpout, "\n%4d  ", j+1);
		for (k = 0; k < nyear; k++) {
			if (map[j][k] >= missing)
				fprintf(fpout, "        ");
			else
				fprintf(fpout, "%8.2f", map[j][k]);
		}
	}
	fprintf(fpout, "\n");

	return 0;
}





/*
 *
 */
void get_interactive_configuration() {

	/* Get input and output file names and open files */

	printf("\n\n\n\n\n%s\n%s\n\n\n%s\n%s",
			"Detrended Kriging Program (DK)",
			"Version 4.8     15 May 2015",
			"Please enter input data file name (q to quit):",
			"   ==> ");
	gets(infile);
	if ((len = strlen(infile)) == 0 || (len == 1 && (infile[0] == 'q' ||
			infile[0] == 'Q')))
		exit(0);
	if ((fpin1 = fopen(infile, "r")) == NULL) {
		printf("\n\nError opening file %s.\nProgram terminated ...\n", infile);
		exit(0);
	}

	printf("\n\n%s\n%s\n%s\n%s\n%s\n%s",
			"What type of data are in this file (q to quit)?",
			"   1) Precipitation",
			"   2) Temperature",
			"   3) Snow water equivalent",
			"   4) Other (no restrictions on detrending)",
			"   ==> ");
	gets(line);
	if ((len = strlen(line)) == 0 || (len == 1 && (line[0] == 'q' ||
			line[0] == 'Q')))
		exit(0);
	if (line[0] == '1') {
		type = 1;
		missing = 9999.8f;
		accum = 8888.7f;
		replace = 8999.9f;
	}
	else if (line[0] == '2') {
		type = 2;
		missing = accum = 9999.8f;
	}
	else if (line[0] == '3') {
		type = 3;
		missing = accum = 9999.8f;
	}
	else if (line[0] == '4') {
		type = 4;
		missing = accum = 9999.8f;
	}
	else {
		printf("\n\nInvalid response ... program terminated ...\n");
		exit(0);
	}

	printf("\n\n%s\n%s\n%s\n%s\n%s\n%s",
			"What is the time step of the data in this file (q to quit)?",
			"   1) Hourly",
			"   2) Daily",
			"   3) Monthly",
			"   4) Yearly",
			"   ==> ");
	gets(line);
	if ((len = strlen(line)) == 0 || (len == 1 && (line[0] == 'q' ||
			line[0] == 'Q')))
		exit(0);
	if (line[0] == '1')
		mtper = 8784;
	else if (line[0] == '2')
		mtper = 366;
	else if (line[0] == '3')
		mtper = 12;
	else if (line[0] == '4')
		mtper = 1;
	else {
		printf("\n\nInvalid response ... program terminated ...\n");
		exit(0);
	}

	while (icoord == 0) {
		printf("\n\n%s%s\n%s\n%s\n%s\n%s\n%s",
				"What is the coordinate system and ",
				"elevation grid format (q to quit)?",
				"   1) Latitude and longitude, column grid format",
				"   2) Easting and northing, column grid format",
				"   3) Easting and northing, GRASS grid format",
				"   4) Easting and northing, ARC/INFO grid format",
				"   ==> ");
		gets(line);
		if (line[0] == 'q' || line[0] == 'Q')
			exit(0);
		if (line[0] == '1')
			icoord = 1;
		else if (line[0] == '2')
			icoord = 2;
		else if (line[0] == '3')
			icoord = 3;
		else if (line[0] == '4')
			icoord = 4;
		else
			printf("\n\nInvalid response ... try again ...\n");
	}

	printf("\n\nPlease enter elevation grid file name (q to quit):\n   ==> ");
	gets(elevfile);
	if ((len = strlen(elevfile)) == 0 || (len == 1 && (elevfile[0] == 'q' ||
			elevfile[0] == 'Q')))
		exit(0);
	if ((fpin2 = fopen(elevfile, "r")) == NULL) {
		printf("\n\nError opening file %s.\nProgram terminated ...\n", elevfile);
		exit(0);
	}

	printf("\n\n%s\n%s\n%s\n%s",
			"Please enter watershed mask file name --",
			"GRASS or ARC/INFO format only, must be same as elevation grid",
			"(q to quit, blank for no mask):",
			"   ==> ");
	gets(maskfile);
	if ((len = strlen(maskfile)) == 1 && (maskfile[0] == 'q' ||
			maskfile[0] == 'Q'))
		exit(0);
	else if (len == 0)
		imask = 0;
	else {
		if ((fpin3 = fopen(maskfile, "r")) == NULL) {
			printf("\n\nError opening file %s.\nProgram terminated ...\n",
					maskfile);
			exit(0);
		}
		imask = 1;
		nmask = 0;
	}

	printf("\n\n%s\n%s\n%s\n%s",
			"Please enter zone grid file name --",
			"GRASS or ARC/INFO format only, must be same as elevation grid",
			"(q to quit, blank for no mask):",
			"   ==> ");
	gets(zonefile);
	if ((len = strlen(zonefile)) == 1 && (zonefile[0] == 'q' ||
			zonefile[0] == 'Q'))
		exit(0);
	else if (len == 0)
		izone = 0;
	else {
		if ((fpin4 = fopen(zonefile, "r")) == NULL) {
			printf("\n\nError opening file %s.\nProgram terminated ...\n",
					zonefile);
			exit(0);
		}
		izone = 1;
	}

	while (iout == 0) {
		printf("\n\n%s\n%s\n%s\n%s\n%s\n%s",
				"Please enter desired output format (q to quit):",
				"   1) Mean areal values -- tabular (one column per year)",
				"   2) Timestep grids -- GRASS format, plus (1) above",
				"   3) Timestep grids -- ARC/INFO format, plus (1) above",
				"   4) Timestep grids -- IPW format, plus (1) above",
				"   ==> ");
		gets(line);
		if (line[0] == 'q' || line[0] == 'Q')
			exit(0);
		if (line[0] == '2') {
			if (icoord != 3) {
				printf("\n\n%s%s\n%s\n",
						"Must have GRASS elevation grid format ",
						"for GRASS output format ...",
						"try again ...");
				continue;
			}
			iout = 2;
			printf("\n\n%s\n%s",
					"Please enter beginning period number for GRASS output",
					"   ==> ");
			gets(line);
			sscanf(line, "%d", &igridout1);
			igridout1--;
			printf("\n\n%s\n%s",
					"Please enter ending period number for GRASS output",
					"   ==> ");
			gets(line);
			sscanf(line, "%d", &igridout2);
			igridout2--;
			printf("\n\n%s\n%s\n%s\n%s\n%s",
					"Please enter desired GRASS output precision (q to quit):",
					"   1) Floating point, one decimal place (default)",
					"   2) Integer",
					"   3) Integer, values multiplied by 10",
					"   ==> ");
			gets(line);
			if (line[0] == 'q' || line[0] == 'Q')
				exit(0);
			if (line[0] == '2')
				igridpr = 2;
			else if (line[0] == '3')
				igridpr = 3;
			else
				igridpr = 1;
		}
		else if (line[0] == '3') {
			if (icoord != 4) {
				printf("\n\n%s%s\n%s\n",
						"Must have ARC/INFO elevation grid format ",
						"for ARC/INFO output format ...",
						"try again ...");
				continue;
			}
			iout = 3;
			printf("\n\n%s\n%s",
					"Please enter beginning period number for ARC/INFO output",
					"   ==> ");
			gets(line);
			sscanf(line, "%d", &igridout1);
			igridout1--;
			printf("\n\n%s\n%s",
					"Please enter ending period number for ARC/INFO output",
					"   ==> ");
			gets(line);
			sscanf(line, "%d", &igridout2);
			igridout2--;
			printf("\n\n%s\n%s\n%s\n%s\n%s",
					"Please enter desired ARC/INFO output precision (q to quit):",
					"   1) Floating point, one decimal place (default)",
					"   2) Integer",
					"   3) Integer, values multiplied by 10",
					"   ==> ");
			gets(line);
			if (line[0] == 'q' || line[0] == 'Q')
				exit(0);
			if (line[0] == '2')
				igridpr = 2;
			else if (line[0] == '3')
				igridpr = 3;
			else
				igridpr = 1;
		}
		else if (line[0] == '4') {
			if (icoord != 3 && icoord != 4) {
				printf("\n\n%s%s\n%s\n",
						"Must have GRASS or ARC/INFO elevation grid format ",
						"for IPW output format ...",
						"try again ...");
				continue;
			}
			iout = 4;
			printf("\n\n%s\n%s",
					"Please enter beginning period number for IPW output",
					"   ==> ");
			gets(line);
			sscanf(line, "%d", &igridout1);
			igridout1--;
			printf("\n\n%s\n%s",
					"Please enter ending period number for IPW output",
					"   ==> ");
			gets(line);
			sscanf(line, "%d", &igridout2);
			igridout2--;
		}
		else
			iout = 1;
	}

	printf("\n\nPlease enter main output file name (q to quit):\n   ==> ");
	gets(line);
	if ((len = strlen(line)) == 0 || (len == 1 && (line[0] == 'q' ||
			line[0] == 'Q')))
		exit(0);
	if ((fpout = fopen(line, "w")) == NULL) {
		printf("\n\nError opening file %s.\nProgram terminated ...\n", line);
		exit(0);
	}

	if (izone == 1) {
		printf("\n\nPlease enter zone output file name (q to quit):\n   ==> ");
		gets(line);
		if ((len = strlen(line)) == 0 || (len == 1 && (line[0] == 'q' ||
				line[0] == 'Q')))
			exit(0);
		if ((fpzone = fopen(line, "w")) == NULL) {
			printf("\n\nError opening file %s.\nProgram terminated ...\n", line);
			exit(0);
		}
	}

	printf("\n\n%s%s\n%s\n%s\n%s",
			"Please select regression method for elevation detrending",
			" (q to quit):",
			"   1) Least squares regression",
			"   2) Least absolute deviations regression",
			"   ==> ");
	gets(line);
	if (line[0] == 'q' || line[0] == 'Q')
		exit(0);
	if (line[0] == '2')
		irmeth = 2;
	else
		irmeth = 1;

	printf("\n\n%s\n%s\n%s\n%s",
			"Please select station weighting method (q to quit):",
			"   1) Distance weighting (linear semivariogram)",
			"   2) Equal weighting",
			"   ==> ");
	gets(line);
	if (line[0] == 'q' || line[0] == 'Q')
		exit(0);
	if (line[0] == '2')
		iwt = 2;
	else
		iwt = 1;

	if (type == 1)
		printf("\n\n%s\n%s\n%s",
				"Please enter the number of time steps per period you desire",
				"or enter 0 or s to use storm periods (q to quit):",
				"   ==> ");
	else
		printf("\n\n%s%s\n%s",
				"Please enter the number of time steps per period you desire",
				" (q to quit):", "   ==> ");
	gets(line);
	if ((len = strlen(line)) == 0 || (len == 1 && (line[0] == 'q' ||
			line[0] == 'Q')))
		exit(0);
	if (line[0] == '0' || line[0] == 's' || line[0] == 'S')
		istorm = 1;
	else {
		dpp = atoi(line);
		nper = mtper / dpp;
	}
}
