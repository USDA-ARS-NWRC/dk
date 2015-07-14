extern double **a;               /* data matrix for solving for kriging
                                    weights (input to m_inv()) */
extern float accum;              /* accumulated precip code */
extern float **ad;               /* matrix of distances between prec/temp
                                    stations for computing kriging weights */
extern float *adata;             /* vector of aggregated data */
extern struct {
   int cols;                     /* number of columns in ARC/INFO raster */
   int rows;                     /* number of rows in ARC/INFO raster */
   double xll;                   /* x-coordinate of lower left corner of
                                    ARC/INFO raster */
   double yll;                   /* y-coordinate of lower left corner of
                                    ARC/INFO raster */
   float cell;                   /* cellsize of ARC/INFO raster */
   float nodata;                 /* nodata value of ARC/INFO raster */
} arc;
extern void arcout();            /* function to write out daily grids in
                                    ARC/INFO format */
extern float **b0, **b1;         /* matrices of regression intercepts
                                    and slopes */
extern float **b02, **b12;       /* matrices of intercepts and slopes for line
                                    between highest zero and lowest nonzero
                                    swe values */
extern double b0dum, b1dum;      /* temporary intercept and slope variables */
extern char dataname[21];        /* name of data type in csv input file
                                    (precip, tmax, or tmin) */
/* extern int dayfrac;              day fraction of data
                                    (beginning of time period) */
extern float **dgrid;            /* matrix of distances between grid cells
                                    and prec/temp stations */
extern float dist_en();          /* function to calculate distances between
                                    stations based on eastings and northings */
extern float dist_ll();          /* function to calculate distances between
                                    stations based on latitude and longitude */
extern double **dmatrix();       /* double matrix space allocation function */
extern int dpp;                  /* days (time steps) per period */
extern int dppl;                 /* days (time steps) in last period */
extern int dstop;                /* stopping day (time step) for storm index */
extern float dum;                /* intermediate calculation variable */
extern double *dvector();        /* double vector space allocation function */
extern int dy_end;               /* ending day of OMS-csv input file */
extern int dy_start;             /* starting day of OMS-csv input file */
extern double exp();             /* exponential function */
extern int *firstday;            /* vector of first day (period) of data for each year */
extern FILE *fopen();            /* file open function */
extern FILE *fpin1, *fpin2, *fpin3, *fpin4;
                                 /* pointers to input files */
extern FILE *fpout, *fpzone;     /* pointers to output files */
extern int getln();              /* function to read line from file */
extern float *gprec;             /* vector of precip at grid cells for
                                    one day */
extern struct {
   double north;                 /* northernmost extent of GRASS raster */
   double south;                 /* southernmost extent of GRASS raster */
   double east;                  /* easternmost extent of GRASS raster */
   double west;                  /* westernmost extent of GRASS raster */
   int rows;                     /* number of rows in GRASS raster */
   int cols;                     /* number of columns in GRASS raster */
   double nsres;                 /* north-south resolution */
   double ewres;                 /* east-west resolution */
} grass;
extern void grassout();          /* function to write out daily grids in
                                    GRASS format */
extern struct {
   float east;                   /* easting (or longitude) of grid cell */
   float north;                  /* northing (or latitude) of grid cell */
   float elev;                   /* elevation of grid cell (thousands) */
   int mask;                     /* 1 = cell is in watershed, 0 = outside */
   int use;                      /* 1 = cell is used (valid elev.), 0 = not used */
   int zone;                     /* zone number */
} grid[];
extern int icoord;               /* coordinate and grid flag
                                    (1 = lat. and long., column format;
                                     2 = easting and northing, column format;
                                     3 = easting and northing, GRASS format)
                                     4 = easting and northing, ARC/INFO format) */
extern int igridout1, igridout2; /* beginning and ending day (period) numbers for
                                    grid output */
extern int igridpr;              /* flag for precision of values in GRASS or
                                    ARC/INFO grid output files (1=float, one
                                    decimal place; 2=integer, 3=integer*10 */
extern int imask;                /* flag indicating if watershed mask is to be
                                    used for calculating spatial averages */
extern int imiss;                /* flag to indicate if one or more stations
                                    have missing data */
extern int iout;                 /* output format (1 = tabular,
                                    2 = GRASS+tabular, 3 = ARC+tabular,
                                    4 = IPW+tabular) */
extern void ipwout();            /* function to write out daily grids in
                                    IPW format */
extern int isleap();             /* determine if given year is a leap year
                                    (1 = leap year, 0 = regular year) */
extern int **iswehz;             /* index of station with highest zero swe */
extern int **isweln;             /* index of station with lowest nonzero swe */
extern int ireg;                 /* flag to request printout of regressions */
extern int irmeth;               /* regression method flag:
                                    1 = least squares regression
                                    2 = least absolute deviations */
extern int *ivector();           /* int vector space allocation function */
extern int iwt;                  /* station weighting flag (1 = distance
                                    weighting; 2 = equal weighting) */
extern int izero;                /* flag indicating a day where all stations
                                    have zero precipitation (izero = 1) */
extern int izone;                /* flag indicating if zones (such as hydrologic 
                                    response units) are to be defined */
extern double *krige();             /* kriging function */
extern int *lastday;             /* vector of last day (period) of data for each year */
extern int len;                  /* string length */
extern char line[501];           /* input line buffer */
extern int luret;                /* return value from lusolv() */
extern int lusolv();             /* linear equation solver - LU decomposition */
extern double mae;               /* mean absolute error */
extern float **map;              /* mean areal prec/temp matrix */
extern float **matrix();         /* float matrix space allocation function */
extern int medfit();             /* least absolute deviations regression
                                    function */
extern float missing;            /* missing data code
                                    (= 99.99 for prec, = 999 for temp) */
extern int mo_end;               /* ending month of OMS-csv input file */
extern int mo_start;             /* starting month of OMS-csv input file */
extern int mtper;                /* maximum number of time periods in a year
                                    (8784 for hourly data; 366 for daily data;
                                    12 for monthly data; 1 for yearly data) */
extern int netcdfout();			 /* NETCDF output function */
extern int ngrid;                /* number of grid cells */
extern int ngriduse;             /* number of grid cells used (non-missing) */
extern int nmask;                /* number of grid cells within watershed mask */
extern int nper;                 /* number of periods */
extern int nperm1;               /* nper minus 1 */
extern int nsta;                 /* number of stations */
extern int nstop;                /* stopping value for loop index n */
extern int nstorm;               /* number of storms */
extern int nyear;                /* number of years of data */
extern int nzone;                /* number of zones */
extern double pow();             /* power function */
extern double r;                 /* correlation coefficient */
extern void readcsv();           /* function to read input data in OMS-csv
                                    format */
extern void readdata();          /* function to read input data */
extern void readgrid();          /* function to read grid data */
extern float replace;            /* code to replace an accumulated precip
                                    value */
extern int ret;                  /* function return code */
extern double se;                /* standard error */
extern float **snolin;           /* snowline */
extern int sreg();               /* simple linear regression function */
extern struct {
   char id[26];                  /* station identifier */
   float elev;                   /* elevation (thousands) */
   float east;                   /* easting (or longitude) of station */
   float north;                  /* northing (or latitude) of station */
   float **data;                 /* data matrix */
} sta[];
extern int *staflg;              /* station use flags */
extern struct {
   int dstart;                   /* index of starting day of storm */
   int ystart;                   /* index of starting year of storm */
   int slen;                     /* storm length (days) */
} storm[];
extern double t;                 /* t-statistic */
extern int type;                 /* data type (1 = prec, 2 = temp, 3 = swe, 
                                    4 = other) */
extern float *vector();          /* float vector space allocation function */
extern double *w;                /* kriging weights */
extern float **wall;             /* kriging weight matrix for all stations */
extern double *x, *y;            /* regression data vectors */
extern int *year;                /* years of data */
extern int yr_end;               /* ending year of OMS-csv input file */
extern int yr_start;             /* starting year of OMS-csv input file */
extern struct {
   int number;                   /* zone number */
   int ncells;                   /* number of grid cells in zone */
   float mean;                   /* mean value of prec/temp/swe within zone */
} zone[];
extern void zoneout();           /* function to compute and write out
                                    zonal means */
extern int zoneseq[];            /* array index number in zone structure
                                    used to produce zone output in numerical
                                    zone order */
