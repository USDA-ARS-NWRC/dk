/*
 *    netcdfout.c
 *
 *    Scott Havens 07/2015
 *
 *    netcdf_create - creates the netcdf file to start writing to
 *    netcdf_write - Write the data to a netcdf file
 *
 *    A little bit about chunking:
 *    http://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_why_it_matters
 *    http://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
 */

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* This is the name of the data file we will create. */
#define DK_TITLE "Created with DK 4.8 tool"
#define VAR_NAME "variable"
#define CONVENTION "CF-1.6 (http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.pdf)"
#define FILE_NAME "dk_out_%i.nc"

/* We are writing 3D data */
#define NDIMS 3
#define CHUNKSIZE 8096		// size of chunk in bytes
#define VALSIZE 4			// size of value in bytes for a 32-bit float

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


/*
 * Create the netcdf file
 */
int netcdf_create(iy, x, y, nx, ny, ncid)
int iy;                          /* year */
//struct grid *grid;							/* grid structure */
float *x;
float *y;
int nx;								/* number of values in x index */
int ny;								/* number of values in y index */
int *ncid;							/* file id for netcdf file */
{

	/* When we create netCDF variables and dimensions, we get back an
	 * ID for each one. */
	int t_dimid, x_dimid, y_dimid, varid[3];
	int dimids[NDIMS];
	int dimids2[2];
	size_t *chunkSizes;
//	int chunkVals;
	int chunk_size();

	/* indexing and error handling. */
	int retval;

	/* File Name */
	char file_name[20];
	snprintf(file_name, sizeof(file_name), "dk_out_%i.nc", iy);


	/* If there isn't a netcdf file, then create one, else just open it
	 * Create the file. NC_NOCOBBER will not overwrite an existing file
	 */
	retval = nc_create(file_name, NC_CLOBBER | NC_NETCDF4, ncid);
	if ((retval == NC_EEXIST))
	{
		/* File exists so open it instead */
		if ((retval = nc_open(file_name, NC_WRITE, ncid)))
			ERR(retval);

	}
	else if ((retval == NC_NOERR))
	{
		/* Define the dimensions. NetCDF will hand back an ID for each. */
		if ((retval = nc_def_dim(*ncid, "Time", NC_UNLIMITED, &t_dimid)))
			ERR(retval);
		if ((retval = nc_def_dim(*ncid, "y", ny, &y_dimid)))
			ERR(retval);
		if ((retval = nc_def_dim(*ncid, "x", nx, &x_dimid)))
			ERR(retval);

		/* The dimids array is used to pass the IDs of the dimensions of
		 * the variable. */
		dimids[0] = t_dimid;
		dimids[1] = dimids2[0] = y_dimid;
		dimids[2] = dimids2[1] = x_dimid;

		/* define x grid location variable */
		if ((retval = nc_def_var(*ncid, "x", NC_FLOAT, 1, &x_dimid, &varid[0])))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[0], "units", strlen("m"), "m")))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[0], "long_name", strlen("x coordinate of projection"), "x coordinate of projection")))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[0], "_CoordinateAxisType", strlen("GeoX"), "GeoX")))
			ERR(retval);

		/* define y grid location variable */
		if ((retval = nc_def_var(*ncid, "y", NC_FLOAT, 1, &y_dimid, &varid[1])))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[1], "units", strlen("m"), "m")))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[1], "long_name", strlen("y coordinate of projection"), "y coordinate of projection")))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[1], "_CoordinateAxisType", strlen("GeoY"), "GeoY")))
			ERR(retval);

		/* Define dk variable */
		if ((retval = nc_def_var(*ncid, VAR_NAME, NC_FLOAT, NDIMS, dimids, &varid[2])))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[2], "units", strlen("units"), "units")))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, varid[2], "long_name", strlen("variable description"), "variable description")))
			ERR(retval);

		/* Set the chunk size */
		retval = chunk_size(nx, ny, &chunkSizes);
		if ((retval = nc_def_var_chunking(*ncid, varid[2], NC_CHUNKED, chunkSizes)))
			ERR(retval);

		/* Add global attributes for file */
		if ((retval = nc_put_att_int(*ncid, NC_GLOBAL, "Water_Year", NC_INT, 1, &iy)))
			ERR(retval);
		if ((retval = nc_put_att_text(*ncid, NC_GLOBAL, "Title", strlen(DK_TITLE), DK_TITLE)))
			ERR(retval);

		char currentDate[100];
		time_t now = time(NULL);
		struct tm *t = localtime(&now);
		strftime(currentDate, sizeof(currentDate)-1, "%d/%m/%Y %H:%M", t);
		if ((retval = nc_put_att_text(*ncid, NC_GLOBAL, "Created", strlen(currentDate), currentDate)))
			ERR(retval);

		if ((retval = nc_put_att_text(*ncid, NC_GLOBAL, "Conventions", strlen(CONVENTION), CONVENTION)))
			ERR(retval);

		retval = nc_enddef(*ncid);

		/* location to put the coordinate data */
		size_t start[] = {0, 0};
		size_t count[] = {ny, nx};

		if ((retval = nc_put_vara_float(*ncid, varid[0], start, count, x)))
			ERR(retval);

		if ((retval = nc_put_vara_float(*ncid, varid[1], start, count, y)))
			ERR(retval);

	}
	else
		ERR(retval);

	return 0;
}


/*
 * Write the data to the netcdf file
 */
int netcdf_write(ncid, ip, data, nx, ny)
int *ncid;						/* file id for netcdf file */
int ip;                          /* period (sequential number beginning Oct 1) */
float *data;						 /* pointer to data matrix */
int nx;								/* number of values in x index */
int ny;								/* number of values in y index */
{

	int varid;		/* variable id */
	int retval;		/* indexing and error handling. */

	/* get variable id */
	if ((retval = nc_inq_varid(*ncid, VAR_NAME, &varid)))
		ERR(retval);

	/* location to put the data */
	size_t start[] = {ip, 0, 0};
	size_t count[] = {1, ny, nx};

	/*
	 * Write the data to the file
	 */

	//	j = nx*ny - 1;   // j will Point to last Element
	//	i = 0;       // i will be pointing to first element
	//	while (i < j) {
	//		temp = data[i];
	//		data[i] = data[j];
	//		data[j] = temp;
	//		i++;             // increment i
	//		j--;          // decrement j
	//	}

	if ((retval = nc_put_vara_float(*ncid, varid, start, count, data)))
		ERR(retval);

	return 0;
}

/*
 * Determine the "optimized" chunk size for the data. This is based
 * off the blog at
 * http://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
 */
int chunk_size(nx, ny, chunkSize)
int nx;			/* number of x dimension */
int ny;			/* number of y dimension */
int *chunkSize;	/* pointer to chunk size */
{
	int t_chunk = 6;						// 'ideal' number of time chunks, favoring spatial slice
	int chunkVals = CHUNKSIZE/VALSIZE; 		// ideal number of values in a chunk
	int numChunks = sqrt(chunkVals/t_chunk);		// ideal number of chuncks
	int cSize[NDIMS];

//	printf("%i\n",t_chunk);
//	printf("%i\n",chunkVals);
//	printf("%i\n",numChunks);

	cSize[0] = t_chunk;
	cSize[1] = numChunks;
	cSize[2] = numChunks;

	chunkSize = cSize;

	return 0;
}

/*
 * Close the netcdf file when done
 */
int netcdf_close(ncid)
int *ncid;						/* file id for netcdf file */
{
	int retval;		/* indexing and error handling. */

	/* Close the file. This frees up any internal netCDF resources
	 * associated with the file, and flushes any buffers. */
	if ((retval = nc_close(*ncid)))
		ERR(retval);

	return 0;
}
