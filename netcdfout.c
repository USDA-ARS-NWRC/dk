/*
 *    netcdfout.c
 *
 *    Scott Havens 07/2015
 *
 *    Write the data to a netcdf file
 */

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <time.h>

/* This is the name of the data file we will create. */
#define DK_TITLE "Created with DK tool"
#define VAR_NAME "variable"

/* We are writing 3D data, a 6 x 12 grid. */
#define NDIMS 3

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int netcdfout(iy, ip, data, nx, ny)
int iy;                          /* year */
int ip;                          /* period (sequential number beginning Oct 1) */
float *data;						 /* pointer to data matrix */
int nx;								/* number of values in x index */
int ny;								/* number of values in y index */
{

	/* When we create netCDF variables and dimensions, we get back an
	 * ID for each one. */
	int ncid, t_dimid, x_dimid, y_dimid, varid;
	int dimids[NDIMS];

	/* indexing and error handling. */
	int i, j, temp, retval;

	/* File Name */
	char file_name[20];
	snprintf(file_name, sizeof(file_name), "dk_out_%i.nc", iy);

	/* location to put the data */
	size_t start[] = {ip, 0, 0};
	size_t count[] = {1, nx, ny};


	/* data size */
//	printf("Year: %/i period: %i filename: %s\n", iy, ip, file_name);

	/* If there isn't a netcdf file, then create one, else just open it
	 * Create the file. NC_NOCOBBER will not overwrite an existing file
	 */
	retval = nc_create(file_name, NC_NOCLOBBER, &ncid);
	if ((retval == NC_EEXIST))
	{
		/* File exists so open it instead */
		if ((retval = nc_open(file_name, NC_WRITE, &ncid)))
			ERR(retval);

		/* Get the variable id */
		if ((retval = nc_inq_varid (ncid, VAR_NAME, &varid)))
			ERR(retval);

	}
	else if ((retval == NC_NOERR))
	{
		/* Define the dimensions. NetCDF will hand back an ID for each. */
		if ((retval = nc_def_dim(ncid, "Time", NC_UNLIMITED, &t_dimid)))
			ERR(retval);
		if ((retval = nc_def_dim(ncid, "y", ny, &y_dimid)))
					ERR(retval);
		if ((retval = nc_def_dim(ncid, "x", nx, &x_dimid)))
			ERR(retval);

		/* The dimids array is used to pass the IDs of the dimensions of
		 * the variable. */
		dimids[0] = t_dimid;
		dimids[1] = y_dimid;
		dimids[2] = x_dimid;

		/* Define variables */
		if ((retval = nc_def_var(ncid, VAR_NAME, NC_DOUBLE, NDIMS, dimids, &varid)))
			ERR(retval);

		/* Add global attributes for file */
		if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "Water_Year", NC_INT, 1, &iy)))
			ERR(retval);
		if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "Title", strlen(DK_TITLE), DK_TITLE)))
			ERR(retval);

		char currentDate[100];
		time_t now = time(NULL);
		struct tm *t = localtime(&now);
		strftime(currentDate, sizeof(currentDate)-1, "%d/%m/%Y %H:%M", t);
		if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "Created", strlen(currentDate), currentDate)))
			ERR(retval);

		retval = nc_enddef(ncid);

	}
	else
		ERR(retval);

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

	if ((retval = nc_put_vara_float(ncid, varid, start, count, data)))
		ERR(retval);



	/* Close the file. This frees up any internal netCDF resources
	 * associated with the file, and flushes any buffers. */
	if ((retval = nc_close(ncid)))
		ERR(retval);

	return 0;
}
