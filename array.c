#include <stdio.h>
#include <stdlib.h>
#include "malloc.h"

/*
 *    vector.c
 *
 *    David Garen  8/89
 *
 *    Allocate a float vector with n elements.
 *
 *    This program is a modified version of one from:
 *    William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 705.
 */

float *vector(n)
int n;
{
   float *v = (float *) NULL;

   if (n > 0)
      v = (float *) malloc(n * sizeof(float));
   if (!v) {
      printf("\n\nAllocation failure in vector().\n");
      exit(0);
   }
   return(v);
}

/*
 *    dvector.c
 *
 *    David Garen  8/89
 *
 *    Allocate a double vector with n elements.
 *
 *    This program is a modified version of one from:
 *    William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 706.
 */

double *dvector(n)
int n;
{
   double *v = (double *) NULL;

   if (n > 0)
      v = (double *) malloc(n * sizeof(double));
   if (!v) {
      printf("\n\nAllocation failure in dvector().\n");
      exit(0);
   }
   return(v);
}

/*
 *    ivector.c
 *
 *    David Garen  8/89
 *
 *    Allocate an integer vector with n elements.
 *
 *    This program is a modified version of one from:
 *    William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 706.
 */

int *ivector(n)
int n;
{
   int *v = (int *) NULL;

   if (n > 0)
      v = (int *) malloc(n * sizeof(int));
   if (!v) {
      printf("\n\nAllocation failure in ivector().\n");
      exit(0);
   }
   return(v);
}

/*
 *    matrix.c
 *
 *    David Garen  8/89
 *
 *    Allocate a float matrix with nr rows and nc columns.
 *
 *    This program is a modified version of one from:
 *    William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 706.
 */

float **matrix(nr, nc)
int nr, nc;
{
   int i;
   float **m = (float **) NULL;

   /* Allocate pointers to rows */

   if (nr > 0)
      m = (float **) malloc(nr * sizeof(float *));
   if (!m) {
      printf("\n\nAllocation failure 1 in matrix()\n");
      exit(0);
   }

   /* Allocate rows and set pointers to them */

   for (i = 0; i < nr; i++) {
      m[i] = (float *) NULL;
      if (nc > 0)
         m[i] = (float *) malloc(nc * sizeof(float));
      if (!m[i]) {
         printf("\n\nAllocation failure 2 in matrix()\n");
         exit(0);
      }
   }

   /* Return pointer to array of pointers to rows */

   return(m);
}

/*
 *    dmatrix.c
 *
 *    David Garen  8/89
 *
 *    Allocate a double matrix with nr rows and nc columns.
 *
 *    This program is a modified version of one from:
 *    William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 706.
 */

double **dmatrix(nr, nc)
int nr, nc;
{
   int i;
   double **m = (double **) NULL;

   /* Allocate pointers to rows */

   if (nr > 0)
      m = (double **) malloc(nr * sizeof(double *));
   if (!m) {
      printf("\n\nAllocation failure 1 in dmatrix()\n");
      exit(0);
   }

   /* Allocate rows and set pointers to them */

   for (i = 0; i < nr; i++) {
      m[i] = (double *) NULL;
      if (nc > 0)
         m[i] = (double *) malloc(nc * sizeof(double));
      if (!m[i]) {
         printf("\n\nAllocation failure 2 in dmatrix()\n");
         exit(0);
      }
   }

   /* Return pointer to array of pointers to rows */

   return(m);
}

/*
 *    imatrix.c
 *
 *    David Garen  8/89
 *
 *    Allocate an integer matrix with nr rows and nc columns.
 *
 *    This program is a modified version of one from:
 *    William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 707.
 */

int **imatrix(nr, nc)
int nr, nc;
{
   int i, **m = (int **) NULL;

   /* Allocate pointers to rows */

   if (nr > 0)
      m = (int **) malloc(nr * sizeof(int *));
   if (!m) {
      printf("\n\nAllocation failure 1 in imatrix()\n");
      exit(0);
   }

   /* Allocate rows and set pointers to them */

   for (i = 0; i < nr; i++) {
      m[i] = (int *) NULL;
      if (nc > 0)
         m[i] = (int *) malloc(nc * sizeof(int));
      if (!m[i]) {
         printf("\n\nAllocation failure 2 in imatrix()\n");
         exit(0);
      }
   }

   /* Return pointer to array of pointers to rows */

   return(m);
}

/*
 *    cube.c
 *
 *    David Garen  12/93
 *
 *    Allocate a float 3-D matrix with nr rows, nc columns,
 *    and nd elements deep.
 */

float ***cube(nr, nc, nd)
int nr, nc, nd;
{
   int i, j;
   float ***m = (float ***) NULL;

   if (nr > 0)
      m = (float ***) malloc(nr * sizeof(float **));
   if (!m) {
      printf("\n\nAllocation failure 1 in cube()\n");
      exit(0);
   }

   for (i = 0; i < nr; i++) {
      m[i] = (float **) NULL;
      if (nc > 0)
         m[i] = (float **) malloc(nc * sizeof(float *));
      if (!m[i]) {
         printf("\n\nAllocation failure 2 in cube()\n");
         exit(0);
      }
   }

   for (i = 0; i < nr; i++) {
      for (j = 0; j < nc; j++) {
         m[i][j] = (float *) NULL;
         if (nd > 0)
            m[i][j] = (float *) malloc(nd * sizeof(float));
         if (!m[i][j]) {
            printf("\n\nAllocation failure 3 in cube()\n");
            exit(0);
         }
      }
   }

   return(m);
}

