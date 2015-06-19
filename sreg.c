/*
 *    sreg.c
 *
 *    David Garen  11/89, 8/91, 11/91, 8/92
 *
 *    Perform simple linear regression.  Computes slope, intercept,
 *    correlation coefficient, standard error, and t-statistic.
 *
 *    Returns 0 for normal completion,
 *            1 if x data are all equal,
 *            3 if both x and y data are all equal.
 */

#include <math.h>

int sreg(x, y, b0, b1, r, se, t, n)
double *x;                       /* independent variable data vector */
double *y;                       /* dependent variable data vector */
double *b0;                      /* intercept */
double *b1;                      /* slope */
double *r;                       /* correlation coefficient */
double *se;                      /* standard error */
double *t;                       /* t-statistic */
int n;                           /* number of observations */
{
   double dum1, dum2, dum3;      /* dummy variables for calculations */
   double err;                   /* error */
   int i;                        /* loop index */
   double sqrt();                /* square root function */
   double sumx = 0;              /* sum of x variable */
   double sumx2 = 0;             /* sum of x variable squared */
   double sumxy = 0;             /* sum of x times y */
   double sumy = 0;              /* sum of y variable */
   double sumy2 = 0;             /* sum of y variable squared */
   int xequal = 0;               /* flag for all x data being equal */
   int yequal = 0;               /* flag for all y data being equal */

   /* Check for all x or y data being equal */

   for (i = 1; i < n; i++)
      if (x[i] != x[0])
         break;
   if (i == n)
      xequal = 1;

   for (i = 1; i < n; i++)
      if (y[i] != y[0])
         break;
   if (i == n)
      yequal = 1;

   if (xequal == 1 && yequal != 1)
      return(1);
   if (xequal == 1 && yequal == 1)
      return(3);

   /* Do regression */

   for (i = 0; i < n; i++) {
      sumx += x[i];
      sumx2 += x[i] * x[i];
      sumxy += x[i] * y[i];
      sumy += y[i];
      sumy2 += y[i] * y[i];
   }
   dum1 = n * sumxy - sumx * sumy;
   dum2 = n * sumx2 - sumx * sumx;
   dum3 = n * sumy2 - sumy * sumy;
   *b1 = dum1 / dum2;
   *b0 = (sumy / n) - (*b1 * sumx / n);
   if (yequal == 1) {
      *r = 1;
      *se = 0;
      *t = 0;
      return(0);
   }
   *r = dum1 / sqrt(dum2 * dum3);
   if (n > 2) {
      dum1 = 0;
      for (i = 0; i < n; i++) {
         err = y[i] - *b0 - *b1 * x[i];
         dum1 += (err * err);
      }
      dum1 /= (n - 2);
      *se = sqrt(dum1);
      *t = (*b1 / *se) * sqrt(dum2 / n);
   }
   else {
      *se = 0;
      *t = 0;
   }
   return(0);
}

