/*
 *    medfit.c
 *
 *    David Garen  9/95
 *
 *    Fits y = a + bx by the criterion of least absolute deviations.  This
 *    is a robust regression procedure from:
 *
 *    W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
 *    Numerical Recipes in C, Cambridge University Press, 1988, pp. 563, 564.
 */

#include <stdlib.h>
#include <malloc/malloc.h>
#include <math.h>
/* Debug
#include <stdio.h>
   End debug */

#define ABS(a) ((a) >= 0 ? (a) : -(a))      /* absolute value operator */

int nt = 0;
double *xt = 0, *yt = 0, aa = 0, abdevt = 0;
/* FILE *debug;
   FILE *fopen(); */

int medfit(x, y, a, b, abdev, n)
double *x;                       /* independent variable data vector */
double *y;                       /* dependent variable data vector */
double *a;                       /* intercept */
double *b;                       /* slope */
double *abdev;                   /* average absolute deviation */
int n;                           /* number of observations */
{
   int j;
   double bb, b1, b2, del, f, f1, f2, sigb, temp;
   double sx = 0, sy = 0, sxy = 0, sxx = 0, chisq = 0;
   double rofunc();
   double sqrt();
   int xequal = 0;               /* flag for all x data being equal */
   int yequal = 0;               /* flag for all y data being equal */

/* debug = fopen("debug.out", "a"); */

   /* Check for all x or y data being equal */

   for (j = 1; j < n; j++)
      if (x[j] != x[0])
         break;
   if (j == n)
      xequal = 1;

   for (j = 1; j < n; j++)
      if (y[j] != y[0])
         break;
   if (j == n)
      yequal = 1;

   if (xequal == 1 && yequal != 1)
      return(1);
   if (xequal == 1 && yequal == 1)
      return(3);

   nt = n;
   xt = x;
   yt = y;

   /* As a first guess for a and b, find the least squares regression line */

   for (j = 0; j < n; j++) {
      sx += x[j];
      sy += y[j];
      sxy += x[j] * y[j];
      sxx += x[j] * x[j];
   }
   del = n * sxx - sx * sx;
   aa = (sxx * sy - sx * sxy) / del;
   bb = (n * sxy - sx * sy) / del;

   /* The standard deviation of the slope will give some idea of how big an
      iteration step to take */

   for (j = 0; j < n; j++) {
      temp = y[j] - (aa + bb * x[j]);
      chisq += (temp * temp);
   }
   sigb = sqrt(chisq / del);
/* Debug
fprintf(debug, "\nmedfit(): aa = %f  bb = %f  chisq = %f  sigb = %f\n",
       aa, bb, chisq, sigb);
   End debug */

   /* Guess bracket as 3-sigb away, in the downhill direction known from f1 */

   b1 = bb;
   f1 = rofunc(b1);
/* Debug
fprintf(debug, "\nmedfit():  b1 = %f   f1 = rofunc(b1) = %f   abdevt = %f\n",
       b1, f1, abdevt);
   End debug */
   b2 = bb + ((f1 > 0.0) ? ABS(3.0 * sigb) : -ABS(3.0 * sigb));
   f2 = rofunc(b2);
/* Debug
fprintf(debug, "\nmedfit():  b2 = %f   f2 = rofunc(b2) = %f   abdevt = %f\n",
       b2, f2, abdevt);
   End debug */

   /* Bracketing */

   while ((f1 * f2) > 0.0) {
      bb = 2.0 * b2 - b1;
      b1 = b2;
      f1 = f2;
      b2 = bb;
      f2 = rofunc(b2);
   }
/* Debug
fprintf(debug, "\nmedfit(): Results after bracketing:\nb1 = %f  f1 = %f  b2 = %f  f2 = %f  abdevt = %f\n",
       b1, f1, b2, f2, abdevt);
   End debug */

   /* Refine until error is a negligible number of standard deviations */

   sigb = 0.01 * sigb;
   while (ABS(b2 - b1) > sigb) {
      bb = 0.5 * (b1 + b2);
      if (bb == b1 || bb == b2)
         break;
      f = rofunc(bb);
      if ((f * f1) >= 0.0) {
         f1 = f;
         b1 = bb;
      }
      else {
         f2 = f;
         b2 = bb;
      }
   }
   *a = aa;
   *b = bb;
   *abdev = abdevt / n;
/* Debug
fprintf(debug, "\nmedfit(): Final results:  a = %f   b = %f   abdev = %f\n",
       *a, *b, *abdev);
fclose(debug);
   End debug */
   return(0);
}

/*
 *    rofunc.c
 *
 *    David Garen  9/95
 *
 *    Intermediate function evaluation for least absolute deviation
 *    regression fit.  From:
 *
 *    W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
 *    Numerical Recipes in C, Cambridge University Press, 1988, p 564.
 */

double rofunc(b)
double b;
{
   int j, n1, nmh, nml;
   double *arr, d, sum = 0;
   double *dvector();
/* Test of indexx() */
   int *inx;
   int *ivector();
   void indexx();
/* void sort(); */

/* Debug
fprintf(debug, "\nbeginning rofunc():  b = %f\n", b);
   End debug */
   arr = dvector(nt);
   inx = ivector(nt);
   n1 = nt + 1;
   nml = n1 / 2 - 1;
   nmh = n1 - nml - 2;
   for (j = 0; j < nt; j++)
      arr[j] = yt[j] - b * xt[j];
/* Debug
fprintf(debug, "\nrofunc():  arr array:\n");
for (j = 0; j < nt; j++)
   fprintf(debug, "%f\n", arr[j]);
   End debug */
   indexx(arr, inx, nt);
   aa = 0.5 * (arr[inx[nml]] + arr[inx[nmh]]);
/* Debug
fprintf(debug, "\nrofunc():  aa = %f\n", aa);
   End debug */
/* sort(nt, arr); */
/* aa = 0.5 * (arr[nml] + arr[nmh]); */
   abdevt = 0;
/* Debug
fprintf(debug, "\nrofunc():  d values:\n");
   End debug */
   for (j = 0; j < nt; j++) {
      d = yt[j] - (aa + b * xt[j]);
/* Debug
fprintf(debug, "%f\n", d);
   End debug */
      abdevt += ABS(d);
      if (d > 0.0001)
         sum += xt[j];
      else if (d < -0.0001)
         sum -= xt[j];
/*    sum += (d > 0.0 ? xt[j] : -xt[j]); */
   }
/* Debug
fprintf(debug, "\nfinishing rofunc():  abdevt = %f   sum = %f\n", abdevt, sum);
   End debug */
   free(arr);
   free(inx);
   return(sum);
}

/*
 *    sort.c
 *
 *    David Garen  9/95
 *
 *    Sorts an array into ascending numerical order using the Heapsort
 *    algorithm.  From:
 *
 *    W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
 *    Numerical Recipes in C, Cambridge University Press, 1988, p 247.
 */

void sort(n, ra)
int n;
double *ra;
{
   int i, ir, j, l;
   double rra;

   l = n / 2;
   ir = n - 1;

   /* The index l will be decremented from its initial value down to 0 during
      the "hiring" (heap creation) phase.  Once it reaches 0, the index ir
      will be decremented from its initial value down to 0 during the
      "retirement-and-promotion" (heap selection) phase. */

   while (1) {
      if (l > 0) {

         /* Still in "hiring" phase */

         rra = ra[--l];
      }
      else {

         /* In "retirement-and-promotion" phase */

         rra = ra[ir];
         ra[ir] = ra[0];
         if (--ir == 0) {
            ra[0] = rra;
            return;
         }
      }

      /* Whether in the hiring or promotion phase, set up to sift down
         element rra to its proper level */

      i = l;
      j = l * 2 + 1;
      while (j <= ir) {
         if (j < ir && ra[j] < ra[j+1])
            ++j;
         if (rra < ra[j]) {
            ra[i] = ra[j];
            j += (i = j);
         }
         else
            j = ir + 1;
      }
      ra[i] = rra;
   }
}

