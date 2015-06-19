/*
 *    index.c
 *
 *    David Garen  8/89
 *
 *    Indexes an array arrin[0..n-1], i.e., outputs the array indx[0..n-1]
 *    such that arrin[indx[j]] is in ascending order for j = 0, 1, ... , n-1.
 *    The input quantities n and arrin are not changed.
 *    
 *    This program is slightly modified from the program "indexx" in:
 *    Press, William H., Brian P. Flannery, Saul A. Teukolsky, and William T.
 *    Vetterling, Numerical Recipes in C:  The Art of Scientific Computing,
 *    Cambridge University Press, 1988, p. 249.
 */

void indexx(arrin, indx, n)
double *arrin;                   /* input vector */
int *indx;                       /* index vector */
int n;                           /* number of elements to sort */
{
   int i, indxt, ir, j, l;
   double q;

   /* Initialize the index vector with consecutive integers */

   for (j = 0; j < n; j++)
      indx[j] = j;

   /* From here on, we just have Heapsort, but with indirect indexing
      through indx in all references to arrin */

   l = n / 2 + 1;
   ir = n;
   while (1) {
      if (l > 1)
         q = arrin[(indxt = indx[--l-1])];
      else {
         q = arrin[(indxt = indx[ir-1])];
         indx[ir-1] = indx[0];
         if (--ir == 1) {
            indx[0] = indxt;
            return;
         }
      }
      i = l;
      j = l * 2;
      while (j <= ir) {
         if (j < ir && arrin[indx[j-1]] < arrin[indx[j]])
            j++;
         if (q < arrin[indx[j-1]]) {
            indx[i-1] = indx[j-1];
            j += (i = j);
         }
         else
            j = ir + 1;
      }
      indx[i-1] = indxt;
   }
}

