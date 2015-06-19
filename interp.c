/*
 *    interp.c
 *
 *    David Garen   9/92
 *
 *    Calculate interpolated value for accumulated daily precipitation
 */

#include <stdio.h>

#include "dk_x.h"

float interp(i, jstart, kstart)
int i;                           /* station index */
int jstart;                      /* day index of start of accum. period */
int kstart;                      /* year index of start of accum. period */
{
   int j, k;                     /* loop indexes */
   int n = 1;                    /* day counter */

   dstop = lastday[kstart] - 1;
   j = jstart;
   k = kstart;
   while (1) {
      if (sta[i].data[j][k] > accum && sta[i].data[j][k] < missing) {
         n++;
         j++;
         if (j > dstop) {
            j = 0;
            k++;
         }
      }
      else
         break;
   }
   dum = sta[i].data[j][k] / n;
   sta[i].data[j][k] = (float) (replace + 0.1);
   return(dum);
}

