/*
 *    isleap.c
 *
 *    David Garen   9/91
 *
 *    Determine if year is a leap year (1 = leap year; 0 = not leap year)
 */

int isleap(year)
int year;
{
   int dum;

   dum = year / 4;
   dum *= 4;
   if ((dum = year - dum) == 0)
      return(1);
   else
      return(0);
}
