#include <stdio.h>

/*
 *    getln.c
 *
 *    David Garen  3/89
 *
 *    Retrieve a line from a file.
 *
 *    Return 0 for normal operation or
 *           EOF if end of file encountered.
 */

int getln(sp, fp)
char *sp;                        /* string pointer */
FILE *fp;                        /* file pointer */
{
   int ch;
   int ch1;

   while((ch = fgetc(fp)) != '\n' && ch != EOF && ch != '\r')
      *sp++ = (char)ch;

   ch1 = fgetc(fp);
   if (ch1 != '\n' && ch1 != '\r')
       ungetc(ch1, fp);

   *sp = '\0';
   return((ch == EOF) ? EOF : 0);
}

