/*----------------------------------------------------------------------------
*
* tspcaips does a quick test of spcaips().  Not part of the official test
* suite.
*
* $Id: tspcaips.c,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <spc.h>

int main()

{
  const char *(ctypes[]) = {"FREQ", "VELO", "FELO"};
  const char *(frames[]) = {"-LSR", "-HEL", "-OBS", "    "};

  char ctype[9], ctypeA[9], specsys[9];
  int  i, j, status, v1, v2, velref;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++) {
      sprintf(ctypeA, "%s%s", ctypes[i], frames[j]);

      for (v1 = 0; v1 <= 8; v1++) {
        velref = v1;
        for (v2 = 0; v2 < 3; v2++) {
          status = spcaips(ctypeA, velref, ctype, specsys);
          printf("'%s'  %3d  %2d  '%s'  '%s'\n", ctypeA, velref,
            status, ctype, specsys);
          velref += 256;
        }
      }

      printf("\n");
    }
  }

  return 0;
}
