/*----------------------------------------------------------------------------
*
* tspctrne does a quick test of spctrne().  Not part of the official test
* suite.
*
* $Id: tspctrne.c,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include <spc.h>
#include <wcserr.h>

int main()

{
  const char ctypeS1[] = "VOPT-F2W";
  const double crvalS1 = 1e6;
  const double cdeltS1 = 1e3;
  const double restfrq = 0.0;
  const double restwav = 0.0;

  int    status;
  char   ctypeS2[9];
  double cdeltS2, crvalS2;
  struct wcserr *err;

  strcpy(ctypeS2, "VRAD-???");

  wcserr_enable(1);
  if (spctrne(ctypeS1, crvalS1, cdeltS1, restfrq, restwav,
              ctypeS2, &crvalS2, &cdeltS2, &err)) {
    wcserr_prt(err, 0x0);
    return err->status;
  }

  printf("'%8s'  %12.6e  %12.6e\n'%8s'  %12.6e  %12.6e\n",
    ctypeS1, crvalS1, cdeltS1, ctypeS2, crvalS2, cdeltS2);

  return 0;
}
