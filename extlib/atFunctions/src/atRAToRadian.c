#include "atFunctions.h"

double                  /* Converting R.A. in hh.mm.ss. to radian */
atRAToRadian(                  /* out: right ascension in radian  */
        AtRightAscension ra)   /* in:  right ascension in hhmmss  */
{                              /*           ver 1.0  92/07/01  ay */
  double  x;
  x = ra.hour*15.0 + ra.min*0.25 + ra.sec*0.25/60.0;
  return (x * DEG2RAD);
}
