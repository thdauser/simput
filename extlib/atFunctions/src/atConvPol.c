#include "atFunctions.h"
#include "atError.h"
#include <math.h>

int             /* Converting R.A. & Dec. from hh.mmss. to radian */
atConvPol(      /*                    ver 1.0  92/07/01    ay     */
        AtPolarVect60 *x,       /* input */
        AtPolarVect   *y)       /* result */
{
  y->lon = atRAToRadian( x->ra );
  y->lat = atDecToRadian( x->dec );
  return (NORMAL_END);
}
