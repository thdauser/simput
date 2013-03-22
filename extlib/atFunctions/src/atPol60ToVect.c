#include "atFunctions.h"

int               /* Coverting R.A. and Dec in hr/deg,min,s to Vector. */
atPol60ToVect(
	    AtPolarVect60 *p,	/* polar vector in hh/deg mm ss */
	    AtVect x) 		/* output: vector */
{                               /*       ver 1.0 92/07/01   ay  */
  AtPolarVect   pv;
  int  rc;
  pv.lon = atRAToRadian( p->ra );
  pv.lat = atDecToRadian( p->dec );
  rc = atPolToVect( &pv, x );
  return ( rc );
}
