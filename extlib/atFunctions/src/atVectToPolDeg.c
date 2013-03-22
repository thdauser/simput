#include "atFunctions.h"

/*
 * covert Vector to R.A. and Dec in degree.
 */
int atVectToPolDeg(
	    AtVect x,		/* input: vector */
	    double *r, 		/* vector length */
	    double *alpha, 	/* R.A. in degree */
	    double *delta) 	/* Dec. in degree */
{
  AtPolarVect  y;
  int code;

  code = atVectToPol(x, &y);
  *r = y.r ;
  *alpha = y.lon * RAD2DEG;
  *delta = y.lat * RAD2DEG;
  return code;
}
