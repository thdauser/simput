#include "atFunctions.h"

/*
 * covert R.A. and Dec in degree to Vector
 */
int
atPolDegToVect(
	    double r, 		/* vector length */
	    double alpha, 	/* R.A. in degree */
	    double delta, 	/* Dec. in degree */
	    AtVect x)
{
	AtPolarVect  y;
	int code;

	y.lon = alpha * DEG2RAD;
	y.lat = delta * DEG2RAD;
	y.r = r;
	code = atPolToVect(&y, x);

	return code;
}
