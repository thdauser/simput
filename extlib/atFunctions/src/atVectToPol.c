#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/* Converting the coordinate from Cartesian to Polar (RECPOL)
	1992/07/01	ver 1.0	ay
	2005/08/06	ver 2.5	Y.ISHISAKI
		modified to use atan2(), which is safer & faster
*/
int
atVectToPol(
	AtVect x,			/* input */
	AtPolarVect *y		/* output: result */
)
{
	double  norm01;

	norm01 = x[0]*x[0] + x[1]*x[1];
	y->r = sqrt( norm01 + x[2]*x[2] );
	if ( 0.0 == y->r ) {
		y->lon = y->lat = 0.0;
		return NULL_VECTOR;
    }

	norm01 = sqrt( norm01 );
	y->lat = atan2(x[2], norm01);	/* -2/PI <= atan2 <= 2/PI, as 0<=norm01 */
	y->lon = atan2(x[1], x[0]);		/* returns -PI <= atan2 <= PI */
	if ( y->lon < 0.0 ) {
		y->lon += TWO_PI;
	}

	return NORMAL_END;
}
