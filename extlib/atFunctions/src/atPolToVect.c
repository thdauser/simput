#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * Converting the coordinate from Polar to Cartesian (POLREC)
 *                           ver 1.0  92/07/01  ay
 */
int
atPolToVect(
	AtPolarVect *x,	/* input */
	AtVect y		/* output: result */
)
{
	y[0] = x->r * cos( x->lat ) * cos( x->lon );
	y[1] = x->r * cos( x->lat ) * sin( x->lon );
	y[2] = x->r * sin( x->lat );
	return NORMAL_END;
}
