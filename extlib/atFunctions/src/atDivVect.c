#include "atFunctions.h"
#include "atError.h"

/*
 * divide scaler for vector (z = (1/d)*x)
 */
int
atDivVect(
	double d,	/* in: divisor for x */
	AtVect x,	/* in: vector */
	AtVect z)	/* out: sum vector of x and y */
{
	if ( 0.0 == d ) {
		return NULL_VECTOR;
	}
    z[0] = x[0] / d;
    z[1] = x[1] / d;
    z[2] = x[2] / d;
    return NORMAL_END;
}
