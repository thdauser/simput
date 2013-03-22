#include "atFunctions.h"
#include "atError.h"

/*
 * subtract vectors (z = x - y)
 */
int
atSubVect(
	AtVect x,	/* in: vector */
	AtVect y,	/* in: vector */
	AtVect z)	/* out: vector of (x - y) */
{
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
    z[2] = x[2] - y[2];
    return NORMAL_END;
}
