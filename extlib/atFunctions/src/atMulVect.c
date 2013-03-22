#include "atFunctions.h"
#include "atError.h"

/*
 * multiply scaler for vector (z = f*x)
 */
int
atMulVect(
	double f,	/* in: multiplicand for x */
	AtVect x,	/* in: vector */
	AtVect z)	/* out: sum vector of x and y */
{
    z[0] = f * x[0];
    z[1] = f * x[1];
    z[2] = f * x[2];
    return NORMAL_END;
}
