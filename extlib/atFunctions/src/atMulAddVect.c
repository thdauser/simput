#include "atFunctions.h"
#include "atError.h"

/*
 * multiply and add vectors (z = f*x + g*y)
 *			           		 ver 1.5  93/01/25  n.kawai
 */
int atMulAddVect(
	double f,	/* in: multiplicand for x */
	AtVect x,	/* in: vector */
	double g,	/* in: multiplicand for y */
	AtVect y,	/* in: vector */
	AtVect z)	/* out: answer*/
{
    z[0] = f*x[0] + g*y[0];
    z[1] = f*x[1] + g*y[1];
    z[2] = f*x[2] + g*y[2];
    return (NORMAL_END);
}
