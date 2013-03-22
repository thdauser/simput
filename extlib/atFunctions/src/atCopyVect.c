#include "atFunctions.h"
#include "atError.h"

/*
 * copy vector (z = x)
 *             			ver 1.5  93/01/26  n.kawai
 */
int atCopyVect(
	AtVect x,	/* in: vector */
	AtVect z)	/* out: copy vector of x*/
{
    z[0] = x[0];
    z[1] = x[1];
    z[2] = x[2];
    return (NORMAL_END);
}
