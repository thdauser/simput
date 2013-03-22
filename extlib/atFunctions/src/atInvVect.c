#include "atFunctions.h"
#include "atError.h"

int                     /* Making an inverse vector        */
atInvVect(              /*          ver 1.0  92/07/01      */
        AtVect x,       /* in: vector */
        AtVect y)       /* out: reversed vector */
{
	y[0] = -x[0];
	y[1] = -x[1];
	y[2] = -x[2];
	return NORMAL_END;
}
