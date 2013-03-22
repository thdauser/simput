/************************************************************************
  atAngDistance.c	calc angular distance (RIKAKU)

	1992/07/01 A.YOSHIDA	version 1.0

	1993/01/26 N.KAWAI		version 1.5

	2005/12/04 Y.ISHISAKI	version 2.7
		use fast macros
************************************************************************/
#include <stdio.h>
#include <math.h>
#include "atFunctions.h"
#include "atError.h"

/*
 * Normalizing a vector
 */
int
atNormVect(
	   AtVect x,       /* input:  vector */
	   AtVect y)       /* output: normalized vector*/
{
	double nrm;

	nrm = ATScalProd(x, x);

	if ( 0.0 == nrm ) {
		y[0] = y[1] = y[2] = 0.0;
		return NULL_VECTOR;

    } else if ( 1.0 == nrm ) {
		ATCopyVect(x, y);

    } else {
		nrm = sqrt(nrm);
		y[0] = x[0] / nrm;
		y[1] = x[1] / nrm;
		y[2] = x[2] / nrm;
    }

	return NORMAL_END;
}
