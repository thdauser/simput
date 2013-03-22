#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert equatorial coordinate systems correcting for precession (SAISA)
 */
int atPrecession(
	double mjd0,	/* input: modified Julian Day */
	AtVect x0,	/* input: vector in equatorial coordiante at mjd0 */
	double mjd,	/* input: modified Julian Day */
	AtVect x)	/* output: vector in equatorial coordiante at mjd1 */
{
    AtRotMat rm;

    atPrecessRM(mjd0, mjd, rm);
    ATRotVect(rm, x0, x);
    return 0;
}
