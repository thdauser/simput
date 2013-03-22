#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * examine earth occultation of specified direction (YEARTH)
 */
int atEarthOccult(
	AtVect satVect,	/* input: satellite pos. in equatorial coord. */
	AtVect xVect,	/* input: direction to be examined. (normalized) */
	AtVect sunVect,	/* input: vector to the sun */
	int *flag, 	/* output: condition 0:sky, 1:dark, 2:bright earth */
	double *el)	/* output: elevation angle of z from the earth edge */
{
    double earthSize, xDist, dot, zCross, satDistance;
    AtRotMat rm;
    AtVect satV, earthVect, xV;

    ATInvVect(satVect, earthVect);
    satDistance = ATNorm(earthVect);
    earthSize = asin(EARTH_RADIUS / satDistance);
    atAngDistance(xVect, earthVect, &xDist);
    *el = xDist - earthSize;

    if (*el <= - EPS) {
	atSetRotMatZX(sunVect, satVect, rm);
	ATRotVect(rm, satVect, satV);
	ATRotVect(rm, xVect, xV);
	dot = ATScalProd(earthVect, xVect);
	zCross = satV[2] + xV[2]* (dot
		- sqrt( EARTH_RADIUS * EARTH_RADIUS
		- satDistance * satDistance + dot * dot));
	if (zCross < 0.) {
	    *flag = 1;		/* Dark Earth */
	} else {
	    *flag = 2;		/* Bright Earth */
	}
    } else {
	*flag = 0;
    }
    return 0;
}
