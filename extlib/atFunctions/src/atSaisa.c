#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * conversion of equatorial coordinate systems
 * correcting for precession (based on SAISA in ASTRO-C csub1)
 * Ver 1.1 (Dec 6 1992)
 */

int atSaisa(
	double mjd0,		/* input: MJD for the original coord*/
	AtPolarVect *pv0,	/* input: original polar coordinate */
	double mjd,		/* input: MJD for the new coord*/
	AtPolarVect *pv)	/* output: Euler Ang to correct precession*/
{
    double cosx3, sinx3, u, x, y, z, u0, x1, x2, x3, cosal0,
	    cosdl0, sindl0;

    u0 = (mjd0 - 15020.) / 36524.22;
    u = (mjd - mjd0) / 36524.22;
    x1 = (u0 * 1.396 + 2304.25 + (u * .018 + .302) * u) * u * DEG2RAD / 3600.;
    x2 = x1 + u * .791 * u * DEG2RAD / 3600.;
    x3 = (2004.682 - u0 * .853 + (-.426 - u * .042) * u) * u * DEG2RAD / 3600.;
    cosal0 = cos(pv0->lon + x1);
    sindl0 = sin(pv0->lat);
    cosdl0 = cos(pv0->lat);
    cosx3 = cos(x3);
    sinx3 = sin(x3);
    x = cosx3 * cosdl0 * cosal0 - sinx3 * sindl0;
    y = cosdl0 * sin(pv0->lon + x1);
    z = cosx3 * sindl0 + sinx3 * cosdl0 * cosal0;
    pv->lon = atan2(y, x);
    x /= cos(pv->lon);
    pv->lat = atan2(z, x);
    pv->lon += x2;
    if (pv->lon < 0.) pv->lon += TWO_PI;
    return 0;
}
