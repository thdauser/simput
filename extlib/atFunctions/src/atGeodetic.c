#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert sidereal equatorial coord to geodetic corrdinate
 */
int
atGeodetic(
	double mjd,	/* input: time in MJD */
	AtVect x,	/* input: vector in sidereal equatorial coordinate */
	AtVect y)	/* output: vector in geodetic coordinate at mjd ,
					i.e. x-axis at (long, lat)=(0,0), and z at N-pole*/
{
    int code1;
    AtRotMat rm;

    code1 = atSetGeoRM(mjd, rm);
    ATRotVect(rm, x, y);

    return code1;
}
