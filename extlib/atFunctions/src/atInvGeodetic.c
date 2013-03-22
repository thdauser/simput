#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert geodetic corrdinate to sidereal equatorial coord
 */
int
atInvGeodetic(
	double mjd,	/* input: time in MJD */
	AtVect x,	/* intput: vector in geodetic coordinate at mjd,
					i.e. x-axis at (long, lat)=(0,0), and z at N-pole*/
	AtVect y)	/* output: vector in sidereal equatorial coordinate */
{
    int code1;
    AtRotMat rm, inv_rm;

    code1 = atSetGeoRM(mjd, rm);
	ATInvRotMat(rm, inv_rm);
    ATRotVect(inv_rm, x, y);

    return code1;
}
