#include "atFunctions.h"

/*
 * Azimuth and Elevation for a vector in Ground (Earth-bound) coord
 */
int
atAzEl(
	AtVect x,	    	/* input: satellite pos. in geodetic coord. */
	AtVect vStation,	/* input: station pos. in geodetic coord. */
	AtRotMat stationRM, /* input: rotation matrix which converts geographic
							coord of tracking station to local coord */
	AtPolarVect *y     /* output: local coordinate:
			    			y.r: distance, y.lon:azimuth, y.lat:elevation */
)
{
    int istat;
    AtVect xLocal, xx;

	xx[0] = x[0] - vStation[0];
	xx[1] = x[1] - vStation[1];
	xx[2] = x[2] - vStation[2];
    ATRotVect(stationRM, xx, xLocal);
    istat = atVectToPol(xLocal, y);
    y->lon = TWO_PI - y->lon;

    return istat;
}
