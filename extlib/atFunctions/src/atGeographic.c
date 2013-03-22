#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert geographic location of Tracking Station on Earth to Geodetic
 */
int
atGeographic(
	AtPolarVect *y,	/* input: geographic tracking station position :
							y.r: altitude from Earth surface
							y.lon, y.lat: longitude and lattitude */
	AtVect z		/* output: geodetic coordinate of station */
)
{
	double sinlat, rc, rq;

	sinlat = sin(y->lat);
	rc = EARTH_RADIUS / sqrt(1. - EARTH_E2 * sinlat * sinlat);

	rq = (rc + y->r) * cos(y->lat);
	z[0] = rq * cos(y->lon);
	z[1] = rq * sin(y->lon);
	z[2] = (rc * (1. - EARTH_E2) +  y->r) * sinlat;

	return NORMAL_END;
}
