#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * Set Up geodetic vector to the tracking station and Rotation Matirix
 * to convert sidefloat coord to Ground (Earth-bound) coord
 *       Based on "ADAZEL"
 *
 * 2005/08/06 v2.5	Y.ISHISAKI
 *		modified to use atGeographic()
 *		remove static for local variables
 */
int
atAzElSet(
	AtPolarVect *vpSt,	/* input: geographic tracking station position :
							vpSt.r: altitude from Earth surface
							vpSt.lon, vpSt.lat: longitude and lattitude */
	AtVect vStation,	/* output: geodetic coordinate of station */
	AtRotMat stationRM) /* output: rotation matrix which converts geographic
							coord of tracking station to local coord */
{
	AtVect np, z;
	AtPolarVect north, zenith;

	atGeographic(vpSt, vStation);

/* set up Rotation matrix  for geodetic to ground coord */
	zenith.r = 1.0;
	zenith.lon = vpSt->lon;
	zenith.lat = vpSt->lat;
	atPolToVect(&zenith, z);

	north.r = 1.0;
	north.lon = vpSt->lon + PI;
	north.lat = PI/2 - vpSt->lat;
	atPolToVect(&north, np);

	atSetRotMatZX(z, np, stationRM);

	return 0;
}
