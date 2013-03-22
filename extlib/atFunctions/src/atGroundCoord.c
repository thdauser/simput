#include "atFunctions.h"

/*
 * Convert from geocentric to Local Coordinate at the Tracking Station
 */
int atGroundCoord(
	AtVect station,	/* input: tracking station in geographic coord */
	AtRotMat rm)		/* output: rotation matrix to local coord*/

{
    static AtVect north = {0., 0., 1.};
    int code;

    code = atSetRotMatZX(station, north, rm);
    return code;
}
