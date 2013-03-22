#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert geographic position on earth surface into geodetic vector
 */
int
atGeodesic(AtPolarVect *pv, AtVect vect)
{
	return atGeographic(pv, vect);
}
