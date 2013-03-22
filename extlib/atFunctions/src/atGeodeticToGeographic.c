#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert geodetic vector into geographic position on earth surface
 */
int
atGeodeticToGeographic(AtVect vect, AtPolarVect *pv)
{
	int istat;
	double latt, heigh;

	atVectToPol(vect, pv);
	istat = atEllipsoid(pv, &latt, &heigh);
	pv->r = heigh;
	pv->lat = latt;

	return istat;
}
