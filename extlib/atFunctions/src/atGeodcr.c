#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert sidereal equatorial coordinate
 * to the geographic position on earth surface(GEODCR)
 *              Ver 1.6 (93/04/07)  Loop logic revised by N.Kawai
 *
 *	2005/08/06	Ver 2.5	Y.ISHISAKI
 *		modified to use atGeodeticToGeographic()
 */
int
atGeodcr(
	double mjd,		/* input: time in MJD */
	AtVect x,		/* input: vector in celestial coordinate */
	double *heigh,	/* output: altitude from the earth surface */
	double *longi,	/* output: longitude on the earth surface */
	double *latt	/* output: latitude on the earth surface */
)
{
/* Local variables */
	double gsttod;
	AtPolarVect xp;
	AtVect vec;

/* convert R.A. -> Longitude */
	atVectToPol(x, &xp);
	atSidereal(mjd, &gsttod);
	xp.lon = xp.lon - gsttod;

/* RS,DS -) HEIGH, LATT  (DISTORTIION OF THE EARTH) */
	atPolToVect(&xp, vec);
	atGeodeticToGeographic(vec, &xp);
	*heigh = xp.r;
	*longi = xp.lon;
	*latt  = xp.lat;

	return NORMAL_END;
}
