#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * calc satellite position (SATPNT) with the earth center as origin
 *			based on  CSUB1 VERSION 2.1   1986/7/12
 *
 *	2005/08/29	v2.5	Y.ISHISAKI
 *		modified for faster calculation
 */
int atSatPos(
	double mjd,	/* input: time in MJD */
	AtVect xyz)	/* output: sidereal vector to the satellite
			 from the earth center*/
{
    double e, u, timed, timem, amean0, omebig, omeltl, semaxs;
    AtVect xys;

    timed = mjd - atElement.mjdz;
    timem = timed * 1440.;
    semaxs = atElement.semiax + atElement.adot * timed;
    omebig = atElement.ragome + atElement.ragdot * timed;
    omeltl = atElement.smaome + atElement.smodot * timed;
    amean0 = atElement.omean0 + (atElement.znbadt*.5*timed
			      + atElement.znbar) *timem;

	e = atElement.eccent;
    atKepler(amean0, e, &u);

	xys[0] = semaxs * ( cos(u) - e );
	xys[1] = semaxs * sqrt(1 - e*e ) * sin(u);
	xys[2] = 0.0;

    atOrbPlane(xys, omeltl, omebig, atElement.aincln, xyz);
    return 0;
}
