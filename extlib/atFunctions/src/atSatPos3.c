/*
 * calc satellite position (SATPNT) with the earth center as origin
 *			based on  CSUB1 VERSION 2.1   1986/7/12
 *
 * 2005/08/29	v2.5	Y.ISHISAKI
 *		modified for faster calculation
 *
 * 2009/07/24	Y.ISHISAKI	version 3.1
 *		modified to use atElement3 in argument
 */

#include <math.h>
#include "atFunctions.h"
#include "atError.h"

int
atSatPos3(
	AtElement3 *el,	/* input: contents of orbital elements to be read */
	double mjd,		/* input: time in MJD */
	AtVect xyz		/* output: sidereal vector to the satellite
			 				   from the earth center */
)
{
    double e, u, timed, timem, amean0, omebig, omeltl, semaxs;
    AtVect xys;

    timed = mjd - el->mjdz;
    timem = timed * 1440.;
    semaxs = el->semiax + el->adot * timed;
    omebig = el->ragome + el->ragdot * timed;
    omeltl = el->smaome + el->smodot * timed;
    amean0 = el->omean0 + (el->znbadt*.5*timed + el->znbar) *timem;

	e = el->eccent;
    atKepler(amean0, e, &u);

	xys[0] = semaxs * ( cos(u) - e );
	xys[1] = semaxs * sqrt(1 - e*e ) * sin(u);
	xys[2] = 0.0;

    atOrbPlane(xys, omeltl, omebig, el->aincln, xyz);
    return 0;
}
