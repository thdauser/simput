#include "atFunctions.h"
#include "atError.h"

/*
 * Find Euler Angles for conversion of equatorial coordinate systems
 * correcting for precession (SAISA)
 */
int atPrecessEuler(
	double mjd0,	/* input: modified Julian Day for the original coord*/
	double mjd,	/* input: modified Julian Day for the new coord*/
	AtEulerAng *ea)	/* output: Euler Ang to correct precession*/
{
#if 0	/* old version */
    double u, u0;

    u0 = (mjd0 - 15020.) / 36524.22;
    u = (mjd - mjd0) / 36524.22;
    ea->phi = -(u0*1.396 + 2304.25 + (u*.018 + .302) *u) * u * DEG2RAD /3600.;
    ea->psi =  ea->phi - u * .791 * u * DEG2RAD / 3600.;
    ea->theta = (2004.682 - u0*.853 + (-.426 - u*.042)* u)* u* DEG2RAD /3600.;
#else
	AtRotMat rm;

	atPrecessRM(mjd0, mjd, rm);
	atRMToEuler(rm, ea);
#endif
    return 0;
}
