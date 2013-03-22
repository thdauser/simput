#include <math.h>
#include "atFunctions.h"
#include "atError.h"

/*
 * Equivalent to atPrecessRM(mjd, MJD_J2000, rm)
 */
void
atPrecessRMJ2000(double mjd, AtRotMat rm)
{
#ifdef ARCSEC2RAD
#undef ARCSEC2RAD
#endif
#define ARCSEC2RAD	(DEG2RAD/60.0/60.0)

	double t, zeta, z, theta;
	double cos_zeta, cos_z, cos_theta;
	double sin_zeta, sin_z, sin_theta;

	t = ( mjd - MJD_J2000 ) / 36525.0;

	zeta  = (2306.2181 + (0.30188 + 0.017998*t)*t)*t * ARCSEC2RAD;
	z     = (2306.2181 + (1.09468 + 0.018203*t)*t)*t * ARCSEC2RAD;
	theta = (2004.3109 - (0.42665 + 0.041833*t)*t)*t * ARCSEC2RAD;

	cos_zeta = cos(zeta);
	sin_zeta = sin(zeta);

	cos_z = cos(z);
	sin_z = sin(z);

	cos_theta = cos(theta);
	sin_theta = sin(theta);

	rm[0][0] =  cos_zeta*cos_theta*cos_z - sin_zeta*sin_z;
	rm[1][0] = -sin_zeta*cos_theta*cos_z - cos_zeta*sin_z;
	rm[2][0] = -sin_theta*cos_z;
	rm[0][1] =  cos_zeta*cos_theta*sin_z + sin_zeta*cos_z;
	rm[1][1] = -sin_zeta*cos_theta*sin_z + cos_zeta*cos_z;
	rm[2][1] = -sin_theta*sin_z;
	rm[0][2] =  cos_zeta*sin_theta;
	rm[1][2] = -sin_zeta*sin_theta;
	rm[2][2] =  cos_theta;
}

/*
 * Find Rotation Matrix for conversion of equatorial coordinate systems
 * correcting for precession (SAISA)
 */
int
atPrecessRM(
	double mjd0,	/* input: modified Julian Day for the original coord*/
	double mjd,		/* input: modified Julian Day for the new coord*/
	AtRotMat rm)	/* output: rotation matrix to correct precession*/
{
	AtRotMat RmAto2000, RmBto2000, Rm2000toB;

	atPrecessRMJ2000(mjd0, RmAto2000);
	atPrecessRMJ2000(mjd,  RmBto2000);
	ATInvRotMat(RmBto2000, Rm2000toB);
	ATRMProd(RmAto2000, Rm2000toB, rm);

	return NORMAL_END;
}
