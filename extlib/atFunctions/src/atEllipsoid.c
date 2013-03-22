#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert polar geodetic coordinate latitude radial distance
 * to the geographic latitude and altitude from the earth surface
 * correcting for the ellipsoidal shape of the earth
 *
 * 2005/08/08	v2.5	Y.ISHISAKI
 *		rewrite to be consistent with atGeographic.c
 *
 * 2005/10/24	v2.6	Y.ISHISAKI
 *		bug fix in calculating *heigh, which was negative in *latt < 0 in v2.5
 */
int
atEllipsoid(
	AtPolarVect *xp,	/* input: vector in celestial coordinate */
	double *latt,		/* output: latitude on the earth surface */
	double *heigh		/* output: altitude from the earth surface */
)
{
	static double eps = 1e-10;

/* Local variables */
	AtVect vec;
	int isign;
	double x, z, b, b2, det, w0, w1, wm, rm, xs, rs;
	double sin_t0, sin_t1, sin_tm, sin_latt;
	double cos_t0, cos_t1, cos_tm, cos_latt;
	double x0, x1, xm, dx0, dx1, dxm;
	double z0, z1, zm, dz0, dz1, dzm;

	atPolToVect(xp, vec);
	x = sqrt(vec[0]*vec[0] + vec[1]*vec[1]) / EARTH_RADIUS;
	z = vec[2] / EARTH_RADIUS;
	b2 = 1.0 - EARTH_E2;
	b = sqrt(b2);

	if ( fabs(z) < eps ) {
		*heigh = EARTH_RADIUS * ( sqrt(x*x + z*z) - 1.0 );
		*latt = 0;
		if ( fabs(x) < eps ) {
			return NULL_VECTOR;
		}
		return NORMAL_END;
	}

	if ( fabs(x) < eps ) {
		*heigh = EARTH_RADIUS * ( sqrt(x*x + z*z) - b );
		if ( 0 < z ) {
			*latt = PI / 2;
		} else {
			*latt = - PI / 2;
		}
		return NORMAL_END;
	}

	isign = 1;
	if ( z < 0.0 ) {
		isign = -1;
		z = -z;
	}

	sin_t0 = 0.0;
	cos_t0 = 1.0;
	x0 = cos_t0;
	z0 = b * sin_t0;
	dx0 = x - x0;
	dz0 = z - z0;
	w0 = fabs( ( - dx0*sin_t0 + b*cos_t0*dz0) / sqrt(dx0*dx0 + dz0*dz0) );

	sin_t1 = 1.0;
	cos_t1 = 0.0;
	x1 = cos_t1;
	z1 = b * sin_t1;
	dx1 = x - x1;
	dz1 = z - z1;
	w1 = fabs( ( - dx1*sin_t1 + b*cos_t1*dz1) / sqrt(dx1*dx1 + dz1*dz1) );

	sin_tm = z / sqrt( b2*x*x + z*z );
	cos_tm = sqrt( 1 - sin_tm * sin_tm );
	xm = cos_tm;
	zm = b * sin_tm;
	dxm = x - xm;
	dzm = z - zm;

	for (;;) {

		det = - sin_tm * dxm + b * cos_tm * dzm;
		rm = sqrt(dxm*dxm + dzm*dzm);

		if ( rm < eps ) {
			break;
		}

		wm = fabs( det / rm );

		if ( 0 < det ) {
			sin_t0 = sin_tm;
			cos_tm = cos_tm;
			x0 = xm;
			z0 = zm;
			dx0 = dxm;
			dz0 = dzm;
			w0 = wm;
		} else {
			sin_t1 = sin_tm;
			cos_t1 = cos_tm;
			x1 = xm;
			z1 = zm;
			dx1 = dxm;
			dz1 = dzm;
			w1 = wm;
		}

		if ( wm < eps ) {
			break;
		}

/*		printf("sin_tm=%.12f, wm=%.12f\n", sin_tm, wm);*/

		sin_tm = (w1 * sin_t0 + w0 * sin_t1) / ( w0 + w1 );
		cos_tm = sqrt( 1 - sin_tm * sin_tm );
		xm = cos_tm;
		zm = b * sin_tm;
		dxm = x - xm;
		dzm = z - zm;

	}

/*
   Here, P = (xm, zm) represents the vertical point (CHOKKA-TEN)
   of the satellite on the ellipsoid of the major axis a = 1.
   The vertical line intersects with x-axis at S = (a - b^2 xm, 0).
   Therefore, tan(latt) = zm / (b^2 xm) = zm / rs.
*/
	xs = b2 * xm;
	rs = sqrt(zm*zm + xs*xs);
	sin_latt = zm / rs;
	cos_latt = xs / rs;
	if ( sin_latt > cos_latt ) {
		*heigh = EARTH_RADIUS * ( z - b * sin_tm ) / sin_latt;
	} else {
		*heigh = EARTH_RADIUS * ( x - cos_tm ) / cos_latt;
	}

	if ( isign < 0 ) {
		zm = - zm;		/* latt < 0 */
	}

	*latt  = atan2(zm, xs);

	return NORMAL_END;
}
