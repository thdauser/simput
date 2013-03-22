/************************************************************************
  atEcliptic.c   convert Equatorial <-> Ecliptic coordinates

  2007/04/07 Y.ISHISAKI
************************************************************************/

#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/* obliquity of the ecliptic at epoch J2000.0, epsilon = 23d26m21.448s
   Astrophysical Formulae Volume II Space, Time, Matter and Cosmology
   3rd Enlarged and Revised Edition (1999)
   Kenneth R. Lang, Springer, ISBN 3-540-64664-7 */

#define EPSILON_J2000	(23.4392911111*DEG2RAD)
#define cos_e			0.91748206206925896487
#define sin_e			0.39777715593173577844

#define TOO_SMALL		1.0e-15

/* convert Equatorial coordinate R.A. Dec. (alpha, delta) [deg] in J2000
   to Ecliptic coordinates (lambda, beta) [deg] */
void
atJ2000ToEcliptic(double alpha, double delta, double *l, double *b)
{
	double sin_a, cos_a, sin_d, cos_d;
	double cos_b_sin_l, cos_b_cos_l, sin_b, cos_b;
	double abs_cos_b_sin_l, abs_cos_b_cos_l;

	alpha *= DEG2RAD;
	delta *= DEG2RAD;
	sin_a = sin(alpha);
	cos_a = cos(alpha);
	sin_d = sin(delta);
	cos_d = cos(delta);

	cos_b_sin_l = cos_d * sin_a * cos_e + sin_d * sin_e;
	cos_b_cos_l = cos_d * cos_a;
	*l = atan2(cos_b_sin_l, cos_b_cos_l);
	if ( fabs(*l) < TOO_SMALL ) {
		*l = 0.0;
	}

	sin_b = sin_d * cos_e - cos_d * sin_a * sin_e;
	if ( fabs(sin_b) < TOO_SMALL ) {
		*l *= RAD2DEG;
		*b = 0.0;
		goto skip;
	}
	abs_cos_b_sin_l = fabs(cos_b_sin_l);
	abs_cos_b_cos_l = fabs(cos_b_cos_l);
	if ( abs_cos_b_sin_l < TOO_SMALL && abs_cos_b_cos_l < TOO_SMALL ) {
		*l = 0.0;
		*b = ( 0.0 < sin_b ) ? 90.0 : -90.0;
		return;
	} else if ( abs_cos_b_sin_l < abs_cos_b_cos_l ) {
		cos_b = cos_b_cos_l / cos(*l);
	} else {
		cos_b = cos_b_sin_l / sin(*l);
	}
	*b = atan2(sin_b, cos_b);

	*l *= RAD2DEG;		/* -180 < l <= 180 */
	*b *= RAD2DEG;		/* -90 <= b <= 90 */

 skip:
	while ( *l < 0.0 ) {
		*l += 360.0;	/* 0.0 <= l < 360.0 */
	}
}

/* convert Ecliptic coordinates (lambda, beta) [deg]
   to Equatorial coordinate R.A. Dec. (alpha, delta) [deg] in J2000 */
void
atEclipticToJ2000(double l, double b, double *alpha, double *delta)
{
	double sin_l, cos_l, sin_b, cos_b;
	double cos_d_sin_a, cos_d_cos_a, sin_d, cos_d;
	double abs_cos_d_sin_a, abs_cos_d_cos_a;

	l *= DEG2RAD;
	b *= DEG2RAD;
	sin_l = sin(l);
	cos_l = cos(l);
	sin_b = sin(b);
	cos_b = cos(b);

	cos_d_cos_a = cos_b * cos_l;
	cos_d_sin_a = cos_b * sin_l * cos_e - sin_b * sin_e;
	*alpha = atan2(cos_d_sin_a, cos_d_cos_a);
	if ( fabs(*alpha) < TOO_SMALL ) {
		*alpha = 0.0;
	}

	sin_d = cos_b * sin_l * sin_e + sin_b * cos_e;
	if ( fabs(sin_d) < TOO_SMALL ) {
		*alpha *= RAD2DEG;
		*delta = 0.0;
		goto skip;
	}
	abs_cos_d_sin_a = fabs(cos_d_sin_a);
	abs_cos_d_cos_a = fabs(cos_d_cos_a);
	if ( abs_cos_d_sin_a < TOO_SMALL && abs_cos_d_cos_a < TOO_SMALL ) {
		*alpha = 0.0;
		*delta = (0.0 < sin_d) ? 90.0 : -90.0;
		return;
	} else if ( abs_cos_d_sin_a < abs_cos_d_cos_a ) {
		cos_d = cos_d_cos_a / cos(*alpha);
	} else {
		cos_d = cos_d_sin_a / sin(*alpha);
	}
	*delta = atan2(sin_d, cos_d);

	*alpha *= RAD2DEG;		/* -180 < alpha <= 180 */
	*delta *= RAD2DEG;		/* -90 <= delta <= 90 */

 skip:
	while ( *alpha < 0.0 ) {
		*alpha += 360.0;	/* 0.0 <= alpha < 360.0 */
	}
}

#ifdef TEST
#include <stdio.h>
#include <stdlib.h>

/*
 * test routine
 */
int
main(int argc, char **argv)
{
	double alpha, delta, l, b, a, d;

	if ( argc < 2 ) {
		fprintf(stderr, "\
usage: atEcliptic alpha delta\n");
		return 1;
	}

	alpha = atof(argv[1]);
	delta = atof(argv[2]);

	atJ2000ToEcliptic(alpha, delta, &l, &b);
	atEclipticToJ2000(l, b, &a, &d);

	printf("   ( alpha , delta ) = ( %10.6f , %+10.6f )\n", alpha, delta);
	printf("-> ( lambda , beta ) = ( %10.6f , %+10.6f )\n", l, b);
	printf("-> ( alpha , delta ) = ( %10.6f , %+10.6f )\n", a, d);

	return 0;
}
#endif
