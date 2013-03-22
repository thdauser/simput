#include "atFunctions.h"
#include "atError.h"
#include <math.h>

static int
moonag_(double mjd, double *ta, double *a, double *b, double *c, double *d,
		double *e, double *g, double *j, double *l, double *m, double *n,
		double *v, double *w)
{
	double tb;

	/* *ta= JULIAN CENTURY FROM 1900 JAN. 0.5
	 *  A= MEAN ANOMALY OF MOON
	 *  B= ARGUMENT LATITUDE OF MOON
	 *  C= MEAN LONGITUDE OF MOON
	 *  D= ELONGATION  OF MOON FROM SUN
	 *  E= LONGITUDE OF EARTH
	 *  G= MEAN ANOMALY OF SUN
	 *  J= MEAN ANOMALY OF JUPITER
	 *  L= MEAN LONGITUDE OF SUN
	 *  M= MEAN ANOMALY OF MARS
	 *  N= LONGITUDE OF ASCENDING NODE OF MOON'S ORBIT
	 *  V= MEAN ANOMALY OF VENUS
	 *  W= MEAN LONGITUDE OF VENUS
	 */

	*ta = (mjd - 15019.5) / 36525.;
	tb = *ta * *ta;

	*a = DEG2RAD*(*ta * 4.77e5 +296.1044608 + *ta * 198.849108 + tb * .009192);
	*b = DEG2RAD*(*ta * 483120. + 11.250889 + *ta * 82.02515 - tb * .003211);
	*c = DEG2RAD*(*ta * 480960. +270.434164 + *ta * 307.883142 - tb * .001133);
	*d = DEG2RAD*(*ta * 444960 + 350.737486 + *ta * 307.114217 - tb * .001436);
	*e = DEG2RAD*(*ta * 35640 + 98.998753 + *ta * 359.372886);
	*g = DEG2RAD*(*ta * 35999.04975 + 358.475833 - tb * 1.5e-4);
	*j = DEG2RAD*(*ta * 2880 + 225.444651 + *ta * 154.906654);
	*l = DEG2RAD*(*ta * 36000.76892 + 279.696678 + tb * 3.03e-4);
	*m = DEG2RAD*(*ta * 19080 + 319.529425 + *ta * 59.8585 + tb * 1.81e-4);
	*n = DEG2RAD*(259.183275 - *ta * 1800 - *ta * 134.142008 + tb * .002078);
	*v = DEG2RAD*(*ta * 58320 + 212.603219 + *ta * 197.803875 + tb * .001286);
	*w = DEG2RAD*(*ta * 58320 + 342.767053 + *ta * 199.211911 * 3.1e-4 * tb);
	return (NORMAL_END);
} /* moonag_ */

static int
moonth_(double ta, double a, double b, double c, double d,
		double e, double g, double j, double l, double m, double n,
		double v, double w, double *mx, double *my, double *mz)
{
	/* .... MOON THETA */
	*mx = sin(a + b - d * 4.) * -.00101;
	*mx -= sin(a - b - d * 4. - n) * .00102;
	*mx -= ta * .00103 * sin(a - b - n);
	*mx -= sin(a - g - b - d * 2. - n) * .00107;
	*mx -= sin(a * 2. - b - d * 4. - n) * .00121;
	*mx += sin(a * 3. + b + n) * .0013;
	*mx -= sin(a + b - n) * .00131;
	*mx += sin(a + b - d + n) * .00136;
	*mx -= sin(g + b) * .00145;
	*mx -= sin(a + g - b - d * 2.) * .00149;
	*mx += sin(g - b + d - n) * .00157;
	*mx -= sin(g - b) * .00159;
	*mx += sin(a - g + b - d * 2. + n) * .00184;
	*mx -= sin(b - d * 2. - n) * .00194;
	*mx -= sin(g - b + d * 2. - n) * .00196;
	*mx += sin(b - d) * .002;
	*mx -= sin(a + g - b) * .00205;
	*mx += sin(a - g - b) * .00235;
	*mx += sin(a - b * 3 - n) * .00246;
	*mx -= sin(a * 2 + b - d * 2) * .00262;
	*mx -= sin(a + g + b - d * 2) * .00283;
	*mx -= sin(g - b - d * 2 - n) * .00339;
	*mx += sin(a - b + n) * .00345;
	*mx -= sin(g - b + d * 2) * .00347;
	*mx -= sin(b + d + n) * .00383;
	*mx -= sin(a + g + b + n) * .00411;
	*mx -= sin(a * 2 - b - d * 2 - n) * .00442;
	*mx += sin(a - b + d * 2) * .00449;
	*mx -= sin(b * 3 - d * 2 + n) * .00456;
	*mx += sin(a + b + d * 2 + n) * .00466;
	*mx += sin(a * 2 - b) * .0049;
	*mx += sin(a * 2 + b) * .00561;
	*mx += sin(a - g + b + n) * .00564;
	*mx -= sin(a + g - b - n) * .00638;
	*mx -= sin(a + g - b - d * 2 - n) * .00713;
	*mx -= sin(g + b - d * 2) * .00929;
	*mx -= sin(a * 2 - b - n) * .00947;
	*mx += sin(a - g - b - n) * .00965;
	*mx += sin(b + d * 2) * .0097;
	*mx += sin(b - d + n) * .01064;
	*mx -= ta * .0125 * sin(b + n);
	*mx -= sin(g + b - d * 2 + n) * .01434;
	*mx -= sin(a + g + b - d * 2 + n) * .01652;
	*mx -= sin(a * 2 + b - d * 2 + n) * .01868;
	*mx += sin(a * 2 + b + n) * .027;
	*mx -= sin(a - b - d * 2) * .02994;
	*mx -= sin(g + b + n) * .03759;
	*mx -= sin(g - b - n) * .03982;
	*mx += sin(b + d * 2 + n) * .04732;
	*mx -= sin(b - n) * .04771;
	*mx -= sin(a + b - d * 2) * .06505;
	*mx += sin(a + b) * .13622;
	*mx -= sin(a - b - d * 2 - n) * .14511;
	*mx -= sin(b - d * 2) * .18354;
	*mx -= sin(b - d * 2 + n) * .20017;
	*mx -= sin(a + b - d * 2 + n) * .38899;
	*mx += sin(a - b) * .40248;
	*mx += sin(a + b + n) * .65973;
	*mx += sin(a - b - n) * 1.96763;
	*mx += sin(b) * 4.95372;
	*mx += sin(b + n) * 23.89684;
/* .... MOON RHO */
	*my = cos(a * 2 + g) * .05491;
	*my += cos(a + d) * .0629;
	*my -= cos(d * 4) * .06444;
	*my -= cos(a * 2 - g) * .06652;
	*my -= cos(g - d * 4) * .07369;
	*my += cos(a - d * 3) * .08119;
	*my -= cos(a + d * 4) * .09261;
	*my += cos(a - b * 2 + d * 2) * .10177;
	*my += cos(a + g + d * 2) * .10225;
	*my -= cos(a + g * 2 - d * 2) * .10243;
	*my -= cos(b * 2) * .12291;
	*my -= cos(a * 2 - b * 2) * .12291;
	*my -= cos(a + g - d * 4) * .12428;
	*my -= cos(a * 3) * .14986;
	*my -= cos(a - g + d * 2) * .1607;
	*my -= cos(a - d) * .16949;
	*my += cos(a + b * 2 - d * 2) * .17697;
	*my -= cos(a * 2 - d * 4) * .18815;
	*my -= cos(g * 2 - d * 2) * .19482;
	*my += cos(b * 2 - d * 2) * .22383;
	*my += cos(a * 3 - d * 2) * .22594;
	*my += cos(a * 2 + g - d * 2) * .24454;
	*my -= cos(g + d) * .31717;
	*my -= cos(a - d * 4) * .36333;
	*my += cos(a - g - d * 2) * .47999;
	*my += cos(g + d * 2) * .63844;
	*my += cos(g) * .8617;
	*my += cos(a - b * 2) * 1.50534;
	*my -= cos(a + d * 2) * 1.67417;
	*my += cos(a + g) * 1.99463;
	*my += cos(d) * 2.07579;
	*my -= cos(a - g) * 2.455;
	*my -= cos(a + g - d * 2) * 2.74067;
	*my -= cos(g - d * 2) * 3.83002;
	*my -= cos(a * 2) * 5.37817;
	*my += cos(a * 2 - d * 2) * 6.60763;
	*my -= cos(d * 2) * 53.97626;
	*my -= cos(a - d * 2) * 68.62152;
	*my -= cos(a) * 395.13669;
	*my += 3649.33705;
/*  ... MOON PHI */
	*mz = sin(a - g - b * 2 - n * 2) * -.001;
	*mz -= sin(a + g - d * 4) * .001;
	*mz += sin(a * 2 - g) * .001;
	*mz += sin(a - g + d * 2) * .00102;
	*mz -= sin(a * 2 - b * 2 - n) * .00106;
	*mz -= sin(a * 2 + n) * .00106;
	*mz -= sin(a + b * 2 - d * 2) * .00109;
	*mz -= sin(b * 2 - d + n * 2) * .0011;
	*mz += sin(d * 4) * .00112;
	*mz -= sin(a * 2 - n) * .00122;
	*mz -= sin(a * 2 + b * 2 + n) * .00122;
	*mz += sin(g + b * 2 - d * 2 + n * 2) * .00149;
	*mz -= sin(a * 2 - d * 4) * .00157;
	*mz += sin(a + g + b * 2 - d * 2 + n * 2) * .00171;
	*mz -= sin(a * 2 + g - d * 2) * .00175;
	*mz -= sin(g * 2 - d * 2) * .0019;
	*mz += sin(a + e * 16 - w * 18) * .00193;
	*mz += sin(a * 2 + b * 2 - d * 2 + n * 2) * .00194;
	*mz += sin(g - d * 2 - n) * .00201;
	*mz += sin(g + b * 2 - d * 2 + n) * .00201;
	*mz -= sin(a + g * 2 - d * 2) * .00207;
	*mz -= sin(g * 2) * .0021;
	*mz -= sin(d * 2 - n) * .00213;
	*mz -= sin(b * 2 + d * 2 + n) * .00213;
	*mz -= sin(a * 3 - d * 2) * .00215;
	*mz -= sin(a - d * 4) * .00247;
	*mz -= sin(a - b * 2 + d * 2) * .00253;
	*mz += ta * .00279 * sin(b * 2 + n * 2);
	*mz -= sin(a * 2 + b * 2 + n * 2) * .0028;
	*mz += sin(a * 3) * .00312;
	*mz -= sin(a + b * 2) * .00317;
	*mz -= sin(a + e * 16 - w * 18) * .0035;
	*mz += sin(g + b * 2 + n * 2) * .0039;
	*mz += sin(g - b * 2 - n * 2) * .00413;
	*mz -= sin(n * 2) * .0049;
	*mz -= sin(b * 2 + d * 2 + n * 2) * .00491;
	*mz += sin(g + d) * .00504;
	*mz += sin(a - d) * .00516;
	*mz -= sin(g + d * 2) * .00621;
	*mz += sin(a - b * 2 - d * 2 - n) * .00648;
	*mz += sin(a - d * 2 + n) * .00648;
	*mz += sin(a - g - d * 2) * .007;
	*mz += sin(a + d * 2) * .01122;
	*mz += sin(a - d * 2 - n) * .0141;
	*mz += sin(a + b * 2 - d * 2 + n) * .0141;
	*mz += sin(a - b * 2) * .01424;
	*mz += sin(a - b * 2 - d * 2 - n * 2) * .01506;
	*mz -= sin(b * 2 - d * 2) * .01567;
	*mz += sin(b * 2 - d * 2 + n * 2) * .02077;
	*mz -= sin(a + g) * .02527;
	*mz -= sin(a - n) * .02952;
	*mz -= sin(a + b * 2 + n) * .02952;
	*mz -= sin(d) * .03487;
	*mz += sin(a - g) * .03684;
	*mz -= sin(d * 2 + n) * .03983;
	*mz += sin(b * 2 - d * 2 + n) * .03983;
	*mz += sin(a + b * 2 - d * 2 + n * 2) * .04037;
	*mz += sin(a * 2) * .04221;
	*mz -= sin(g - d * 2) * .04273;
	*mz -= sin(a * 2 - d * 2) * .05566;
	*mz -= sin(a + g - d * 2) * .05697;
	*mz -= sin(a + b * 2 + n * 2) * .06846;
	*mz -= sin(a - b * 2 - n) * .08724;
	*mz -= sin(a + n) * .08724;
	*mz -= sin(b * 2) * .11463;
	*mz -= sin(g) * .18647;
	*mz -= sin(a - b * 2 - n * 2) * .20417;
	*mz += sin(d * 2) * .59616;
	*mz += sin(n) * 1.07142;
	*mz -= sin(b * 2 + n) * 1.07447;
	*mz -= sin(a - d * 2) * 1.28658;
	*mz -= sin(b * 2 + n * 2) * 2.4797;
	*mz += sin(a) * 6.32962;
	return (NORMAL_END);
} /* moonth_ */

/*
 * calc position of the moon from the earth center in J2000 coordinate.
 *
 *  modified from CSUB1.FORT77(SOLARSYS) in the Ginga library
 * 					n.kawai 93.01.13
 *
 *	2005/08/06	v2.5	Y.ISHISAKI
 *		static declarations of moonag_(), moonth_()
 */
int
atMoon(
	double mjd,	/* input: time in MJD */
	AtVect pos,	/* output: vector to the moon (km) */
	double *size,	/* output: their visual size (radian)*/
	double *phase,	/* output: phase (radian 0:new, pi:full)*/
	double *distan)	/* output: distance to moon (km)*/
{
	double ta, a, b, c, d, e, g, j, l, m, n, v, w, mx, my, mz;
	double r_xy, sin_delta, cos_delta, sin_c, cos_c;
	AtVect x_tod;

	moonag_(mjd, &ta, &a, &b, &c, &d, &e, &g, &j, &l, &m, &n, &v, &w);
	moonth_(ta, a, b, c, d, e, g, j, l, m, n, v, w, &mx, &my, &mz);

	r_xy = sqrt(my - mx * mx);
	sin_delta = mz / r_xy;
	cos_delta = sqrt(1. - sin_delta * sin_delta);
	sin_c = sin(c);
	cos_c = cos(c);

	/* R.A. of moon = mean longitude (c) + delta */

	*distan = EARTH_RADIUS * sqrt(my);
	x_tod[0] = EARTH_RADIUS * r_xy *( cos_delta * cos_c - sin_delta * sin_c );
	x_tod[1] = EARTH_RADIUS * r_xy *( sin_delta * cos_c + cos_delta * sin_c );
	x_tod[2] = EARTH_RADIUS * mx;

	*size = atan(MOON_RADIUS / *distan);
	*phase = fmod(d, TWO_PI);
	atPrecession(mjd, x_tod, MJD_J2000, pos);
	return (NORMAL_END);
} /* moon_ */
