/************************************************************************
  atPlanet.c   calculate position of the nine planets

  2007/04/07 Y.ISHISAKI
	- modified to use atPrecessRMJ2000() for faster calculation
************************************************************************/

#include "atFunctions.h"
#include "atError.h"
#include <math.h>

static void
plelms_(double mjd, double eccen[9], double mean[9], double peri[9],
		double node[9], double incli[9], double axis[9],
		double radi[9], double magni[9])
	   	/*   SET ORBITAL ELEMENTS OF 9 PLANETS */
{
	int j, k;
	double t, t1, t2;

	t = mjd - 33281.923;
	t1 = t * (t * 1.26013e-17 + 2.737909288e-5);
	t2 = t1 * t1;

	/* ... MERCURY */
	eccen[0] = t1 * 2.042e-5 + .20562411 - t2 * 3e-8;
	mean[0] = t1 * 2608.731797 - .729180963 + t2 * 8.242e-8;
	peri[0] = t1 * .00498049095 + .505094601 + t2 * 1.33324e-6;
	node[0] = .833197496 - t1 * .0021919881 - t2 * 1.57564e-6;
	incli[0] = .122238982 - t1 * 1.0510761e-4 + t2 * 1.454e-8;
	axis[0] = .3870986011;
	radi[0] = 2439.;
	magni[0] = 1.16;

	/* ... VENUS */
	eccen[1] = .00679676 - t1 * 4.773e-5 - t2 * 9e-8;
	mean[1] = t1 * 1021.306597 - .848503209 + t2 * 2.259232e-5;
	peri[1] = t1 * .00499862298 + .953579254 - t2 * 2.077911e-5;
	node[1] = 1.330454065 - t1 * .00485584535 - t2 * 1.78896e-6;
	incli[1] = .059237299 - t1 * 1.798659e-5 - t2 * 5.6723e-7;
	axis[1] = .72333162836;
	radi[1] = 6052.;
	magni[1] = -4.00;

	/* ... EARTH   (RADI =  SUN'S) */
	eccen[2] = .01673012 - t1 * 4.192e-5 - t2 * 1.3e-7;
	mean[2] = t1 * 628.288592 - .036149841 - t2 * 2.876525e-6;
	peri[2] = t1 * .0097864 - 1.262568078 - t2 * 2.55497e-6;
	node[2] = 3.044140013 - t1 * .00421225519 + t2 * 2.0847e-7;
	incli[2] = t1 * 2.2713521e-4 - t2 * 2.618e-7;
	axis[2] = 1.00000023;
	radi[2] = 695989.;
	magni[2] = -26.55;		/* magnitude of the sun */

	/* ... MARS */
	eccen[3] = t1 * 9.056e-5 + .09335426 - t2 * 7e-8;
	mean[3] = t1 * 334.0465218 - 3.326348457 + t2 * 3.08342e-6;
	peri[3] = t1 * .01288077229 + 4.991059739 + t2 * 7.98973e-6;
	node[3] = .858200646 - t1 * .00514920611 - t2 * 1.107314e-5;
	incli[3] = .0032288058 - t1 * 1.4539562e-4 - t2 * 3.9755e-7;
	axis[3] = 1.523688174;
	radi[3] = 3397.;
	magni[3] = -1.30;

	/* ... JUPITER */
	eccen[4] = t1 * 4.7756e-5 + .04827062 - t2 * 2.2676e-5;
	mean[4] = t1 * 52.96844784 + 5.28677165 + t2 *
		9.996858e-5;
	peri[4] = -1.50945736 - t1 * 1.7239975e-4 + t2 * 1.885925e-5;
	node[4] = t1 * 3.204618e-5 + 1.741507757 + t2 * 9.40539e-6;
	incli[4] = t1 * 7.8055e-7 + .0228307 + t2 * 3.6991e-7;
	axis[4] = 5.202833481;
	radi[4] = 71398.;
	magni[4] = -8.93;

	/* ... SATURN */
	eccen[5] = .05604508 - t1 * 2.5595e-5 - t2 * 1.6172e-5;
	mean[5] = t1 * 21.328205912 + 1.165262433 - t2 *
		5.4580324e-4;
	peri[5] = t1 * 4.2707237e-4 - .38321152 + t2 * 2.6548394e-4;
	node[5] = t1 * 3.005845e-5 + 1.980742073 - t2 * 5.657776e-5;
	incli[5] = t1 * 8.69756e-6 + .043422822 + t2 * 3.56623e-6;
	axis[5] = 9.538762055;
	radi[5] = 6e4;
	magni[5] = -8.68;

	/* ... URANUS */
	eccen[6] = .04613734 - t1 * 4.8118e-5 - t2 * 1.5396e-5;
	mean[6] = t1 * 7.479637598 - 1.273611261 + t2 *
		8.2743151e-4;
	peri[6] = 1.716582467 - t1 * .00237485982 - t2 * 8.1424458e-4;
	node[6] = t1 * 6.429559e-4 + 1.28641873 + t2 * 3.975476e-6;
	incli[6] = .013501673 - t1 * 1.72933e-5 - t2 * 8.7412e-7;
	axis[6] = 19.19139128;
	radi[6] = 25400.;
	magni[6] = -6.85;

	/* ... NEPTUNE */
	eccen[7] = t1 * .001095407 + .00971449 + t2 * 3.62034e-4;
	mean[7] = t1 * 3.994465294 + 2.724754505 + t2 *
		.0483557044;
	peri[7] = -1.62194737 - t1 * .18123685121 - t2 * .048352097;
	node[7] = t1 * 4.467073e-5 + 2.290559396 - t2 * 1.844231e-5;
	incli[7] = .03096505 - t1 * 3.001e-6 + t2 * 3.6216e-7;
	axis[7] = 30.06106906;
	radi[7] = 24300.;
	magni[7] = -7.05;

	/* ... PLUTO */
	eccen[8] = t1 * 4.97082e-4 + .24824802 + t2 * 5.63208e-4;
	mean[8] = t1 * 2.553487672 - .999328454 + t2 *
		.00686515565;
	peri[8] = 1.977072713 - t1 * .01828324506 - t2 *
		.00676126008;
	node[8] = t1 * 1.027805e-5 + 1.913508742 + t2 * 5.6820159e-5;
	incli[8] = t1 * 1.5198909e-4 + .29929226 + t2 * 5.269925e-5;
	axis[8] = 39.52940243;
	radi[8] = 1500.;
	magni[8] = -1.;			/* iikagen */

	for (j = 0; j < 9; ++j) {
		k = (int) (mean[j] / TWO_PI);
		mean[j] -= k * TWO_PI;
		k = (int) (peri[j] / TWO_PI);
		peri[j] -= k * TWO_PI;
		k = (int) (node[j] / TWO_PI);
		node[j] -= k * TWO_PI;
		k = (int) (incli[j] / TWO_PI);
		incli[j] -= k * TWO_PI;
	}
}

static void
ephcns_(double peri, double node, double incli, AtVect f, AtVect q)
{
	static double k4 = .9174369375;
	static double k5 = .3978812206;

	double r1, r2, r3, r4, r5, r6, r7, r8, r9;

	r1 = sin(peri);
	r2 = sin(node);
	r3 = sin(incli);
	r4 = cos(peri);
	r5 = cos(node);
	r6 = cos(incli);
	f[0] = r5 * r4 - r1 * r6 * r2;
	q[0] = -r1 * r5 - r4 * r6 * r2;
	r7 = r6 * r5 * k4;
	r8 = r2 * k4;
	r9 = r3 * k5;
	f[1] = r1 * r7 + r4 * r8 - r1 * r9;
	q[1] = r4 * r7 - r1 * r8 - r4 * r9;
	r7 = r3 * k4;
	r8 = r6 * r5 * k5;
	r9 = r2 * k5;
	f[2] = r1 * r7 + r1 * r8 + r4 * r9;
	q[2] = r4 * r7 + r4 * r8 - r1 * r9;
} /* ephcns_ */

static void
hcent_(double axis, double eccen, double e, AtVect f, AtVect q, AtVect x)
{
	double b, w1, w2;

	b = axis * sqrt(1. - eccen * eccen);
	w1 = cos(e) - eccen;
	w2 = sin(e);
	x[0] = axis * f[0] * w1 + b * q[0] * w2;
	x[1] = axis * f[1] * w1 + b * q[1] * w2;
	x[2] = axis * f[2] * w1 + b * q[2] * w2;
}

/*
 * calc position of the planets from the earth center in J2000 coordinate.
 * array index 0-8 corresponds to the planets except 2 for the sun.
 *
 *  modified from CSUB1.FORT77(SOLARSYS) in the Ginga library
 * 					n.kawai 93.01.13
 *
 *	2005/08/06	v2.5	Y.ISHISAKI
 *		static declarations of plelms_(), ephcns_(), hcent_()
 */
int
atPlanet(
	double mjd,		/* input: time in MJD */
	AtVect pos[9],	/* output: vector to the planets and sun (A.U.) */
	double size[9],	/* output: their visual size (radian)*/
	double mag[9])	/* output: their visual magnitude*/
{
	static int first_time = 1;
	static AtRotMat rm1950_2000;

	double e, d_earth, d_sun;
	AtVect f, q, x, y;
	int j;

	double eccen[9], mean[9], peri[9], node[9], incli[9], axis[9],
		radi[9], magni[9];

	if ( first_time ) {
		atPrecessRMJ2000(MJD_B1950, rm1950_2000);
		first_time = 0;
	}

	plelms_(mjd, eccen, mean, peri, node, incli, axis, radi, magni);

	j = 2;
	ephcns_(peri[j], node[j], incli[j], f, q);
	atKepler(mean[j], eccen[j], &e);
	hcent_(axis[j], eccen[j], e, f, q, x);
	ATInvVect(x, y);
	ATRotVect(rm1950_2000, y, pos[j]);
	size[j] = atan(radi[j] / (ATNorm(x) * AU));
	mag[j] = magni[j];

	for (j = 0; j < 9; j++) {
		if (j != 2) {
			ephcns_(peri[j], node[j], incli[j], f, q);
			atKepler(mean[j], eccen[j], &e);
			hcent_(axis[j], eccen[j], e, f, q, x);
			ATRotVect(rm1950_2000, x, y);
			ATAddVect(y, pos[2], pos[j]);
			d_sun = ATNorm(x);
			d_earth = atNorm(pos[j]);
			size[j] = atan(radi[j] / (d_earth * AU));
			mag[j] = magni[j] + 5. * log10(d_earth * d_sun);
		}
	}

	return NORMAL_END;
}
