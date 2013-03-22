/************************************************************************
  atGeomag.c	calc Geomagnetic Field (GA)

	Modified from "GA & GASET" in CSUB1 for Ginga

	2005/12/04 Y.ISHISAKI	version 2.7
		change float -> double, remove unnecessary static declarations

 VERSION 2.1  2003/12/20  by Prof. I. Kondo
    ALTERED TO CONVERT SCHIMDT NORMALIZED TO GAUSE NORIMALIZED
    IGRF 2000 WITH ITS DERIVATIVE; MAX. ORDER 10 IS USED.
 VERSION 2.2  2005/12/15  by S. Yamauchi
    COEFFICIENTS OF IGRF 2000 WITH MAX. ORDER OF 13 IS USED.
 VERSION 2.3  2006/1/26  by S. Yamauchi
    COEFFICIENTS OF IGRF 2005 WITH MAX. ORDER OF 13 IS USED.

	2006/02/03 Y.ISHISAKI	version 2.8
		change float -> double
		remove unused variables, unnecessary static declarations
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "atFunctions.h"
#include "atError.h"

static struct {
	double cnm[105], knm[105], pa[105], pb[105];
	int nm0[14];
	int nmax1;
}  atConsnm;

int
atGeomag(
	AtPolarVect *xp,/* input: polar vector for geodetic position */
	AtVect xs,		/* input: sidereal vector to the satellite */
	AtVect field)	/* output: sidereal direction of field line */
{
/* Initialized data */
	static AtVect south = { 0.0, 0.0, -1.0 };

/* Local variables */
	int n, m, nm, nmax, nm1, nm2, nm3;
	double aras, bcos, bsin, cosq, sinq, sumx, sumy, sumz;
	double x, y, z, arass, qq, qx, qy;
	double psaisg , tethas, acs, asn, tn11, tn21;
	double fm, fn;
	double p[105], dp[105];

	AtVect v;
	AtRotMat rm, rmi;

	tethas = xp->lat;
	qx = sin(tethas);
	qy = cos(tethas);
	psaisg = xp->lon;
	if ( psaisg < 0. ) {
		psaisg += TWO_PI;
	}
	p[0] =  1.0;
	p[1] =   qx;
	p[2] =   qy;
	dp[0] = 0.0;
	dp[1] = -qy;
	dp[2] =  qx;
	if (qy < EPS) return -11;

	nm = 0;		/* for gcc warning */
	nmax = atConsnm.nmax1 - 1;
	for (n = 2; n <= nmax; n++) {
 		fn = n;
 		tn11 = fn * 2. - 2.;
 		tn21 = sqrt(fn * 4. - 2.);
 		if (n == 2) {
	 		tn21 = sqrt(fn * 8. - 4.);
 		}
 		for (m = 0; m < n; m++) {
	 		nm1 = m + atConsnm.nm0[n - 1];
	 		nm2 = m + atConsnm.nm0[n - 2];
	 		nm  = m + atConsnm.nm0[n];
	 		fm  = m;
	 		p[nm] = p[nm1] * qx - atConsnm.knm[nm] * p[nm2];
	 		dp[nm] = qx * dp[nm1] - qy * p[nm1] - atConsnm.knm[nm] * dp[nm2];
		}
/* n = m */
		nm += 1;
		nm3 = nm - n -1;
		p[nm] = qy * p[nm3];
		dp[nm] = qy * dp[nm3] + qx * p[nm3];
	}
/*  conpute mag field.  */
	x = 0.;
	y = 0.;
	z = 0.;
	aras = EARTH_RADIUS / xp->r;
	arass = aras * aras;
	for (n = 1; n <= nmax; n++) {
 		arass *= aras;
 		fn = n;
 		tn11 = fn * (fn + 1.) / 2.;
 		sumx = 0.;
 		sumy = 0.;
 		sumz = 0.;
 		for (m = 0; m <= n; m++) {
	 		nm = m + atConsnm.nm0[n];
	 		fm = m;
	 		qq = fm * psaisg;
	 		cosq = cos(qq);
	 		sinq = sin(qq);
	 		acs = atConsnm.pa[nm] * cosq;
	 		asn = atConsnm.pa[nm] * sinq;
	 		bsin = atConsnm.pb[nm] * sinq;
	 		bcos = atConsnm.pb[nm] * cosq;
	 		sumx -= (acs + bsin) * dp[nm];
	 		sumy += fm * (asn - bcos) * p[nm] / qy;
	 		sumz += (acs + bsin) * p[nm];
/* L500: */
 		}
/* X: SOUTH, Y: EAST, Z: UPWARD COMPONENTS OF GEOMAG FIELD. */
		 x += arass * sumx;
		 y += arass * sumy;
 		 z += (fn + 1.) * arass * sumz;
/* L1000: */
	}

/*	   COMPUTE COMPONENTS OF GEOMAG. FIELD IN EQUATORIAL COORDINATES */
/*	   XYZ(1): COMPONENT TO VERNAL EQUINOX (R.A.=0) */
/*	   XYZ(2): COMPONENT TO Y-AXIS (R.A.=90) */
/*	   XYZ(3): COMPONENT PARALLEL TO EARTH'S ROTATION AXIS. */

	v[0] = x;
	v[1] = y;
	v[2] = z;
	atSetRotMatZX(xs, south, rm);
	atInvRotMat(rm, rmi);
	atRotVect(rmi, v, field);

	return 0;
} /* ga_ */


/* Factorial computation routine. */
static double
fact(long n)
{
	long nn, n1;
	double fac2;

	nn = n;
	n1 = nn -1;
	fac2 = nn;
	if (n == 0) return (1.0);
	if (n == 1) return (1.0);
	if (n <= -1) return (-999.9);
	for (nn = n1 ; nn > 1; nn--) {
		fac2 = fac2 * (double)(nn);
	}
	return fac2;
}

/*
 * Set up Table of Geomagnetic Field (GASET)
 */
int
atGeomagSet(
	double mjd,	/* input: time in MJD */
	int nmax)	/* input: maximum order of spherical harmonics (not used) */
{
/* Initialized data */
	static double g[105] = {
		0,
		-.295568, -.016718,
		-.023405,  .030470,  .016569,
		 .013357, -.023053,  .012468,  .006744,
		 .009198,  .007982,  .002115, -.003795,  .001002,
		-.002276,  .003544,  .002088, -.001366, -.001683, -.000141,
		 .000729,  .000696,  .000766, -.001511, -.000150,  .000147, -.000864,
		 .000798, -.000744, -.000014,  .000386,  .000123,  .000094,  .000055,  .000020,
		 .000248,  .000077, -.000114, -.000068, -.000180,  .000100,  .000094, -.000114, -.000050,
		 .000056,  .000098,  .000036, -.000070,  .000050, -.000108, -.000013,  .000087, -.000067, -.000092,
		-.000022, -.000063,  .000016, -.000025, -.000001,  .000030,  .000003,  .000021,  .000039, -.000001, -.000022,
		 .000029, -.000016, -.000017,  .000015, -.000002,  .000002, -.000007,  .000005,  .000018,  .000001,  .000010,  .000041,
		-.000022, -.000003,  .000003,  .000009, -.000004,  .000010, -.000004,  .000005, -.000003, -.000004,  .000000, -.000004,  .000000,
		-.000002, -.000009,  .000003,  .000003, -.000004,  .000012, -.000004,  .000007, -.000003,  .000004, -.000001,  .000004, -.000001, -.000003
	};

	static double h[105] = {
		0,
		0,	 .050800,
		0,	-.025949, -.005167,
		0,	-.002004,  .002693, -.005245,
		0,	 .002814, -.002258,  .001457, -.003047,
		0,	 .000427,  .001798, -.001230, -.000195,  .001036,
		0,	-.000202,  .000547,  .000637, -.000634,  .000000,  .000503,
		0,	-.000614, -.000225,  .000069,  .000254,  .000109, -.000264, -.000048,
		0,	 .000112, -.000210,  .000097, -.000198,  .000161,  .000077, -.000128, -.000001,
		0,	-.000201,  .000129,  .000127, -.000067, -.000081,  .000081,  .000029, -.000079,  .000059,
		0,	 .000024,  .000002,  .000044,  .000047, -.000065, -.000010, -.000034, -.000009, -.000023, -.000080,
		0,	 .000003,  .000014, -.000007, -.000024,  .000009, -.000006, -.000027, -.000010, -.000015, -.000020, -.000014,
		0,	-.000005,  .000003,  .000023, -.000027,  .000006,  .000004, 0.000000, 0.000000,  .000003, -.000008, -.000004,  .000010,
		0,	-.000007,  .000003,  .000017, -.000005, -.000010,  .000000,  .000007,  .000002,  .000006,  .000004, -.000002, -.000005, -.000010
	};

	static double dg[45] = {
		0,
		 .000088,  .000108,
		-.000150, -.000069, -.000010,
		-.000003, -.000031, -.000009, -.000068,
		-.000025,  .000028, -.000071,  .000059, -.000032,
		-.000026,  .000004, -.000030, -.000012,  .000002, -.000006,
		-.000008,  .000002, -.000002,  .000021, -.000021, -.000004,  .000013,
		-.000004,  .000000, -.000002,  .000011,  .000006,  .000004, -.000005,  .000009,
		-.000002,  .000002, -.000002,  .000002, -.000002,  .000002,  .000005, -.000007,  .000005
	};

   static double dh[45] = {
	   0,
	   0,	-.000213,
	   0,	-.000233, -.000140,
	   0,	 .000054, -.000065, -.000020,
	   0,	 .000020,  .000018,  .000056,  .000000,
	   0,	 .000001,  .000018,  .000020,  .000045, -.000010,
	   0,	-.000004, -.000019, -.000004, -.000004, -.000002,  .000009,
	   0,	 .000008,  .000004,  .000001,  .000002, -.000009, -.000003,  .000003,
	   0,	-.000002,  .000002,  .000002,  .000004,  .000002, -.000003,  .000005,  .000004
   };

	static AtTimeD ito = { 2005, 1, 1, 0, 0, 0, 0.0 };

/* System generated locals */
	int  m, n, nm, nm1;
	double year;
	double fm, fn, tfn, fn12, twp;
	double mjd0, fac2n, twopn1, facn1, fsq1, tfn13, fac3n;

/*
 SUBROUTINE GASET TO SET UP CONSTANT TABLE FOR THE SUBROUTINE GA PRAMETERS:
  MJD (IN )  : MODIFIED JULIAN DAY
  NMAX (IN ) : MAXIMUM NUMBER OF THE SPHERICAL HARMONICS
 VERSION 1.1  1986/7/12
	 SCHMIDT COEFFICIENTS FOR GEOMAG. FIELD. (IGRF 2000, IN GAUSS)
	 SECULAR VARIATION OF SPHERICAL HARMONICS. (GAUSS/YEAR)
	 PREPARE CONSTANT TABLE FOR SPHERICAL HARMONICS EQUATIONS.
	 SET MAX ORDER OF HARMONICS TO 10
  VERSION 2.1  2003/12/20  by Prof. I. Kondo
	 ALTERED TO CONVERT SCHIMDT NORMALIZED TO GAUSE NORIMALIZED
	 IGRF 2000 WITH ITS DERIVATIVE; MAX. ORDER 10 IS USED.
  VERSION 2.2  2005/12/15  by S. Yamauchi
	 COEFFICIENTS OF IGRF 2000 WITH MAX. ORDER OF 13 IS USED.
  VERSION 2.3  2006/01/26  by S. Yamauchi
	 COEFFICIENTS OF IGRF 2005 WITH MAX. ORDER OF 13 IS USED.
*/

	nmax = 13;
	atConsnm.nmax1 = nmax + 1;
	atConsnm.nm0[0] = 0;
	atConsnm.nm0[1] = 1;
	atConsnm.nm0[2] = 3;
	atConsnm.knm[0] = 0;
	atConsnm.knm[1] = 0;
	atConsnm.knm[2] = 1.0;
	atConsnm.knm[3] = 1.0/3.0;
	atConsnm.knm[4] = 0;
	atConsnm.knm[5] = -1.0;
	atConsnm.cnm[0] = 1.0;
	atConsnm.cnm[1] = 1.0;
	atConsnm.cnm[2] = 1.0;
	atConsnm.cnm[3] = 3.0/2.0;
	atConsnm.cnm[4] = sqrt(3.0);
	atConsnm.cnm[5] = sqrt(3.0)/2.0;
	twopn1 = 2.0;
	facn1  = 1.0;
	fsq1 = sqrt (2.0 / 3.0 / 1.0);
	twp = 2.0;
	nm = 5;
/*  L1000  start  */
	for (n = 3; n <= nmax; n++) {
/*	 set coeffs for m = 0  */
		nm += 1;
		atConsnm.nm0[n] = nm;
		nm1 = nm - n;
		fn = n;
		tfn = 2.0 * fn;
		twp = 2.0 * twp;
		fn12 = (fn - 1.0) * (fn - 1.0);
		fac2n = (tfn - 1.0) * (tfn - 2.0);
		fac3n = fact(2*n - 1) / twp / fact(n-1) ;
		tfn13 = (tfn - 1.0 ) * ( tfn - 3.0) ;
		atConsnm.cnm[nm] = fac3n / fact(n) ;
		atConsnm.knm[nm] = fn12 /(tfn -1.0) /(tfn - 3.0) ;
/*  L1100  start  */
		for (m = 1; m <= n; m++) {
			fm = m ;
			nm  += 1;
			nm1 += 1;
			fsq1 = sqrt(2.0 / fact(n + m) / fact(n - m)) ;
			atConsnm.cnm[nm] = fsq1 * fac3n ;
			atConsnm.knm[nm] = (fn12 - fm * fm) / tfn13 ;
		}
/*  L1100  end  */
	}
/*  L1000  end  */
/* L2000:  start */
/*	 SET UP SCHMIDT NORMALIZED COEFFICIENTS FOR GEOMAGNETIC FIELD. */
/*	 COMPUTE FOR ANY DATE FROM IGRF(2005). */
	atMJulianD(&ito, &mjd0);
/*	 GET DIFFERENCE IN YEAR OF THE DATE FROM EPOCH DAY. */
	year = mjd - mjd0;
	year = year / 365.25;
/*	 COMPUTE GEOMAG COEFFICIENTS. */
/*	 DERIVATIVES ARE UP TO ORDER = 13  */
	nm = 0;
	for (n = 1; n <= 8; ++n) {
   		for (m = 0; m <= n; ++m) {
	 		atConsnm.pa[nm] = atConsnm.cnm[nm] * (g[nm] + year * dg[nm]);
	 		atConsnm.pb[nm] = atConsnm.cnm[nm] * (h[nm] + year * dh[nm]);
	 		++nm;
	 	}
/* L3005: */
	}
	for (n = 9; n <= 13; ++n) {
		for (m = 1; m <= n; ++m) {
			atConsnm.pa[nm] = atConsnm.cnm[nm] * g[nm] ;
			atConsnm.pb[nm] = atConsnm.cnm[nm] * h[nm] ;
	 		++nm;
	 	}
	}
/* L3000: */
	return 0;
}
