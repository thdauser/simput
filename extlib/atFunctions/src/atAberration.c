/************************************************************************
  atAberration.c   correct annual aberration

  2006/02/08 R.Fujimoto
	- bug fixed at calculating lambda (celestial longitude) of the Sun
	  in atAberration and atInvAberration, which was pointed out by
	  the ASTRO-F team. Thanks.
  2007/04/07 Y.ISHISAKI
	- update atAberration() and atInvAberration() for faster calculation
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "atFunctions.h"

#define EPSILON		(23.440527*DEG2RAD)
#define COS_EPSILON	0.91747348167110063813
#define SIN_EPSILON	0.39779694623049654754
#define KAPPA		9.936508e-5

static void
setAberrMat(AtVect x, AtRotMat rm)
{
	rm[0][0] = 1.0 - x[0]*x[0];
	rm[1][1] = 1.0 - x[1]*x[1];
	rm[2][2] = 1.0 - x[2]*x[2];

	rm[0][1] = rm[1][0] = -x[0]*x[1];
	rm[0][2] = rm[2][0] = -x[0]*x[2];
	rm[1][2] = rm[2][1] = -x[1]*x[2];
}

void
atAberration(
	double mjd,    /* input: modified Julian Day */
	AtVect x0,     /* input: vector in equatorial coordiante at mjd */
	AtVect x)      /* output: vector in equatorial coordiante at mjd */
{
	AtVect xPrecess, aberr, xAberr, earth;
	AtVect posSun, posSunPrecess, posSunEcliptic;
	AtRotMat rm;
	double lambda, cos_lmd, cos_eps, sin_eps;
	AtRotMat precession_rm, precession_inv_rm;

	atPrecessRMJ2000(mjd, precession_inv_rm);
	ATInvRotMat(precession_inv_rm, precession_rm);

/*	atPrecession(MJD_J2000, x0, mjd, xPrecess);*/
    ATRotVect(precession_rm, x0, xPrecess);

	atSun(mjd, posSun);
/*	atPrecession(MJD_J2000, posSun, mjd, posSunPrecess);*/
    ATRotVect(precession_rm, posSun, posSunPrecess);

	cos_eps = COS_EPSILON;	/*cos(EPSILON);*/
	sin_eps = SIN_EPSILON;	/*sin(EPSILON);*/
	posSunEcliptic[0] =  posSunPrecess[0];
	posSunEcliptic[1] =  posSunPrecess[1]*cos_eps + posSunPrecess[2]*sin_eps;
	posSunEcliptic[2] = -posSunPrecess[1]*sin_eps + posSunPrecess[2]*cos_eps;

	lambda = atan2(posSunEcliptic[1], posSunEcliptic[0]);

	cos_lmd = cos(lambda);
	earth[0] =  sin(lambda);
	earth[1] = - cos_lmd * cos_eps;
	earth[2] = - cos_lmd * sin_eps;

	setAberrMat(xPrecess, rm);

	ATRotVect(rm, earth, aberr);
	aberr[0] *= KAPPA;
	aberr[1] *= KAPPA;
	aberr[2] *= KAPPA;
	ATAddVect(xPrecess, aberr, xAberr);

/*	atPrecession(mjd, xAberr, MJD_J2000, x);*/
    ATRotVect(precession_inv_rm, xAberr, x);
}

void
atInvAberration(
	double mjd,     /* input: modified Julian Day */
	AtVect x0,      /* input: vector in equatorial coordiante at mjd */
	AtVect x)       /* output: vector in equatorial coordiante at mjd */
{
	AtVect xPrecess, aberr, xAberr, earth;
	AtVect posSun, posSunPrecess, posSunEcliptic;
	AtRotMat rm;
	double lambda, cos_lmd, cos_eps, sin_eps;
	AtRotMat precession_rm, precession_inv_rm;

	atPrecessRMJ2000(mjd, precession_inv_rm);
	ATInvRotMat(precession_inv_rm, precession_rm);

/*	atPrecession(MJD_J2000, x0, mjd, xPrecess);*/
    ATRotVect(precession_rm, x0, xPrecess);

	atSun(mjd, posSun);
/*	atPrecession(MJD_J2000, posSun, mjd, posSunPrecess);*/
    ATRotVect(precession_rm, posSun, posSunPrecess);

	cos_eps = COS_EPSILON;	/*cos(EPSILON);*/
	sin_eps = SIN_EPSILON;	/*sin(EPSILON);*/
	posSunEcliptic[0] =  posSunPrecess[0];
	posSunEcliptic[1] =  posSunPrecess[1]*cos_eps + posSunPrecess[2]*sin_eps;
	posSunEcliptic[2] = -posSunPrecess[1]*sin_eps + posSunPrecess[2]*cos_eps;

	lambda = atan2(posSunEcliptic[1], posSunEcliptic[0]);

	cos_lmd = cos(lambda);
	earth[0] =  sin(lambda);
	earth[1] = - cos_lmd * cos_eps;
	earth[2] = - cos_lmd * sin_eps;

	setAberrMat(xPrecess, rm);

	ATRotVect(rm, earth, aberr);
	aberr[0] *= KAPPA;
	aberr[1] *= KAPPA;
	aberr[2] *= KAPPA;
/*	ATInvVect(aberr, invaberr);
	ATAddVect(xPrecess, invaberr, xAberr);*/
	ATSubVect(xPrecess, aberr, xAberr);

/*	atPrecession(mjd, xAberr, MJD_J2000, x);*/
    ATRotVect(precession_inv_rm, xAberr, x);
}

#ifdef TEST
/*
 * test main routine
 */
int
main(void)
{
	AtPolarVect60 star;
	AtVect posSun, posSun2, posStar, posStar2;
	double mjd;
	AtTime attime;
	double lambda;
	double r, alpha, delta;
	AtVect earth, aberr;
	AtRotMat rm;
	int i;

	printf("time = ");
	scanf("%d %d %d %d %d %d %f",
		&attime.yr, &attime.mo, &attime.dy, &attime.hr,
		&attime.mn, &attime.sc, &attime.ms);
	atMJulian(&attime, &mjd);
	printf("mjd = %lf\n", mjd);

	printf("position (alpha) = ");
	scanf("%lf", &alpha);
	printf("position (delta) = ");
	scanf("%lf", &delta);
	atPolDegToVect(1.0, alpha, delta, posStar);
	printf("Direction of Star -> (%lf, %lf, %lf)\n",
		posStar[0], posStar[1], posStar[2]);

	atPrecession(mjd, posStar, MJD_J2000, posStar2);
	atAberration(mjd, posStar2, posStar);
	atPrecession(MJD_J2000, posStar, mjd, posStar2);

	atVectToPolDeg(posStar2, &r, &alpha, &delta);
	printf("alpha = %lf delta = %lf\n", alpha, delta);

	return 0;
}
#endif
