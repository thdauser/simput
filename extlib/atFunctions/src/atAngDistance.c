/************************************************************************
  atAngDistance.c	calc angular distance (RIKAKU)

	1992/07/01 A.YOSHIDA	version 1.0

	2005/12/04 Y.ISHISAKI	version 2.7
		modified to use atan2(), [atAngDistance3]
		because acos() does not have enough precision in small angle
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "atFunctions.h"
#include "atError.h"

#define atAngDistance3	atAngDistance

#ifdef atAngDistance1
int                 /* Finding angular distance                 */
atAngDistance1(     /*                   ver 1.0   92/07/01  ay */
        AtVect x,   /* input: */
        AtVect y,   /* input: */
        double *r)  /* output: angular distance between x and y in radian */
{
  double  d1;
  d1 = ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] ) / atNorm(x) / atNorm(y) ;
  if ( d1 > 1.0 ) { d1 = 1.0; }
  if ( d1 < -1.0 ) { d1 = -1.0; }
  *r = acos( d1 );
  printf("AngDistance: d1=%.8f, acos=%.12f\n", d1, *r);
  if ( *r < 0.0 || *r > PI ) { return (-1); }
  return (NORMAL_END);
}
#endif

#ifdef atAngDistance2
int                 /* Finding angular distance */
atAngDistance2(
        AtVect x,   /* input: */
        AtVect y,   /* input: */
        double *r)  /* output: angular distance between x and y in radian */
{
	AtVect z;
	double lx, ly, lz, d1;
	lx = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	ly = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	z[0] = x[0]/lx - y[0]/ly;
	z[1] = x[1]/lx - y[1]/ly;
	z[2] = x[2]/lx - y[2]/ly;
	lz = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
	if ( 1.5 < lz ) {
		d1 = ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] ) / lx / lz;
		if ( d1 < -1.0 ) d1 = -1.0;
/*		printf("AngDistance2: d1=%.8f, acos=%.12f\n", d1, *r); */
		*r = acos(d1);
	} else {
		*r = 2 * asin(lz/2);
/*		printf("AngDistance2: lz=%.8f, asin=%.12f\n", lz, *r); */
	}
	return 0;
}
#endif

#ifdef atAngDistance3
int                 /* Finding angular distance */
atAngDistance3(
        AtVect x,   /* input: */
        AtVect y,   /* input: */
        double *r)  /* output: angular distance between x and y in radian */
{
#define DOT(v1,v2)	(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
#define NORM(v)		sqrt(DOT(v,v))
	AtVect nvx, nvy, xy, z;
	double lx, ly, lxy, lz, acos_2, asin_2;

	lx = ATNorm(x);
	ly = ATNorm(y);
	if ( 0.0 == lx || 0.0 == ly ) {
		*r = 0.0;
		return NULL_VECTOR;
	}

	nvx[0] = x[0] / lx;
	nvx[1] = x[1] / lx;
	nvx[2] = x[2] / lx;

	nvy[0] = y[0] / ly;
	nvy[1] = y[1] / ly;
	nvy[2] = y[2] / ly;

	xy[0] = nvx[0] + nvy[0];
	xy[1] = nvx[1] + nvy[1];
	xy[2] = nvx[2] + nvy[2];

	lxy = ATNorm(xy);
	if ( 0.0 == lxy ) {					/* x, y is opposite direction */
		*r = PI;
		return 0;
	}
	acos_2 = ATScalProd(nvx, xy) / lxy;

	z[0] = nvx[0] - nvy[0];
	z[1] = nvx[1] - nvy[1];
	z[2] = nvx[2] - nvy[2];
	lz = ATNorm(z);
	asin_2 = lz / 2;

	*r = 2 * atan2(asin_2, acos_2);		/* since 0 <= lz, 0 <= *r <= PI */

#undef NORM
#undef DOT

	return 0;
}
#endif
