/*************************************************************************
  atSpline.c

	2005/12/04 Y.ISHISAKI	v2.7
		add atSplineP, atSplintP, atSplineD, atSplintD,
	    	atSpline2D,atSplint2D, atSplint2P
**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atFunctions.h"
#include "atSpline.h"
#include "atError.h"

/* uncomment following line for debug message with 1st derivatives */
/* #define DEBUG_Y1 */

#ifdef DEBUG_Y1
/* calculate 1st derivative */
static int
atSpline_y1(double xa[], double ya[], double y2a[], int n, double y1a[])
{
	int i;
	double a, b, h;

	a = 1.0;
	b = 0.0;
	h = 1.0;
	for (i = 0; i < n-1; i++) {
		h = xa[i+1] - xa[i];
		y1a[i] = (ya[i+1] - ya[i]) / h +
			( (3*b*b-1)*y2a[i+1] - (3*a*a-1)*y2a[i] ) * h / 6.0;
	}
	a = 0.0;
	b = 1.0;
	y1a[n-1] = (ya[i] - ya[i-1]) / h +
			( (3*b*b-1)*y2a[i] - (3*a*a-1)*y2a[i-1] ) * h / 6.0;

	return 0;
}
#endif

/* Spline interpolation for periodic data, y[0] must be y[n-1] */
int
atSplineP(double x[], double y[], int n, double z[], double work[/* 3*n */])
{
	int i;
	double t, *h, *d, *w;

	h = work;
	d = &work[n];
	w = &work[2*n];

	if ( y[0] != y[n-1] ) {
		return -1;
	}

	if ( n < 3 ) {
		for (i = 0; i < n; i++) {
			z[i] = 0.0;
		}
		return 0;
	}

	h[0] = x[1] - x[0];
	if ( h[0] <= 0.0 ) return -1;
	w[0] = (y[1] - y[0]) / h[0];
	for (i = 1; i < n-1; i++) {
		h[i] = x[i+1] - x[i];
		if ( h[i] <= 0.0 ) return -1;
		w[i] = ( y[i+1] - y[i] ) / h[i];
		d[i] = 2 * ( x[i+1] - x[i-1] );
		z[i] = w[i] - w[i-1];
	}
	w[n-1] = w[0];
	d[n-1] = 2 * ( h[n-2] + h[0] );
	z[n-1] = w[n-1] - w[n-2];

	w[1] = h[0];
	for (i = 2; i < n-2; i++) {
		w[i] = 0;
	}
	w[n-2] = h[n-2];
	w[n-1] = d[n-1];

	for (i = 1; i < n-1; i++) {
		t = h[i] / d[i];
		z[i+1] = z[i+1] - z[i] * t;
		d[i+1] = d[i+1] - h[i] * t;
		w[i+1] = w[i+1] - w[i] * t;
	}

	w[0] = w[n-1];
	z[0] = z[n-1];
	for (i = n - 3; 0 <= i; i--) {
		t = h[i] / d[i+1];
		z[i] = z[i] - z[i+1] * t;
		w[i] = w[i] - w[i+1] * t;
	}

	t = z[0] / w[0];
	z[0] = 6*t;
	for (i = 1; i < n-1; i++) {
		z[i] = 6 * ( z[i] - w[i]*t ) / d[i];
	}
	z[n-1] = 6*t;

	return 0;
}

int
atSplintP(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo, khi, k;
	double period, h, b, a;

	period = xa[n-1] - xa[0];
	while ( xa[n-1] < x ) x -= period;
	while ( x < xa[0] ) x += period;

	klo = 0;
	khi = n-1;
	while ( 1 < khi - klo ) {
		k = (khi+klo) >> 1;
		if ( x < xa[k] ) {
			khi = k;
		} else {
			klo = k;
		}
	}

	h = xa[khi]-xa[klo];
	if ( 0.0 == h ) return -1;

	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;

	*y = a*ya[klo] + b*ya[khi] +
		((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.;

	return 0;
}

/* Numerical Recipes in C, pp.96-97, slightly modified */
int
atSplineD(double x[], double y[],
	int n, double yp1, double ypn, double y2[], double work[/* n */])
{
	int i;
	double h0, h1, hh, p, qn, sig, un, *u;

	if ( n < 3 ) {
		for (i = 0; i < n; i++) {
			y2[i] = 0.0;
		}
		return 0;
	}

/*	u = malloc((unsigned) n*sizeof(*y2));*/
	u = work;

	h0 = x[1] - x[0];
	if ( h0 <= 0.0 ) return -1;			/* must be positive value */

	if ( 0.99e30 < yp1 ) {
		y2[0] = u[0] = 0.0;
	} else {
		y2[0] = -0.5;
		u[0] = (3.0 / h0) * ( (y[1]-y[0]) / h0 - yp1 );
	}

	for (i = 1; i < n-1; i++) {
		h1 = x[i+1] - x[i];
		if ( h1 <= 0.0 ) return -1;		/* must be positive value */
		hh = x[i+1] - x[i-1];
		sig = h1 / hh;
		p = sig*y2[i-1] + 2.0;
		y2[i] = (sig - 1.) / p;
		u[i] = (y[i+1]-y[i]) / h1 - (y[i]-y[i-1]) / h0;
		u[i] = ( 6.*u[i]/hh - sig*u[i-1] ) / p;
		h0 = h1;
	}

	if ( 0.99e30 < ypn ) {
		qn = un = 0.0;
	} else {
		qn = 0.5;
		un = (3.0 / h0) * ( ypn - (y[n-1]-y[n-2]) / h0 );
	}

	y2[n-1] = ( un - qn*u[n-2] ) / ( qn*y2[n-2]+1.0 );
	for (i = n-2; 0 <= i; i--) {
		y2[i] = y2[i] * y2[i+1] + u[i];
	}

	return 0;
}

int
atSplintD(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo, khi, k;
	double h, b, a;

	klo = 0;
	khi = n-1;
	while ( 1 < khi - klo ) {
		k = (khi+klo) >> 1;
		if ( x < xa[k] ) {
			khi = k;
		} else {
			klo = k;
		}
	}

	h = xa[khi]-xa[klo];
	if ( 0.0 == h ) return -1;

	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;

	*y = a*ya[klo] + b*ya[khi] +
		((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.;

	return 0;
}

/* Two-dimensional Spline interpolation
   Numerical Recipes in C, pp.109-110, slightly modified

   Given a tabulated funcion ya[0..m-1][0..n-1], and tabulated independent
   variables x1a[0..m-1] and x2a[0..n-1], this routine constructs
   one-dimensional natural cubic splines of the rows of ya and returns
   the second-derivativbes in the array y2a[0..m-1][0..n-1]. */
int
atSpline2D(double x1a[], double x2a[], double ya[],
	int m, int n, double y2a[], double work[/* n */])
{
	int j, code;

	for (j = 0; j < m; j++) {
		code = atSplineD(x2a, ya+j*n, n, 1.e30, 1.e30, y2a+j*n, work);
		if ( code ) return code;
	}

	return 0;
}

/* Given x1a, x2a, ya, m, n as described in "splie2", and y2a as
   produced by that routine; and given a desired interpolating point
   x1, x2; this routine returns an interpolated function value "y" by
   bicubic spline interpolation */
int
atSplint2D(double x1a[], double x2a[], double ya[], double y2a[],
	int m, int n, double x1, double x2, double *y, double work[/* 3*m */])
{
	int j, code;
	double *ytmp, *y2tmp;

/*	ytmp = malloc((unsigned) m*sizeof(double));
	y2tmp = malloc((unsigned) m*sizeof(double));
*/
	ytmp = &work[m];
	y2tmp = &work[2*m];

	for (j = 0; j < m; j++) {
		code = atSplintD(x2a, ya+j*n, y2a+j*n, n, x2, &ytmp[j]);
		if ( code ) return code;
	}
	atSplineD(x1a, ytmp, m, 1.0e30, 1.0e30, y2tmp, work);
	code = atSplintD(x1a, ytmp, y2tmp, m, x1, y);

#ifdef DEBUG_Y1
	atSpline_y1(x1a, ytmp, y2tmp, m, work);
	for (j = 0; j < m; j++) {
		printf("%5.1f %6.3f %+.8f %+.8f\n", x1a[j],ytmp[j],work[j],y2tmp[j]);
	}
#endif

/*	free(y2tmp);
	free(ytmp);*/

	return code;
}

/* Two-dimensional Spline interpolation, with periodic x-axis */
int
atSplint2P(double x1a[], double x2a[], double ya[], double y2a[],
	int m, int n, double x1, double x2, double *y, double work[/* 5*m */])
{
	int j, code;
	double *ytmp, *y2tmp;

	ytmp = &work[3*m];
	y2tmp = &work[4*m];

	for (j = 0; j < m; j++) {
		code = atSplintD(x2a, ya+j*n, y2a+j*n, n, x2, &ytmp[j]);
		if ( code ) return code;
	}
	atSplineP(x1a, ytmp, m, y2tmp, work);
	code = atSplintP(x1a, ytmp, y2tmp, m, x1, y);

#ifdef DEBUG_Y1
	atSpline_y1(x1a, ytmp, y2tmp, m, work);
	for (j = 0; j < m; j++) {
		printf("%5.1f %6.3f %+.8f %+.8f\n", x1a[j],ytmp[j],work[j],y2tmp[j]);
	}
#endif

	return code;
}
