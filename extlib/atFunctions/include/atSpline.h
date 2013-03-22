/*************************************************************************
  atSpline.h

  2005/12/04 Y.ISHISAKI	v2.7
	add atSplineP, atSplintP, atSplineD, atSplintD,
	    atSpline2D,atSplint2D, atSplint2P
**************************************************************************/

#ifdef __cplusplus
extern "C"
{
#endif

/* Spline interpolation for non-periodic data.
   y2[0] = y2[n-1] = 0.0  if ( 0.99e30 < yp1 ) .and. ( 0.99e30 < ypn ),
   otherwise, first derivatives of y1[0] = yp1, y1[n-1] = ypn.
   atSplineD() must be called before atSplintD() */
int atSplineD(
	double x[/* n */],	/* input: x values, must be increase order */
	double y[/* n */],	/* input: y values */
	int n,			/* input: dimension of data */
	double yp1,		/* input: 1st derivative at x[0] if < 1e30 */
	double ypn,		/* input: 1st derivative at x[n-1] if < 1e30 */
	double y2[/* n */],	/* output: 2nd derivatives, for atSplintD() */
	double work[/* n */]	/* input/output: required work area */
);

int atSplintD(
	double xa[/* n */],	/* input: x values, must be increase order */
	double ya[/* n */],	/* input: y values */
	double y2a[/* n */],	/* input: 2nd derivatives, by atSplineD() */
	int n,			/* input: dimension of data */
	double x,		/* input: where to interpolate */
	double *y		/* output: interpolated value */
);

/* Spline interpolation for periodic data, y[0] must be y[n-1].
   atSplineP() must be called before atSplintP() */
int atSplineP(
	double x[/* n */],	/* input: x values, must be increase order */
	double y[/* n */],	/* input: y values */
	int n,			/* input: dimension of data */
	double y2[/* n */],	/* output: 2nd derivatives, for atSplintP() */
	double work[/* 3*n */]	/* input/output: required work area */
);

int atSplintP(
	double xa[/* n */],	/* input: x values, must be increase order */
	double ya[/* n */],	/* input: y values */
	double y2a[/* n */],	/* input: 2nd derivatives, by atSplineP() */
	int n,			/* input: dimension of data */
	double x,		/* input: where to interpolate */
	double *y		/* output: interpolated value */
);

/* Two-dimensional Spline interpolation
   Numerical Recipes in C, pp.109-110, slightly modified.
   Given a tabulated funcion ya[0..m-1][0..n-1], and tabulated independent
   variables x1a[0..m-1] and x2a[0..n-1], this routine constructs
   one-dimensional natural cubic splines of the rows of ya and returns
   the second-derivativbes in the array y2a[0..m-1][0..n-1].
   atSpline2D() must be called before atSplint2D() or atSplint2P() */
int atSpline2D(
	double x1a[/* m */],	/* input: x values, must be increase order */
	double x2a[/* n */],	/* input: y values, must be increase order */
	double ya[/* m * n */],	/* input: z values */
	int m,			/* input: dimension of x1a[] */
	int n,			/* input: dimension of x2a[] */
	double y2a[/* m * n */],/* output: 2nd derivatives, for atSplint2D() */
	double work[/* n */]	/* input/output: required work area */
);

/* Given x1a, x2a, ya, m, n as described in "splie2", and y2a as
   produced by that routine; and given a desired interpolating point
   x1, x2; this routine returns an interpolated function value "y" by
   bicubic spline interpolation */
int atSplint2D(
	double x1a[/* m */],	/* input: x values, must be increase order */
	double x2a[/* n */],	/* input: y values, must be increase order */
	double ya[/* m * n */],	/* input: z values */
	double y2a[/* m * n */],/* input: 2nd derivatives, by atSpline2D() */
	int m,			/* input: dimension of x1a[] */
	int n,			/* input: dimension of x2a[] */
	double x1,		/* input: where to interpolate in 1st axis */
	double x2,		/* input: where to interpolate in 2nd axis */
	double *y,		/* output: interpolated value */
	double work[/* 3*m */]	/* input/output: required work area */
);

/* Two-dimensional Spline interpolation, with periodic x-axis */
int atSplint2P(
	double x1a[],		/* input: x values, must be increase order */
	double x2a[],		/* input: y values, must be increase order */
	double ya[],		/* input: z values */
	double y2a[],		/* input: 2nd derivatives, by atSpline2D() */
	int m,			/* input: dimension of x1a[] */
	int n,			/* input: dimension of x2a[] */
	double x1,		/* input: where to interpolate in 1st axis */
	double x2,		/* input: where to interpolate in 2nd axis */
	double *y,		/* output: interpolated value */
	double work[/* 5*m */]	/* input/output: required work area */
);

/***********************************************************************
   The following functions are obsolete (before atFunctions-2.7)
   and not recommended to use.  They uses malloc() internal.
***********************************************************************/

/* Numerical Recipes in C, pp.96-97, slightly modified
   atSpline() must be called before atSplint() */
int atSpline(float x[],float y[], int n, float yp1, float ypn, float y2[]);

int atSplint(float xa[], float ya[], float y2a[], int n, float x, float *y);

/* Two-dimensional Spline interpolation
   Numerical Recipes in C, pp.109-110, slightly modified

   Given a tabulated funcion ya[0..m-1][0..n-1], and tabulated independent
   variables x1a[0..m-1] and x2a[0..n-1], this routine constructs
   one-dimensional natural cubic splines of the rows of ya and returns
   the second-derivativbes in the array y2a[0..m-1][0..n-1].

   atSplie2() must be called before atSplin2() */
int atSplie2(float x1a[], float x2a[], float ya[], int m, int n, float y2a[]);

/* Given x1a, x2a, ya, m, n as described in "splie2", and y2a as
   produced by that routine; and given a desired interpolating point
   x1, x2; this routine returns an interpolated function value "y" by
   bicubic spline interpolation */
int atSplin2(float x1a[], float x2a[], float ya[], float y2a[], int m, int n,
	float x1, float x2, float *y);

#ifdef __cplusplus
}
#endif
