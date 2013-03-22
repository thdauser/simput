#include <stdio.h>
#include <stdlib.h>

/* Numerical Recipes in C, pp.96-97, slightly modified */
int
atSpline(float x[],float y[], int n, float yp1, float ypn, float y2[])
{
    int i, k;
    float p, qn, sig, un, *u;

    u = malloc( n * sizeof(float) );

    if (yp1 > 0.99e30) {
		y2[0] = u[0] = 0.0;
    } else {
		y2[0] = -0.5;
		u[0] = (3.0/ (x[1]-x[0]) )*( (y[1]-y[0]) / (x[1]-x[0]) -yp1);
	}
    for (i=1; i<n-1; i++) {
		sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
		p = sig*y2[i-1] + 2.0;
		y2[i] = (sig - 1.)/p;
		u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
		u[i] = (6.*u[i]/ (x[i+1]-x[i-1]) - sig*u[i-1])/p;
    }
    if (ypn>0.99e30) {
		qn = un = 0.0;
    } else {
		qn = 0.5;
		un = (3./ (x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
    }
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.);
	for (k=n-2; k>=0; k--) {
		y2[k] = y2[k] * y2[k+1] + u[k];
	}

    free(u);
    return 0;
}

int
atSplint(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
    int klo, khi, k;
    float h, b, a;

    klo = 0;
    khi = n-1;
    while (khi - klo > 1) {
		k = (khi+klo) >> 1;
		if (xa[k] > x) {
			khi = k;
		} else {
			klo = k;
		}
	}
	h = xa[khi]-xa[klo];

    if (h==0.) return -20;

	a = (xa[khi] - x) /h;
	b = (x - xa[klo]) /h;
	*y = a*ya[klo]+b*ya[khi] +
		((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.;

    return 0;
}

/*
   Two-dimensional Spline interpolation
   Numerical Recipes in C, pp.109-110, slightly modified

   Given a tabulated funcion ya[0..m-1][0..n-1], and tabulated
   independent variables x1a[0..m-1] and x2a[0..n-1], this routine
   constructs one-dimensional natural cubic splines of the rows of ya
   and returns the second-derivativbes in the array
   y2a[0..m-1][0..n-1]. */
int
atSplie2(float x1a[], float x2a[], float ya[], int m, int n, float y2a[])
{
    int j;

    for (j=0; j<m; j++) {
		atSpline(x2a, ya+j*n, n, 1.e30, 1.e30, y2a+j*n);
	}

	return 0;
}


/*
   Given x1a, x2a, ya, m, n as described in "splie2", and y2a as
   produced by that routine; and given a desired interpolating point
   x1, x2; this routine returns an interpolated function value "y" by
   bicubic spline interpolation */
int
atSplin2(float x1a[], float x2a[], float ya[], float y2a[], int m, int n,
	float x1, float x2, float *y)
{
    int j, code;
    float *ytmp, *yytmp;

    ytmp = malloc( m * sizeof(float) );
    yytmp = malloc( m * sizeof(float) );

    for (j=0; j<m; j++) {
		atSplint(x2a, ya+j*n, y2a+j*n, n, x2, &yytmp[j]);
	}
	atSpline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
    code = atSplint(x1a, yytmp, ytmp, m, x1, y);

    free(yytmp);
    free(ytmp);
    return code;
}
