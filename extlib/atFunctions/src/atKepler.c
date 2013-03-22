#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * solve Kepler equation (KEPLER)  g + e sin E = E
 */
int atKepler(
	double g,	/* input: mean anomaly */
	double eccent,	/* input: eccentricity */
	double *e)	/* output: eccentric anomaly */
{
    static double eps = 1e-15;
	static int imax = 50;

	int i;
	static double error, deltae, d__1;

	*e = g;
	if (g == 0.) return 0;

	for (i=0; i<imax; i++) {
		deltae = (g - *e + eccent * sin(*e)) / (1. - eccent * cos(*e));
		*e += deltae;
		error = (d__1 = deltae / *e, fabs(d__1));
		if (error < eps) return 0;
	}
	return NOT_CONVERGED;
}
