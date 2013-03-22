#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * calc Greenwich sidereal time (KOSEJI)

   v2.4 March 31, 2005,  S. Yamauchi
		Reference: Kenneth R. Lang, "Astrophysical Formulae"
              vol 2, Sec. 5.3.7, formula (5.268)

   2005/08/06 v2.5	Y.ISHISAKI
   		modified to use floor() instead of (long)mjd, for negative value
		remove static for local variables
 */
int
atSidereal(
	double mjd,		/* input: time in MJD */
	double *gsttod)	/* output: Greenwich sidereal time (radian)
							   at mjd true of date */
{
/* Local variables */
	double d, m, x, a0, a1, a2, a3, g;

	a3 = -6.2e-6 / 3600.0;
	a2 = 0.093104 / 3600.0;
	a1 = 8640184.812866 / 3600.0 ;
	a0 = 24110.54841 / 3600.0;

	d = floor(mjd);
	x = (d - 51544.5) / 36525.;
	m = (mjd - d) * 24. * 60.;
	g = (((a3 * x + a2) * x + a1) * x + a0)*15.0 + m * .25068447;

	while ( 360.0 <= g ) {
		g -= 360.0;
	}

	while ( g < 0.0 ) {
		g += 360.0;
	}

	*gsttod = g * DEG2RAD;

	return NORMAL_END;
}
