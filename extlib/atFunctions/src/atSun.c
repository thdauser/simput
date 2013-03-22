/************************************************************************
  atSun.c   calculate position of the Sun

  2007/04/07 Y.ISHISAKI
	- update atSun() for faster calculation
************************************************************************/

#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * calc position of the sun in J2000 coordinate
 */
int
atSun(double mjd,	/* input: time in MJD */
	  AtVect pos)	/* output: vector to the sun in A.U.*/
{
	static int first_time = 1;
	static AtRotMat rm1950_2000;

    double l, m, r, t;
	double sin_2m, sin_l;
    AtVect x;

	if ( first_time ) {
		atPrecessRMJ2000(MJD_B1950, rm1950_2000);
		first_time = 0;
	}

    t = mjd - 4.5e4;
    m = ( fmod(t * .985600267, 360.0) + 27.26464 ) * DEG2RAD;
	sin_2m = sin(2*m);

    l = ( fmod(t * .985609104, 360.0) - 50.55138
	      + sin(m) * 1.91553 + sin_2m * .0201 ) * DEG2RAD;
	sin_l = sin(l);

    r = 1.00014 - cos(m) * .01672 - sin_2m * 1.4e-4;

    x[0] = r * cos(l);
    x[1] = r * .91744 * sin_l;
    x[2] = r * .39788 * sin_l;

/*	atPrecession(MJD_B1950, x, MJD_J2000, pos);*/
    ATRotVect(rm1950_2000, x, pos);

    return NORMAL_END;
}
