#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"

void
atDegToDec(double deg, AtDeclination *dec)
{
	double i_min;
	if ( deg < 0.0 ) {
		dec->sign = -1;
		deg = -deg;
    } else {
    	dec->sign = 1.;
    }
	dec->sec = 60 * modf(60*deg+EPS, &i_min);
	dec->deg = i_min / 60 + EPS;
	dec->min = i_min - 60*dec->deg + EPS;
}

void
atRadianToDec(double radian, AtDeclination *dec)
{
	atDegToDec(radian*RAD2DEG, dec);
}
