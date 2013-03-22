#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"

void
atDegToRA(double deg, AtRightAscension *ra)
{
	double i_min;
	ra->sec = 60 * modf(4*deg+EPS, &i_min) ;
	ra->hour = i_min / 60 + EPS;
	ra->min = i_min - 60*ra->hour + EPS;
}

void
atRadianToRA(double radian, AtRightAscension *ra)
{
	atDegToRA(radian*RAD2DEG, ra);
}
