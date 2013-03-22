#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"

double
atParseRAToRadian(char *expression)
{
	char *p, *q;
	AtRightAscension ra;
	ra.hour = ra.min = 0;
	ra.sec = 0.0;
	p = q = expression;
	while ( *q++ ) {
		if ( 'h' == *q || 'H' == *q ) {
			ra.hour = atoi(p);
			p = q + 1;
		} else if ( 'm' == *q || 'M' == *q ) {
			ra.min = atoi(p);
			p = q + 1;
		}
	}
	ra.sec = atof(p);
	if ( p == expression ) return DEG2RAD * ra.sec;
	return atRAToRadian(ra);
}
