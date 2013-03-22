#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"

double
atParseDecToRadian(char *expression)
{
	char *p, *q;
	AtDeclination dec;
	dec.sign = 1;
	dec.deg = dec.min = 0;
	dec.sec = 0.0;
	p = q = expression;
	if ( '-' == *p ) dec.sign = -1;
	while ( *q++ ) {
		if ( 'd' == *q || 'D' == *q ) {
			dec.deg = abs(atoi(p));
			p = q + 1;
		} else if ( 'm' == *q || 'M' == *q ) {
			dec.min = atoi(p);
			p = q + 1;
		}
	}
	dec.sec = atof(p);
	if ( p == expression ) return DEG2RAD * dec.sec;
	return atDecToRadian(dec);
}
