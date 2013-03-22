#include "atFunctions.h"

/*
 * Converting Declination in deg.min.sec to radian
 *			  	ver 1.0  92/07/01  ay
 *             			ver 1.5  93/01/26  n.kawai
 */

double atDecToRadian(	        /* out: declination in radian  */
	AtDeclination dec) 	/* in:  declination in ddmmss */
{                               /*        ver 1.5 93/01/26 n.kawai */
    double  x;
    x = dec.sign * (dec.deg + dec.min/60. + dec.sec/3600.);
    return (x* DEG2RAD);
}
