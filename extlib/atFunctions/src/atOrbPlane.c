#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * coordinate conversion from orbital plane to celestial (SATINA)
 */
int atOrbPlane(
	AtVect x,	/* input: vector in orbital plane */
	double ol,	/* input: Angle from Ascending Node */
	double ob, 	/* input: Right Ascention of the Ascending Node */
	double ai,	/* input: inclination of the orbit */
	AtVect y)	/* output: vector in celestial coordinate at mjd */
{
    double q1, q2, q3, q4, q5, q6, qq1, qq2;

    q1 = cos(ob);
    q2 = sin(ob);
    q3 = cos(ai);
    q4 = sin(ai);
    q5 = cos(ol);
    q6 = sin(ol);
    qq1 = x[0] * q5 - x[1] * q6;
    qq2 = x[0] * q6 + x[1] * q5;
    y[0] = q1 * qq1 - q2 * (qq2 * q3 - q4 * x[2]);
    y[1] = q2 * qq1 + q1 * (qq2 * q3 - q4 * x[2]);
    y[2] = q4 * qq2 + q3 * x[2];
    return (NORMAL_END);
}
