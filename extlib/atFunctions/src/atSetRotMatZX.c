#include "atFunctions.h"

/*
 * ROTATION MATRIX defined by New Z-axis and a vector in (+)X-Z half plane.
 */
int
atSetRotMatZX(
	AtVect zAxis,	/* input: vector defining new z-axis */
	AtVect xAxis,	/* input: vector in new +X-Z half plane */
	AtRotMat rm)	/* output: rotation matrix */
{
    int code;
    AtVect x, y, z, yAxis;

    ATVectProd(zAxis, xAxis, yAxis);
    if ((code = atNormVect(zAxis, z))!=0) return code;
    if ((code = atNormVect(yAxis, y))!=0) return code;
    ATVectProd(y, z, x);

	rm[0][0] = x[0];
	rm[0][1] = x[1];
	rm[0][2] = x[2];

	rm[1][0] = y[0];
	rm[1][1] = y[1];
	rm[1][2] = y[2];

	rm[2][0] = z[0];
	rm[2][1] = z[1];
	rm[2][2] = z[2];

    return 0;
}
