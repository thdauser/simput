#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * CALC MATRIX FOR Rotation of coordinate with AXIS AND ROLL ANGLE.
 */
int
atSetRotMat(
	AtVect axis,	/* input: rotation axis of coord, should be non zero*/
	double roll,	/* input: roll angle around axis (radian) */
	AtRotMat rm)	/* output: rotation matrix */
{
    AtVect normAxis;
    double cosRoll, c1, sinRoll;
    int i, i1, i2, cond;

    cond = atNormVect(axis, normAxis);
    if (cond == NORMAL_END) {
		cosRoll = cos(roll);
		c1 = 1-cosRoll;
		sinRoll = sin(roll);
		for (i=0; i<3; i++) {
			i1 = (i+1)%3;
			i2 = (i+2)%3;
			rm[i][i]  = cosRoll + normAxis[i]*normAxis[i] *c1;
			rm[i1][i] = normAxis[i]*normAxis[i1]*c1 - normAxis[i2]*sinRoll;
			rm[i2][i] = normAxis[i]*normAxis[i2]*c1 + normAxis[i1]*sinRoll;
		}
    }
    return cond;
}
