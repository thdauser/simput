#include "atFunctions.h"
#include "atError.h"

/*
 * ROTATE A VECTOR WITH ROTATION MATRIX.
 */
int
atRotVect(
	AtRotMat rm,	/* input: rotation matrix */
	AtVect x,		/* input: vector */
	AtVect y)		/* output: vector */
{
	return ATRotVect(rm, x, y);
}
