/************************************************************************
  atInvRotMat.c

	1993/01/13 N.KAWAI		version 1.0

	2005/12/04 Y.ISHISAKI	version 2.7
		use fast macro
************************************************************************/
#include "atFunctions.h"
#include "atError.h"

/*
 * product of two rotation matrices rm2 = rm1 rm0
 * no checking for consistency for unitarity is done.
 */
int
atRMProd(
	AtRotMat rm0,	/* input: rotation matrix to be multiplied*/
	AtRotMat rm1,	/* input: rotation matrix to multiply*/
	AtRotMat rm2)	/* output: product */
{
	return ATRMProd(rm0, rm1, rm2);
}
