/************************************************************************
  atInvRotMat.c		Making an Inverse Rotation Matrix

	1992/07/01 A.YOSHIDA	version 1.0

	2005/12/04 Y.ISHISAKI	version 2.7
		use fast macro
************************************************************************/

#include "atFunctions.h"
#include "atError.h"

/*
 * Making an Inverse Rotation Matrix
 */
int
atInvRotMat(
	    AtRotMat rm,	/* input: rotation matrix */
	    AtRotMat rm2)	/* output: inversed rotation matrix */
{
	return ATInvRotMat(rm, rm2);
}
