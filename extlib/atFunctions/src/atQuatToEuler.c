#include "atFunctions.h"
#include "atError.h"

/*
 * convert quaternion parameter set to Euler angles
 */
int
atQuatToEuler(
	AtQuat q,			/* input: quaternion */
	AtEulerAng *ea		/* output: z-y-z Euler Angle (radian) */
)
{
	int istat;
	AtRotMat rm;

	atQuatToRM(q, rm);
	istat = atRMToEuler(rm, ea);

	return istat;
}
