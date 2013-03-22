#include "atFunctions.h"
#include "atError.h"

/*
 * convert Euler angles to quaternion parameter set
 */
int
atEulerToQuat(
	AtEulerAng *ea,		/* output: z-y-z Euler Angle (radian) */
	AtQuat q			/* input: quaternion */
)
{
	AtRotMat rm;

	atEulerToRM(ea, rm);
	atRMToQuat(rm, q);

	return NORMAL_END;
}
