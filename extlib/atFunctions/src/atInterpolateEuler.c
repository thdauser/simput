#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"
#include "atError.h"

/*
 * interpolation between two sets of Euler angles
 */
int
atInterpolateEuler(
	double t0, AtEulerAng *ea0,	/* input: q-parameter q0 at time t0 */
	double t1, AtEulerAng *ea1,	/* input: q-parameter q0 at time t0 */
	double t,  AtEulerAng *ea	/* input: time t
								   output: interpolated Euler angles */
)
{
	AtQuat q0, q1, q;

	atEulerToQuat(ea0, q0);
	atEulerToQuat(ea1, q1);
	atInterpolateQuat(t0, q0, t1, q0, t, q);
	atQuatToEuler(q, ea);

	return NORMAL_END;
}
