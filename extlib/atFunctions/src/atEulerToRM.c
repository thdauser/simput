#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/* Finding a Rotation-matrix for given Euler angles

			Revised with correct definition  N.Kawai 93.01.12
 */

int atEulerToRM(
        AtEulerAng *ea,  /* input: z-y-z Euler Angle (radian) */
        AtRotMat rm)     /* output: rotation matrix */
{
	double  c[3], s[3];
	c[0]=cos( ea->phi );   c[1]=cos( ea->theta );   c[2]=cos( ea->psi );
	s[0]=sin( ea->phi );   s[1]=sin( ea->theta );   s[2]=sin( ea->psi );
	rm[0][0] =  c[2]*c[1]*c[0] - s[2]*s[0];
	rm[0][1] =  c[2]*c[1]*s[0] + s[2]*c[0];
	rm[0][2] =                 - c[2]*s[1];
	rm[1][0] = -s[2]*c[1]*c[0] - c[2]*s[0];
	rm[1][1] = -s[2]*c[1]*s[0] + c[2]*c[0];
	rm[1][2] =                   s[2]*s[1];
	rm[2][0] =  s[1]*c[0];
	rm[2][1] =  s[1]*s[0];
	rm[2][2] =  c[1];
	return NORMAL_END;
}
