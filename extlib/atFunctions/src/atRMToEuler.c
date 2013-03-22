#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * Finding Euler angles from a given Rotation-matrix
 *                             ver 1.0  92/07/01  ay
 *	case for theta=0 or 180 deg	ver 1.1 93.1.11 n.kawai
 *	revised with correct definition	ver 1.2 93.1.12 n.kawai
 */
int
atRMToEuler(
        AtRotMat rm,            /* input: rotation matrix */
        AtEulerAng *ea)         /* output: z-y-z Euler Angle (radian) */
{
    if ( rm[2][2] > 1.-EPS ) {	/* only rotation about z */
		ea->theta = 0;
		ea->phi = atan2(rm[0][1], rm[0][0]);
		ea->psi = 0;
    } else if ( rm[2][2] < -1.+EPS ) { /* flip w/ y and rot about z */
		ea->theta = PI;
		ea->phi = atan2(-rm[0][1], -rm[0][0]);
		ea->psi = 0;
    } else {
		ea->theta = acos( rm[2][2] );

		if ( fabs(rm[2][0]) > EPS || fabs(rm[2][1]) > EPS ) {
			ea->phi = atan2(rm[2][1], rm[2][0]);
		} else {
			ea->phi = 0;
			return INCONSISTENT_RM;
		}

		if ( fabs(rm[0][2]) > EPS || fabs(rm[1][2]) > EPS ) {
			ea->psi = atan2(rm[1][2], -rm[0][2]);
		} else {
			ea->psi = 0;
			return INCONSISTENT_RM;
		}
    }
    if ( ea->phi < 0.) {
		ea->phi += TWO_PI;
	}
    if ( ea->psi < 0.) {
		ea->psi += TWO_PI;
	}
    return NORMAL_END;
}
