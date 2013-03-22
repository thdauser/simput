#include "atFunctions.h"
#include "atError.h"

/*
 * product of two quaternion: q2 = q0 q1
 * 	(in matrix representation: rm(q2) = rm(q1) rm(q0).
 *	note the inverse order of quaternion as compared to matrix )
 * no checking for consistency for unitarity is done.
 *						n.kawai 93.01.13
 */
int
atQuatProd(
	AtQuat q0,	/* input: quaternion to be multiplied*/
	AtQuat q1,	/* input: quaternion to multiply*/
	AtQuat q2)	/* output: product */
{
	q2[0] =  q1[3]*q0[0] + q1[2]*q0[1] - q1[1]*q0[2] + q1[0]*q0[3];
	q2[1] = -q1[2]*q0[0] + q1[3]*q0[1] + q1[0]*q0[2] + q1[1]*q0[3];
	q2[2] =  q1[1]*q0[0] - q1[0]*q0[1] + q1[3]*q0[2] + q1[2]*q0[3];
	q2[3] = -q1[0]*q0[0] - q1[1]*q0[1] - q1[2]*q0[2] + q1[3]*q0[3];
    return NORMAL_END;
}
