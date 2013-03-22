#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * Checking consistency (unitarity) of a rotation matrix
 *
 * if the matrix is not unitary, matrix is modified to be unitary
 * using quaternion conversion algorithm with a "INCONSISTENT_RM" error code.
 *
 * If there is a rank deficiency, all elements of rm is set to zero
 * with a "NULL_VECTOR" return code.
 *						n.kawai 93.01.13
 */

int atRMCheck(
	AtRotMat rm)	/* input/output: rotation matrix to be checked*/
{
    static AtRotMat identity = {
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	};

    AtRotMat rm_inv, rm_prod;
    AtQuat q;
    int i, j;
    double norm, diff;

    norm = 0.;
    ATInvRotMat(rm, rm_inv);
    ATRMProd(rm, rm_inv, rm_prod);
    for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			diff =  rm_prod[i][j] - identity[i][j];
			norm += (diff * diff);
		}
    }

    if (norm < EPS) return (NORMAL_END); /* rm is unitary */

    atRMToQuat(rm, q);
    norm = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] +q[3]*q[3];

    if (fabs(norm) < EPS) return (NULL_VECTOR); /* determinant null */

    norm = sqrt( norm );
    for (i=0; i<3; i++) {		/* normalize quaternion */
		q[i] /=norm;
    };
    atQuatToRM(q, rm);			/* rm is normalized to unitary */
    return (INCONSISTENT_RM);
}
