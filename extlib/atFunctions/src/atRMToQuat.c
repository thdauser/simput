#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * convert a Rotation Matrix to Quaternion.		N.Kawai 92.11.20
 *
 *	Check for special cases added			N.Kawai 93.01.11
 *
 * 	Replaced with coord library implementation	B.Wiegand 2003.7.14
 * 		(which is more numerically stable)
 */
int atRMToQuat(
	AtRotMat rm,	/* input: rotation matrix */
	AtQuat q)		/* output: quaternion */
{
	double diag_sum[4];
	double recip;
	int i, maxi;

/****************************************************************
* there are four equivalent ways of calculating this.
* we pick the one with the greatest numerical accuracy
* for this matrix
*****************************************************************/
	diag_sum[0] = 1 + rm[0][0] - rm[1][1] - rm[2][2];
	diag_sum[1] = 1 - rm[0][0] + rm[1][1] - rm[2][2];
	diag_sum[2] = 1 - rm[0][0] - rm[1][1] + rm[2][2];
	diag_sum[3] = 1 + rm[0][0] + rm[1][1] + rm[2][2];

	maxi = 0;
	for (i = 1; i < 4; ++i) {
		if (diag_sum[i] > diag_sum[maxi]) {
			maxi = i;
		}
	}

	q[maxi] = 0.5 * sqrt(diag_sum[maxi]);
	recip = 1. / (4. * q[maxi]);

	if (maxi == 0) {
		q[1] = recip * (rm[0][1] + rm[1][0]);
		q[2] = recip * (rm[2][0] + rm[0][2]);
		q[3] = recip * (rm[1][2] - rm[2][1]);
	}
	else if (maxi == 1) {
		q[0] = recip * (rm[0][1] + rm[1][0]);
		q[2] = recip * (rm[1][2] + rm[2][1]);
		q[3] = recip * (rm[2][0] - rm[0][2]);
	}
	else if (maxi == 2) {
		q[0] = recip * (rm[2][0] + rm[0][2]);
		q[1] = recip * (rm[1][2] + rm[2][1]);
		q[3] = recip * (rm[0][1] - rm[1][0]);
	}
	else if (maxi == 3) {
		q[0] = recip * (rm[1][2] - rm[2][1]);
		q[1] = recip * (rm[2][0] - rm[0][2]);
		q[2] = recip * (rm[0][1] - rm[1][0]);
	}

	return (NORMAL_END);
}
