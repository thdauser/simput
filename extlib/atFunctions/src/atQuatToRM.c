#include "atFunctions.h"
#include "atError.h"

/*
 * convert quaternion parameter set to a rotation matrix.  N.Kawai 92.11.20
 */

int atQuatToRM(
	AtQuat q,		/* input: quaternion */
	AtRotMat rm)		/* output: rotation matrix */

{
  rm[0][0] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  rm[0][1] = 2*(q[0]*q[1] + q[2]*q[3]);
  rm[0][2] = 2*(q[0]*q[2] - q[1]*q[3]);
  rm[1][0] = 2*(q[0]*q[1] - q[2]*q[3]);
  rm[1][1] = - q[0]*q[0] + q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  rm[1][2] = 2*(q[1]*q[2] + q[0]*q[3]);
  rm[2][0] = 2*(q[0]*q[2] + q[1]*q[3]);
  rm[2][1] = 2*(q[1]*q[2] - q[0]*q[3]);
  rm[2][2] = - q[0]*q[0] - q[1]*q[1] + q[2]*q[2] + q[3]*q[3];

  return (NORMAL_END);
}
