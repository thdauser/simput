#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/*
 * Finding crossing points of cones		ver 1.0 92/07/01  ay
 *             					ver 1.5  93/01/26  n.kawai
 */

int
atCrossPts(
        AtVect x1,	/* input:  */
        double r1,	/* input: angular distance from x (radian) */
        AtVect x2,	/* input:  */
        double r2,	/* input: angular distance from y (radian) */
        AtVect z[2])	/* output: two crossing points of the two cones */
{
    double  cos1, cos2, v_norm, c_norm, chord;
    AtVect  v, c;

    ATVectProd(x1, x2, v);
    v_norm = ATNorm(v);
    if (v_norm < EPS) return NULL_VECTOR;
	/* Two vectors pointing the same direction */

    cos1 = cos(r1);
    cos2 = cos(r2);

    c[0] = (cos2*(v[1]*x1[2]-v[2]*x1[1])+cos1*(v[2]*x2[1]-v[1]*x2[2]))/v_norm;
    c[1] = (cos2*(v[2]*x1[0]-v[0]*x1[2])+cos1*(v[0]*x2[2]-v[2]*x2[0]))/v_norm;
    c[2] = (cos2*(v[0]*x1[1]-v[1]*x1[0])+cos1*(v[1]*x2[0]-v[0]*x2[1]))/v_norm;

    c_norm = ATNorm(c);
    if (c_norm > 1.+EPS) {
		return NOT_CONVERGED;	/* There is no crossing point */
    } else if (c_norm > 1.- EPS){
		atNormVect(c, z[0]);
		atCopyVect(z[0], z[1]);	/* Degenerate solution- single point */
    } else {
		chord = sqrt( (1.- c_norm)/v_norm );
		atMulAddVect(1., c,  chord, v, z[0]);
		atMulAddVect(1., c, -chord, v, z[1]); /* Double crossing points */
    }

    return NORMAL_END;
}
