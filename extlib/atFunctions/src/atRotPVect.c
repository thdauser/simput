#include "atFunctions.h"
#include "atError.h"

int               /* Rotating a PolarVector "x" with "ea" to "y" (ROTAT2) */
atRotPVect(       /*                            ver 2.0  93/01/11  n.kawai   */
        AtEulerAng *ea,         /* input */
        AtPolarVect *x,         /* input */
        AtPolarVect *y)         /* output: result */
{
    AtRotMat rm;
    AtVect u, v;
    atEulerToRM( ea, rm);
    atPolToVect(x, u);
    ATRotVect(rm, u, v);
    atVectToPol(v, y);
    return (NORMAL_END);
}
