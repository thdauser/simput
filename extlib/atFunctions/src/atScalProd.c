#include "atFunctions.h"

/*
 * calc Scalar Product (Naiseki, Dot-prooduct) of two vectors
 */
double
atScalProd(		/* output: scaler (inner) product of x and y */
	AtVect x, 	/* input */
	AtVect y)	/* input */
{
	return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}
