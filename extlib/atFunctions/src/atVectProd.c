/************************************************************************
  atVectProd.c		Making a vector prodact

	1992/07/01 A.YOSHIDA	version 1.0

	2005/12/04 Y.ISHISAKI	version 2.7
		use fast macro
************************************************************************/

#include "atFunctions.h"
#include "atError.h"

/*
 * Making a vector prodact
 */
int
atVectProd(
        AtVect x,       /* input */
        AtVect y,       /* input */
        AtVect z)       /* output: vector (outer) product of x and y */
{
	return ATVectProd(x, y, z);
}
