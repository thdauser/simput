/*
 * get time of orbital elements, with mjdz <= mjd0 < mjdn, see atSetElement3()
 *
 * 2009/07/24	Y.ISHISAKI	version 3.1
 *		created for atElement3
 */

#include "atFunctions.h"

int
atElementTime3(
	AtElement3 *el,	/* input: contents of orbital elements to be read */
	double *mjdz,	/* output: current Epoch of orbital elements in MJD */
	double *mjdn	/* output: next Epoch of orbital elements in MJD */
)
{
	*mjdz = el->mjdz;
	*mjdn = el->mjdn;

	return 0;
}
