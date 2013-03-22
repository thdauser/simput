/************************************************************************
  atSISBrazil.c		calc if the point is in the (SAA) for ASCA SIS

	original code "dp10SISBrazil.c"   by C.Otani   1994.6.27
	return: =1 if in Deep SAA, =0 if Not Deep SAA or outside

  2006/02/08 Y.ISHISAKI	version 2.8
	- add comments that longitude & latitude are in geodetic coordinate,
	  which does not consider the ellipticity of the Earth.
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"

#define SIS_SAA_QU  -0.0042677
#define SIS_SAA_LI  -0.1666
#define SIS_SAA_CO  -11.33

int
atSISBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag)	/* output: =1 if in SAA, =0 if outside */
{
	double dlat, dlon;

	dlon = RAD2DEG * lon;
	dlat = RAD2DEG * lat;

	if ( dlon >= 180.) dlon -= 360.;

	if ( dlat <= ( SIS_SAA_QU*dlon*dlon + SIS_SAA_LI*dlon + SIS_SAA_CO ) ) {
		*flag = 1;
	} else {
		*flag = 0;
	}

	return 0;
}
