/************************************************************************
  atSTTBrazil.c		calc if the point is in the (SAA) for ASCA STT

	original code "dp10SttBrazil.c"   by C.Otani
	return: =1 if in Deep SAA, =0 if Not Deep SAA or outside

  2006/02/08 Y.ISHISAKI	version 2.8
	- add comments that longitude & latitude are in geodetic coordinate,
	  which does not consider the ellipticity of the Earth.
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atFunctions.h"

int
atSTTBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag)	/* output: =1 if in SAA, =0 if outside */
{
	double dlat, dlon;

	dlon = RAD2DEG * lon;
	dlat = RAD2DEG * lat;

	if (dlon >= 180.) dlon -= 360.;

	if ( -60.0 <= dlon && dlon <= -10.0 && dlat <= -20. &&
		 dlat <= dlon*2./3.+10. && dlat <= dlon*(-2.)-50. ) {
		*flag = 1;
	} else {
		*flag = 0;
	}

	return 0;
}
