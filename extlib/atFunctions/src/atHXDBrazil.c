/************************************************************************
  atHXDBrazil.c		calc if the point is in the (SAA) for Suzaku HXD

  2005/09/28 Y.Fukazawa	version 2.5

  2006/02/08 Y.ISHISAKI	version 2.8
	- add comments that longitude & latitude are in geodetic coordinate,
	  which does not consider the ellipticity of the Earth.
************************************************************************/

#include "atFunctions.h"
#include "atError.h"

int
atHXDBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag)	/* output: =1 if in SAA, =0 if outside */
{
	double dlat, dlon;

	dlon = RAD2DEG * lon;
	dlat = RAD2DEG * lat;

	if (dlon >= 180.) dlon -= 360.;

	if ( -85.0 <= dlon && dlon <= 10. && dlat < 2.0 &&
		 dlat <= 0.4*dlon+24. && dlat <= -0.56*dlon-14.4 ) {
	  *flag = 1;
	} else {
	  *flag = 0;
	}

	return 0;
}
