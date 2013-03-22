/************************************************************************
  atBrazil.c	calc if the point is in the "Brazil Anomaly" (SAA)

  2006/02/08 Y.ISHISAKI	version 2.8
	- slight modification of SAA area to be continuous at connecting longitude,
	  but still concave at (lon, lat) = (30, -77), possible short-time SAA-out
	  may appear for high lattitude (lat < -77) satellites.
	- add comments that longitude & latitude are in geodetic coordinate,
	  which does not consider the ellipticity of the Earth.

  2006/06/02 Y.ISHISAKI	version 2.8.1
	- bug fix of atBrazil() in v2.8, "-" sign was missing,
	     (-34.0 <= dlon && dlon < 30.0 && 64*dlat <  17*(dlon + 34) )
	  -> (-34.0 <= dlon && dlon < 30.0 && 64*dlat < -17*(dlon + 34) )
	  pointed out by the Swift team. Thanks.
************************************************************************/

#include "atFunctions.h"
#include "atError.h"

int
atBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag)	/* output: =1 if in SAA, =0 if outside */
{
	double r1, r2, r3, r4;
	double dlat, dlon;

	dlon = RAD2DEG * lon;
	dlat = RAD2DEG * lat;

	if (dlon >= 180.) dlon -= 360.;

	r1 = dlon + 34.;
	r2 = dlat + 70.;
	r3 = dlon - 30.;
	r4 = dlat + 47.;

    if (  -83.0 < dlon && dlon < 60.0 &&
	    ( ( dlon < -34.0 && r1*r1 + r2*r2 < 4900.0 ) ||
	      (-34.0 <= dlon && dlon < 30.0 && 64*dlat < -17*(dlon + 34) ) ||
	      ( 30.0 <= dlon && r3*r3 + r4*r4 < 900.0 )
	    ) ) {
		*flag = 1;
	} else {
		*flag = 0;
	}

	return 0;
}
