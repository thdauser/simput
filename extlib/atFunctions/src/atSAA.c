/************************************************************************
  atSAA.c	calc if the point is in the "Brazil Anomaly" (SAA)

  2006/02/08 Y.ISHISAKI	version 2.8
	added new, internally calls atBrazil/atSISBrazil/atSTTBrazil/atHXDBrazil
************************************************************************/

#include "atFunctions.h"
#include "atError.h"

int
atSAA(
	int saa_type,	/* input: 0=nominal, 1=SIS, 2=STT, 3=HXD */
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	int *flag)		/* output: =1 if in SAA, =0 if outside */
{
	double lon, lat;

	lon = x->lon;
	lat = x->lat;

	if ( SAA_NOMINAL == saa_type ) {
		return atBrazil(lon, lat, flag);
	} else if ( SAA_ASCA_SIS == saa_type ) {
		return atSISBrazil(lon, lat, flag);
	} else if ( SAA_ASCA_STT == saa_type ) {
		return atSTTBrazil(lon, lat, flag);
	} else if ( SAA_SUZAKU_HXD == saa_type ) {
		return atHXDBrazil(lon, lat, flag);
	}

	*flag = 0;

	return -1;		/* saa_type not defined */
}
