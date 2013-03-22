#include "atFunctions.h"
#include <math.h>

/* Coverting Vector to R.A. and Dec in hr/deg,min,s */
/*             ver 1.0  92/07/01  ay    */
/*             ver 1.5  93/01/25  n.kawai    */
int atVectToPol60(
	    AtVect x,		/* input:  vector */
	    AtPolarVect60 *p) 	/* output: polar vector in hh/deg mm ss */
{
	AtPolarVect   pv;
	int  rc;
	double i_min;

	rc = atVectToPol( x, &pv );

	(p->ra).sec =  60. * modf( pv.lon * RAD2DEG * 4.+ EPS, &i_min) ;
	(p->ra).hour = i_min / 60 + EPS;
	(p->ra).min = i_min - 60*(p->ra).hour + EPS;

	if (pv.lat < 0.) {
		(p->dec).sign = -1;
		pv.lat *= -1.;
    } else {
    	(p->dec).sign = 1.;
    }
	(p->dec).sec =  60. * modf( pv.lat * RAD2DEG * 60.+ EPS, &i_min) ;
	(p->dec).deg = i_min / 60 +EPS;
	(p->dec).min = i_min - 60*(p->dec).deg +EPS;

	return ( rc );
}
