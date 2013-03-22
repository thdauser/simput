#include "atFunctions.h"
#include "atError.h"
#include <math.h>

/* convert *vectPtr1 in old co-odinate system to *vectPtr2 in
   new system defined by the Z-Y-Z Euler angles, eaRad.

	2005/08/06	v2.5	Y.ISHISAKI
		modified to use atan2() instead of atan()
*/
int
atRotEuler2(AtEulerAng *eaRad, AtPolarVect *vectPtr1, AtPolarVect *vectPtr2)
{
	double A, S, C, C2, S2, CDCA, CDSA, SD, lon;

	A = vectPtr1->lon - eaRad->phi;
	S = sin(vectPtr1->lat);
	C = cos(vectPtr1->lat);
	C2 = cos(eaRad->theta);
	S2 = sin(eaRad->theta);

	CDCA = C*cos(A);
	CDSA = C*sin(A);
	SD = S2*CDCA + C2*S;
	CDCA = C2*CDCA - S2*S;

	if ( (1.0 - EPS) <= SD ) {
		vectPtr2->lat = PI*0.5;
		vectPtr2->lon = 0.;
	} else if ( SD <= (-1.0 + EPS) ) {
	   	vectPtr2->lat = -PI*0.5;
	   	vectPtr2->lon = 0.;
	} else {
		vectPtr2->lat = asin(SD);
		S = cos(vectPtr2->lat) - CDCA;
		lon = 2.0*atan2(CDSA, S) - eaRad->psi;
		while ( TWO_PI <= lon ) {
			lon -= TWO_PI;
		}
		while ( lon < 0.0 ) {
			lon += TWO_PI;
		}
		vectPtr2->lon = lon;
	}

	return NORMAL_END;
}
