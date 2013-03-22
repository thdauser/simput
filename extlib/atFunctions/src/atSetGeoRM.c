#include "atFunctions.h"
#include "atError.h"

/*
 * calc rotation matrix for converting celestial position (J2000)
 *	to the geographic position on earth
 */
int atSetGeoRM(
	double mjd,	/* input: time in MJD */
	AtRotMat rm)	/* output: rotation matrix  from J2000 to geographic,
			i.e. x-axis at (long, lat)=(0,0), and z at N-pole*/

{
    double gsttod;
    int code1;
    static AtVect z_axis = {0.,0.,1.};

    code1 = atSidereal(mjd, &gsttod);
    atSetRotMat(z_axis, gsttod, rm);
    return code1;
}
