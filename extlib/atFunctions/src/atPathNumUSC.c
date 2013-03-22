#include "atFunctions.h"
#include "atError.h"
#include <math.h>
#include <stdio.h>
#define PATNUMLEN 8

/*
 * calc JAXA style path number "YYMMDDnnnn" at USC 34m station
 *		atSetElement() must be called before using this function
 *
 * 2005/10/09	v2.5	Y.ISHISAKI
 *		modified for new path number definition for Suzaku
 */
int
atPathNumUSC(
	double mjd,		/* input: time in MJD */
	char path[11]	/* output: path name, "YYMMDDnnnn" */
)
{
	static int init_flag = 0;
	static AtVect usc34m_vect = {	/* geodetic location of USC 34m station */
		-3586228.639 * 1e-3,
		4114103.918 * 1e-3,
		3290224.771 * 1e-3
	};
	static AtRotMat usc34m_rm;

	int i, num, num_at_dtmin;
	double t, mjdz, dt, dtmin, sat_elev, last_elev;
	AtPolarVect y;
	AtVect sidereal_sat_vec, geodetic_sat_vec;
	AtTimeD attime;

	if ( 0 == init_flag ) {
		atGroundCoord(usc34m_vect, usc34m_rm);
		init_flag = 1;
	}

	mjdz = floor(mjd);
	num_at_dtmin = -1;

 again:

	num = 0;
	dtmin = 999.0;
	last_elev = 999.0;	/* dummy initialization for gcc warning */

	for (i = 0; i <= 24*60; i++) {	/* divide 1 day into 1 min step */
		t = mjdz + i / (24*60.0);
		atSatPos(t, sidereal_sat_vec);
		atGeodetic(t, sidereal_sat_vec, geodetic_sat_vec);
		atAzEl(geodetic_sat_vec, usc34m_vect, usc34m_rm, &y);
		sat_elev = y.lat;	/* elevation in radian, -PI <= sat_elev <= PI */
		if ( 0.0 <= sat_elev ) {
			if ( 0 == i ) {
				last_elev = sat_elev;
			}
			if ( last_elev < 0.0 ) {
				num++;
			}
			dt = fabs(mjd - t);
			if ( dt < dtmin ) {
				dtmin = dt;
				num_at_dtmin = num;
			}
		}
		last_elev = sat_elev;
	}

	if ( 0 == num_at_dtmin ) {	/* belongs to the previous day */
		mjdz = mjdz - 1.0;
		goto again;
	}

	if ( -1 == num_at_dtmin ) {	/* no (0 < elev) pass or always (0 < elev) */
		num_at_dtmin = 0;
	}

	atMJDateD(mjdz, &attime);
	sprintf(path, "%02d%02d%02d%02d00",
		attime.yr%100, attime.mo, attime.dy, num_at_dtmin%100);

	return NORMAL_END;
}
