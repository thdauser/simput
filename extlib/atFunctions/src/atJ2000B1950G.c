#include <stdio.h>
#include <math.h>
#include "atFunctions.h"

/* James Peachey, HEASARC/GSFC/NASA, Raytheon STX, 10 June, 1998
   Installing this version of atFunctions into Ftools following
   the release of Ftools 4.1.
   The following #define is commented out so that the C SLA library
   will not be used. The use of the C SLA library would be a copyright
   infringement.
#define SLALIB
*/

#ifdef SLALIB

#include "f2c.h"

static
#include "d_mod.c"
static
#include "d_sign.c"

#include "../slalib/slalib.c"

static double B1950	= 1950.0;

void
atJ2000toB1950(double lo2000, double la2000, double *lo1950, double *la1950)
{
    double dlo1950, dla1950;

    lo2000 *= DEG2RAD;
    la2000 *= DEG2RAD;
    sla_fk54z__(&lo2000, &la2000, &B1950, lo1950, la1950, &dlo1950, &dla1950);
   *lo1950 *= RAD2DEG;
   *la1950 *= RAD2DEG;
}

void
atB1950toJ2000(double lo1950, double la1950, double *lo2000, double *la2000)
{
    lo1950 *= DEG2RAD;
    la1950 *= DEG2RAD;
    sla_fk45z__(&lo1950, &la1950, &B1950, lo2000, la2000);
   *lo2000 *= RAD2DEG;
   *la2000 *= RAD2DEG;
}

void
atJ2000toGal(double ra, double dec, double *gl, double *gb)
{
    ra  *= DEG2RAD;
    dec *= DEG2RAD;
    sla_eqgal__(&ra, &dec, gl, gb);
   *gl  *= RAD2DEG;
   *gb  *= RAD2DEG;
}

void
atGaltoJ2000(double gl, double gb, double *ra, double *dec)
{
    gl  *= DEG2RAD;
    gb  *= DEG2RAD;
    sla_galeq__(&gl, &gb, ra, dec);
   *ra  *= RAD2DEG;
   *dec *= RAD2DEG;
}

void
atB1950toGal(double ra, double dec, double *gl, double *gb)
{
    ra  *= DEG2RAD;
    dec *= DEG2RAD;
    sla_eg50__(&ra, &dec, gl, gb);
   *gl  *= RAD2DEG;
   *gb  *= RAD2DEG;
}

void
atGaltoB1950(double gl, double gb, double *ra, double *dec)
{
    gl  *= DEG2RAD;
    gb  *= DEG2RAD;
    sla_ge50__(&gl, &gb, ra, dec);
   *ra  *= RAD2DEG;
   *dec *= RAD2DEG;
}

#else

static AtRotMat B2J = {
	0.9999256782, -0.0111820611, -0.0048579477,
	0.0111820610,  0.9999374784, -0.0000271765,
	0.0048579479, -0.0000271474,  0.9999881997
};

void
atJ2000toB1950(double lo2000, double la2000, double *lo1950, double *la1950)
{
	static int first_time = 1;
	static AtRotMat J2B;
	AtVect v1950, v2000;
	AtPolarVect pv1950, pv2000;
	if ( first_time ) {
		first_time = 0;
		ATInvRotMat(B2J, J2B);
	}
	pv2000.r = 1.0;
	pv2000.lon = lo2000 * DEG2RAD;
	pv2000.lat = la2000 * DEG2RAD;
	atPolToVect(&pv2000, v2000);
	ATRotVect(J2B, v2000, v1950);
	atVectToPol(v1950, &pv1950);
	*lo1950 = pv1950.lon * RAD2DEG;
	*la1950 = pv1950.lat * RAD2DEG;
}

void
atB1950toJ2000(double lo1950, double la1950, double *lo2000, double *la2000)
{
	AtVect v1950, v2000;
	AtPolarVect pv1950, pv2000;
	pv1950.r = 1.0;
	pv1950.lon = lo1950 * DEG2RAD;
	pv1950.lat = la1950 * DEG2RAD;
	atPolToVect(&pv1950, v1950);
	ATRotVect(B2J, v1950, v2000);
	atVectToPol(v2000, &pv2000);
	*lo2000 = pv2000.lon * RAD2DEG;
	*la2000 = pv2000.lat * RAD2DEG;
}

/* these routine is from att by Hatsukade */
/* convert R.A. Dec. (2000) coodinate to galactic coodinate */
void
atJ2000toGal(double ra, double dec, double *gl, double *gb)
{
	static AtEulerAng eu = {
		192.85 * DEG2RAD,
		62.86666666666667 * DEG2RAD,
		57.06666666666667 * DEG2RAD
	};
	AtPolarVect x,y;
	x.r = 1.0;
	x.lon = ra*DEG2RAD;
	x.lat = dec*DEG2RAD;
	atRotPVect( &eu, &x, &y );
	*gl = y.lon * RAD2DEG;
	*gb = y.lat * RAD2DEG;
}

/* convert galactic coodinate to R.A. Dec. (2000) coodinate */
void
atGaltoJ2000(double gl, double gb, double *ra, double *dec)
{
	static AtEulerAng eu = {
		-57.06666666666667 * DEG2RAD,
		-62.86666666666667 * DEG2RAD,
		-192.85 * DEG2RAD
	};
	AtPolarVect x,y;
	x.r = 1.0;
	x.lon = gl * DEG2RAD;
	x.lat = gb * DEG2RAD;
	atRotPVect(&eu, &x, &y);
	*ra = y.lon * RAD2DEG;
	*dec = y.lat * RAD2DEG;
}

void
atB1950toGal(double ra, double dec, double *gl, double *gb)
{
	atB1950toJ2000(ra, dec, &ra, &dec);
	atJ2000toGal(ra, dec, gl, gb);
}

void
atGaltoB1950(double gl, double gb, double *ra, double *dec)
{
	atGaltoJ2000(gl, gb, ra, dec);
	atJ2000toB1950(*ra, *dec, ra, dec);
}

#endif
