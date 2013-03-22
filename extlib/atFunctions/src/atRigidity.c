/************************************************************************
  atRigidity.c

  atRigidityD(AtPolarVect *x, double *rig)
  atRigidity(AtPolarVect *x, float *rig)	[obsolete]
	calc Cut-off Rigidity for given location on Earth (RIGIDY)

  atRigSet(char *filename) [obsolete, not required]
	Set up Cut-off Rigidity Table (RIGSET), not required since atFunctions-2.7

	2005/12/04 Y.ISHISAKI	version 2.7
		add atRigidityD(), in which atRigSet() is not required now
		close file in atRigSet()
		use periodic spline function by default
		remove global variable atRigidTab

************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atFunctions.h"
#include "atSpline.h"
#include "atError.h"

/* -15 <= lon <= 375 deg, 15 deg step, including 3 duplicated lines */
#define RIG_NX 27
/* -35 <= lat <= +35 deg, 5 deg step */
#define RIG_NY 15

/* define this to use periodic spline interpolation for longitude */
#define USE_PERIODIC_SPLINE

#ifdef USE_PERIODIC_SPLINE
#define MAX_WORK	(RIG_NY < 5*RIG_NX ? 5*RIG_NX : RIG_NY)
#else
#define MAX_WORK	(RIG_NY < 3*RIG_NX ? 3*RIG_NX : RIG_NY)
#endif

static int nx, ny;
static double slong[RIG_NX], slat[RIG_NY];
static double coeff[RIG_NX][RIG_NY];
static double rigidi[RIG_NX][RIG_NY];
static double work[MAX_WORK];

/* altitude [km] for normalization 1/r**2 dependence */
static double rs500 = 6878.142;

/* 3 lines are duplicated for periodicity */
static double rigidity_data[RIG_NX*RIG_NY] = {
/* nlat=15, nlon=27
 --> latitude |south north|
 |
 v  longitude |west|
              |east|
*/
	 5.45,  6.35,  7.65,  8.75,  9.85, 10.65, 11.25, 11.75, 12.05, 12.15, 11.95, 11.55, 10.75,  9.05,  7.45,
	 4.85,  6.05,  6.85,  8.35,  9.55, 10.55, 11.35, 11.95, 12.35, 12.55, 12.45, 11.95, 11.15,  9.65,  7.95,
	 4.25,  5.65,  6.65,  8.15,  9.55, 10.75, 11.55, 12.25, 12.75, 12.95, 12.85, 12.35, 11.65,  9.85,  8.55,
	 3.85,  5.35,  6.75,  8.25,  9.95, 11.05, 11.95, 12.65, 13.15, 13.35, 13.25, 12.75, 11.95, 10.35,  8.75,
	 3.85,  5.25,  6.95,  8.55, 10.25, 11.55, 12.45, 13.15, 13.65, 13.75, 13.65, 13.15, 12.35, 11.15,  9.05,
	 3.65,  5.05,  7.15,  8.85, 10.65, 11.95, 12.95, 13.65, 14.15, 14.25, 14.05, 13.55, 12.75, 11.15,  9.55,
	 3.35,  4.45,  6.85,  9.05, 10.95, 12.35, 13.35, 14.15, 14.55, 14.65, 14.55, 13.95, 13.15, 11.85,  9.55,
	 3.25,  4.55,  6.25,  9.05, 11.25, 12.55, 13.65, 14.35, 14.75, 14.95, 14.75, 14.25, 13.35, 12.15,  9.75,
	 3.15,  4.45,  6.45,  9.15, 11.35, 12.75, 13.75, 14.35, 14.75, 14.85, 14.75, 14.25, 13.45, 12.25, 10.05,
	 2.85,  4.65,  6.35,  9.45, 11.25, 12.65, 13.65, 14.25, 14.55, 14.65, 14.45, 13.95, 13.25, 12.25, 10.25,
	 3.15,  4.65,  6.85,  9.35, 11.25, 12.65, 13.45, 13.95, 14.25, 14.25, 14.05, 13.65, 12.95, 12.05, 10.45,
	 3.65,  5.05,  7.05,  9.05, 11.55, 12.55, 13.35, 13.75, 13.95, 13.85, 13.65, 13.15, 12.45, 11.65,  9.95,
	 4.15,  5.75,  8.15, 10.15, 11.75, 12.65, 13.25, 13.55, 13.65, 13.45, 13.15, 12.65, 11.95, 11.05,  9.25,
	 4.95,  6.75,  9.35, 10.95, 12.05, 12.75, 13.25, 13.35, 13.35, 13.15, 12.75, 12.15, 11.35, 10.25,  8.15,
	 5.75,  8.05,  9.55, 11.35, 12.15, 12.75, 13.05, 13.15, 13.05, 12.75, 12.25, 11.65, 10.85,  9.15,  8.05,
	 6.85,  8.05, 10.35, 11.55, 12.25, 12.65, 12.95, 12.95, 12.75, 12.45, 11.95, 11.15, 10.25,  8.45,  6.65,
	 8.05,  9.25, 10.95, 11.65, 12.15, 12.55, 12.75, 12.65, 12.45, 12.05, 11.45, 10.65,  9.45,  7.75,  5.25,
	 8.45, 10.25, 11.05, 11.65, 12.05, 12.35, 12.45, 12.35, 12.05, 11.65, 10.95,  9.85,  8.35,  5.95,  4.55,
	 9.65, 10.45, 11.05, 11.55, 11.85, 12.05, 12.05, 11.95, 11.55, 10.95,  9.65,  8.35,  6.65,  4.85,  3.65,
	 9.75, 10.35, 10.85, 11.25, 11.55, 11.65, 11.65, 11.45, 11.05,  9.85,  8.15,  6.85,  5.45,  4.05,  2.95,
	 9.55, 10.05, 10.55, 10.95, 11.25, 11.35, 11.35, 11.15, 10.75,  9.65,  7.45,  6.15,  4.95,  3.65,  2.55,
	 8.75,  9.65, 10.15, 10.65, 10.95, 11.25, 11.35, 11.25, 10.95, 10.45,  9.15,  7.35,  5.65,  4.25,  3.25,
	 7.85,  8.65,  9.35, 10.15, 10.75, 11.15, 11.35, 11.55, 11.45, 11.15, 10.65,  9.55,  8.05,  5.85,  4.25,
	 6.55,  7.65,  8.85,  9.55, 10.35, 10.95, 11.35, 11.65, 11.85, 11.75, 11.45, 10.75,  9.85,  8.25,  5.65,
	 5.45,  6.35,  7.65,  8.75,  9.85, 10.65, 11.25, 11.75, 12.05, 12.15, 11.95, 11.55, 10.75,  9.05,  7.45,
	 4.85,  6.05,  6.85,  8.35,  9.55, 10.55, 11.35, 11.95, 12.35, 12.55, 12.45, 11.95, 11.15,  9.65,  7.95,
	 4.25,  5.65,  6.65,  8.15,  9.55, 10.75, 11.55, 12.25, 12.75, 12.95, 12.85, 12.35, 11.65,  9.85,  8.55
};


#if 0
/* conv-rigidity.pl */
#!/usr/local/bin/perl

$nlat = 0;
$nlon = 0;
while (<>) {
	$line = $_;
	(@linedata) = (split);
	if ( $#linedata < 0 ) {
		last;
	}
	if ( 0 == $nlon ) {
		$nlon = $#linedata + 1;
	} elsif ( $nlon != $#linedata + 1 ) {
		print stderr "invalid size of data ignored at:\n";
		print stderr "$line";
		next;
	}
	$data[$nlat] = $line;
	$nlat++;
}

print <<EOF
/* nlat=$nlat, nlon=$nlon
 --> latitude |south north|
 |
 v  longitude |west|
              |east|
*/
EOF
;

for ($i = $0; $i < $nlon; $i++) {
	print "\t";
	for ($j = $nlat - 1; 0 <= $j; $j--) {
		$_ = $data[$j];
		(@linedata) = (split);
		printf "%5s", $linedata[$i];
		if ( $i != $nlon - 1 || $j != 0 ) {
			print ",";
		}
		if ( $j != 0 ) {
			print " ";
		} else {
			print "\n";
		}
	}
}

exit 0;

#endif


static int
atRigInitD(void)
{
	int i, code;

#ifdef USE_PERIODIC_SPLINE
	nx = RIG_NX - 2;
#else
	nx = RIG_NX;
#endif
	ny = RIG_NY;

	for (i = 0; i < nx; i++) {
#ifdef USE_PERIODIC_SPLINE
		slong[i] = i * 15.0;
#else
		slong[i] = (i-1) * 15.0;
#endif
	}

	for (i = 0; i < ny; i++) {
		slat[i] = i*5.0 - 35.0;
	}

#ifdef USE_PERIODIC_SPLINE
	code = atSpline2D(slong, slat, &rigidi[1][0], nx, ny, &coeff[1][0], work);
#else
	code = atSpline2D(slong, slat, &rigidi[0][0], nx, ny, &coeff[0][0], work);
#endif

	return code;
}

/*
 * Set up Cut-off Rigidity Table (RIGSET), not required since atFunctions-2.7
 */
int
atRigSet(
	char *filename)  /* input: path name of the cut-off rigidity file */
{
	FILE *fp;
	int i, j, code;

	/*   INPUT CUTOFF RIGIDITY VALUES FROM filename */

	fp = fopen(filename, "r");
	if ( NULL == fp ) return OPEN_ERROR;

	for (j = RIG_NY-1; 0 <= j; j--) {
		for (i = 0; i < RIG_NX; i++) {
			if ( 1 != fscanf(fp, "%lf", &rigidi[i][j]) ) {
				fclose(fp);
				return OPEN_ERROR;
			}
		}
	}

	fclose(fp);
	code = atRigInitD();
	return code;
}

/*
 * calc Cut-off Rigidity for given location on Earth (RIGIDY)
 */
int
atRigidityD(
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	double *rig)	/* output: cut-off rigidity (GeV/c) */
{
	int code;
	double lon, lat, r;

	if ( 0 == nx || 0 == ny ) {
		nx = RIG_NX-2;
		ny = RIG_NY;
		memcpy(rigidi, rigidity_data, sizeof(rigidi));
		code = atRigInitD();
		if ( code ) return code;
	}

	lon = RAD2DEG * x->lon;
	while ( lon > 360.) lon -= 360.;
	while ( lon < 0.) lon += 360.;

	lat = RAD2DEG * x->lat;
	if ( lat < -35.0 ) {
		lat = -35.0;
	} else if ( 35.0 < lat ) {
		lat = 35.0;
	}

#ifdef USE_PERIODIC_SPLINE
	code = atSplint2P(slong, slat, &rigidi[1][0], &coeff[1][0],
		nx, ny, lon, lat, rig, work);
#else
	code = atSplint2D(slong, slat, &rigidi[0][0], &coeff[0][0],
		nx, ny, lon, lat, rig, work);
#endif

	r = rs500 / x->r ;
	*rig *= (r*r);

	return code;
}

int
atRigidity(
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	float *rig)		/* output: cut-off rigidity (GeV/c) */
{
	int code;
	double cor;
	code = atRigidityD(x, &cor);
	*rig = cor;
	return code;
}
