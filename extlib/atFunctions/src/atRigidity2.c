/************************************************************************
  atRigidity2.c

  atRigSet2(AtRigData **rdpp, char *filename)
	Set up Cut-off Rigidity Table in FITS format
  atRigFree2(AtRigData2 *rdp)
	Free allocated memory by atRigSet2()
  atRigidity2(AtRigData2 *rdp, AtPolarVect pvSatG, double *rig)
	Calculate Cut-off Rigidity

	2007/04/16 Y.ISHISAKI	version 2.9
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "atFunctions.h"
#include "atSpline.h"
#include "atError.h"

/* #define TEST */

#ifdef TEST
typedef struct {
	char *filename;
	int nx;
	int ny;
	double altitude;	/* altitude in km from the center of the Earth */
	double *slon;		/* slon[nx] in geocentric coordinates */
	double *slat;		/* slat[ny] in geocentric coordinates */
	double *coeff;		/* coeff[nx][ny] */
	double *data;		/* data[nx][ny] */
	double work[1];		/* work[5*nx] */
} AtRigData2;
#endif

int
atRigSet2(AtRigData2 **rdpp, char *filename)
{
	fitsfile *fp;
	int ix, iy, nx, ny, crpix1, crpix2, anul;
	double crval1, crval2, cdelt1, cdelt2, altitude;
	AtRigData2 *rdp;

	int istat = 0;

	if (
fits_open_file(&fp, filename, READONLY, &istat) ||
		 0 ) {
		return istat;
	}

	if (
fits_read_key(fp, TINT, "NAXIS1", &nx, NULL, &istat) ||
fits_read_key(fp, TINT, "NAXIS2", &ny, NULL, &istat) ||
fits_read_key(fp, TINT, "CRPIX1", &crpix1, NULL, &istat) ||
fits_read_key(fp, TINT, "CRPIX2", &crpix2, NULL, &istat) ||
fits_read_key_dbl(fp, "CRVAL1", &crval1, NULL, &istat) ||
fits_read_key_dbl(fp, "CRVAL2", &crval2, NULL, &istat) ||
fits_read_key_dbl(fp, "CDELT1", &cdelt1, NULL, &istat) ||
fits_read_key_dbl(fp, "CDELT2", &cdelt2, NULL, &istat) ||
fits_read_key_dbl(fp, "ALTITUDE", &altitude, NULL, &istat) ||
		 0 ) {
		return istat;
	}

	if ( 360.0 != (nx-1) * cdelt1 ) {
		return BAD_INDEX_KEY;
	}

	rdp = malloc( sizeof(*rdp) +
				  sizeof(double) * (6*nx + ny + 2*nx*ny - 1) +
				  strlen(filename) + 1 );
	if ( NULL == rdp ) {
		return MEMORY_ALLOCATION;
	}

	rdp->slon = rdp->work + 5*nx;
	rdp->slat = rdp->slon + nx;
	rdp->coeff = rdp->slat + ny;
	rdp->data = rdp->coeff + nx*ny;
	rdp->filename = (char *)(rdp->data + nx*ny);

	strcpy(rdp->filename, filename);
	rdp->nx = nx;
	rdp->ny = ny;
	rdp->altitude = altitude;

	for (ix = 0; ix < nx; ix++) {
		rdp->slon[ix] = crval1 + cdelt1 * (ix + 1 - crpix1);
	}

	for (iy = 0; iy < ny; iy++) {
		rdp->slat[iy] = crval2 + cdelt2 * (iy + 1 - crpix2);
	}

	if (
fits_read_img_dbl(fp, 1, 1, nx*ny, 0.0, rdp->coeff, &anul, &istat) ||
fits_close_file(fp, &istat) ||
		 0 ) {
		free(rdp);
		return istat;
	}

	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			rdp->data[ix*ny + iy] = rdp->coeff[iy*nx+ix];
		}
		if ( rdp->coeff[iy*nx + 0] != rdp->coeff[iy*nx + (nx-1)] ) {
			free(rdp);
			return BAD_DATATYPE;
		}
	}

	istat = atSpline2D(rdp->slon, rdp->slat, rdp->data, nx, ny,
				rdp->coeff, rdp->work);

	if ( istat ) {
		free(rdp);
		return istat;
	}

	*rdpp = rdp;

	return 0;
}

int
atRigFree2(AtRigData2 *rdp)
{
	free(rdp);
	return 0;
}

int
atRigidity2(
	AtRigData2 *rdp,	/* input: rigidity data pointer by atRigSet2() */
	AtPolarVect *pvSatG,/* input: location [km] in geocentric coordinates */
	double *rig)		/* output: cut-off rigidity (GeV/c) */
{
	int nx, ny, istat;
	double rkm, rr, lon, lat, sval;

	nx = rdp->nx;
	ny = rdp->ny;

	rkm = pvSatG->r;
	lon = pvSatG->lon * RAD2DEG;
	lat = pvSatG->lat * RAD2DEG;

	while ( lon < rdp->slon[0] ) {
		lon += 360.0;
	}

	while ( rdp->slon[nx-1] < lon ) {
		lon -= 360.0;
	}

	if ( lat < rdp->slat[0] ) {
		lat = rdp->slat[0];
	} else if ( rdp->slat[ny-1] < lat ) {
		lat = rdp->slat[ny-1];
	}

	istat = atSplint2P(rdp->slon, rdp->slat, rdp->data, rdp->coeff,
		nx, ny, lon, lat, &sval, rdp->work);

	if ( istat ) {
		return istat;
	}

	rr = rdp->altitude / rkm ;
	*rig = sval * rr * rr;

	return 0;
}

#ifdef TEST

static char pname[] = "atRigidity2";

int
main(int argc, char **argv)
{
	AtPolarVect pvSatG;
	AtRigData2 *rdp;
	double altitude, lon, lat, cor, cor2;
	char *filename;
	int istat;

	if ( argc < 2 ) {
		fprintf(stderr, "\
usage: %s FITS-filename\n", pname);
		return 1;
	}

	filename = argv[1];
	istat = atRigSet2(&rdp, filename);
	if ( istat ) {
		fprintf(stderr, "\
%s: atRigSet2() failed (%d)\n", pname, istat);
		return istat;
	}

	printf("rigidity=%s\n", rdp->filename);

	altitude = EARTH_RADIUS + 570.0;
	for (lat = -35.0; lat <= 35.0; lat += 1.0) {
		for (lon = -180.0; lon <= 180.0; lon += 1.0) {
			pvSatG.r = altitude;
			pvSatG.lon = lon * DEG2RAD;
			pvSatG.lat = lat * DEG2RAD;
			atRigidityD(&pvSatG, &cor);
			atRigidity2(rdp, &pvSatG, &cor2);
			printf("\
lat: %+5.1f    lon: %6.1f    cor: %8.3f    cor2: %8.3f    cor%scor2\n",
				lat, lon, cor, cor2, fabs(cor-cor2) < 1e-12 ? "==" : "!=");
		}
	}

	return 0;
}

#endif
