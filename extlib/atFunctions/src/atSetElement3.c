/*
 * set orbital elements, with mjdz <= mjd0 < mjdn, see atElementTime3()
 *        originally by Dr. A Yoshida/Dr. Kawai
 *        modified to read the ASTRO-D FITS orbital file by M. Itoh
 *                                                         9208
 *        modified to read new orbital file including time derivatives of
 *        eccentricity and inclination.                    9301
 *
 *        James Peachey, HEASARC/GSFC/NASA, Raytheon STX, 26 Dec 1997
 *        modified to use correct cfitsio calls instead of old fitsio
 *        wrapper macros.
 *
 *        Yoshitaka Ishisaki, Tokyo Metro-U, 25 May 1998
 *        close FITS file on error.
 *        remove unused variables.
 *
 * 2004/09/22	Y.ISHISAKI	version 2.3
 *		use AtTimeD instead of AtTime
 *
 * 2009/07/24	Y.ISHISAKI	version 3.1
 *		modified to use AtElement3 in argument
 *
 * 2012/11/22	C.BALUTA & Y.ISHISAKI	version 3.3
 *		- Bug fix in the equation calculating "zzbar",
 *		  the first time derivatitive of the mean anomaly
 *		- The following constants have been updated/corrected:
 *
 *                             Original       Revised/corrected
 *                             --------------------------------
 *      static double gravc =    .0743666;        .0743668
 *      static double j2 =      .00108264;        .001082626
 *
 *      static double dycon =      398599.;       398600.
 *      static double aj2 =    .001082628;        .001082626
 *      static double aj4 =      -2.12e-6;        -1.616e-6
 *      static double eradi =     6378.16;        6378.137
 *
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "atFunctions.h"
#include "atError.h"
#include "fitsio.h"

static char pname[] = "atSetElement3";

static int
readtext(AtElement3 *el, char *filename, double mjd0)
{
	int i, jblk, ibuff[6];
	double mjdz, mjd1, dbuff[13];
	AtTimeD attimez, attime1;

	FILE *fp;
	int lnum;

	int istat = 0;

	mjdz = 0.0;		/* initialization for gcc warning */
	atMJDateD(mjdz, &attimez);

    fp = fopen(filename, "r");
    if ( NULL == fp ) {
		return OPEN_ERROR;
	}

	for (jblk = 1; ; jblk++) {
/* Should be able to read all 6 numbers into the short int ibuff array */
		if ( 7 != fscanf(fp, "%d%d%d%d%d%d%d", &lnum,
						 &ibuff[0], &ibuff[1], &ibuff[2],
						 &ibuff[3], &ibuff[4], &ibuff[5]) ) {
			if ( 1 == jblk ) {	/* first line */
				istat = FILE_FORMAT_ERROR;
			}
			break;
		}
		attime1.yr = ibuff[0];
		attime1.mo = ibuff[1];
		attime1.dy = ibuff[2];
		attime1.hr = ibuff[3];
		attime1.mn = ibuff[4];
		attime1.sc = ibuff[5];
		attime1.ss = 0.0;
		atMJulianD(&attime1, &mjd1);
		if ( jblk == 1 || mjd1 <= mjd0 ) {
/* Should be able to read all 11 numbers into the double dbuff array */
			for (i = 0; i < 11; i++) {
				if ( 1 != fscanf(fp, "%lf", &dbuff[i]) ) {
					istat = FILE_FORMAT_ERROR;
					break;
				}
			}
			if ( 0 != istat ) break;
			mjdz = mjd1;
			attimez = attime1;
		}
		if ( mjd0 < mjd1 ) break;
	}

	if ( 0 != istat ) {
		goto quit;
	}

	el->mjdz = mjdz;
	el->itz = attimez;
	el->mjdn = mjd1;
	el->itn = attime1;
	el->semiax = dbuff[0];
	el->eccent = dbuff[1];
	el->aincln = dbuff[2] * DEG2RAD;
	el->ragome = dbuff[3] * DEG2RAD;
	el->smaome = dbuff[4] * DEG2RAD;
	el->omean0 = dbuff[5] * DEG2RAD;
	el->adot   = dbuff[6];
	el->eccdot = 0.0;		/* change of i and e not supported */
	el->aindot = 0.0;
	el->ragdot = dbuff[7] * DEG2RAD;
	el->smodot = dbuff[8] * DEG2RAD;

 quit:
	fclose(fp);
	return istat;
}

static int
readfits(AtElement3 *el, char *filename, double mjd0)
{
	int i, jblk, ibuff[6];
	double mjdz, mjd1, dbuff[13];
	AtTimeD attimez, attime1;

	fitsfile *fp;
	int naxis2, anyf, hdutype;
	int nmove = 1;
	int inull = 0;
	int colnum = 1, felem = 1, nelem = 1, timeelements = 6;
	double dnull = 0.0;

	int istat = 0;

	mjdz = 0.0;		/* initialization for gcc warning */
	atMJDateD(mjdz, &attimez);

	fits_open_file(&fp, filename, READONLY, &istat);
	if ( 0 != istat ) {		/* if specified file is not a proper FITS */
/*		fprintf( stderr, "\
%s: Specified file is not a FITS orbital file; trying old style\n", pname);*/
		istat = readtext(el, filename, mjd0);
		return istat;
	}

	if (
fits_movrel_hdu(fp, nmove, &hdutype, &istat) ||
fits_read_key(fp, TINT, "NAXIS2", &naxis2, NULL, &istat) ||
		 0 ) {
		goto quit;
	}

	for (jblk = 1; jblk <= naxis2; jblk++) {
/* Should be able to read all 6 numbers into the short int ibuff array */
		fits_read_col_int(fp, colnum, jblk, felem, timeelements, inull,
						  ibuff, &anyf, &istat );
		if ( 0 != istat ) break;
		attime1.yr = ibuff[0];
		attime1.mo = ibuff[1];
		attime1.dy = ibuff[2];
		attime1.hr = ibuff[3];
		attime1.mn = ibuff[4];
		attime1.sc = ibuff[5];
		attime1.ss = 0.0;
		atMJulianD(&attime1, &mjd1);
		if ( jblk == 1 || mjd1 <= mjd0 ) {
/* Should be able to read all 13 numbers into the double dbuff array */
			for (i = 0; i < 13; i++) {
				fits_read_col_dbl(fp, i+2, jblk, felem, nelem, dnull,
								  &dbuff[i], &anyf, &istat);
				if ( 0 != istat ) break;
			}
			if ( 0 != istat ) break;
			mjdz = mjd1;
			attimez = attime1;
		}
		if ( mjd0 < mjd1 ) break;
	}

	if ( 0 != istat ) {
		goto quit;
	}

	el->mjdz = mjdz;
	el->itz = attimez;
	el->mjdn = mjd1;
	el->itn = attime1;
	el->semiax = dbuff[0];
	el->eccent = dbuff[1];
	el->aincln = dbuff[2] * DEG2RAD;
	el->ragome = dbuff[3] * DEG2RAD;
	el->smaome = dbuff[4] * DEG2RAD;
	el->omean0 = dbuff[5] * DEG2RAD;
	el->adot   = dbuff[6];
	el->eccdot = dbuff[7];
	el->aindot = dbuff[8] * DEG2RAD;
	el->ragdot = dbuff[9] * DEG2RAD;
	el->smodot = dbuff[10] * DEG2RAD;

 quit:
	fits_close_file(fp, &istat);
	return istat;
}

int
atSetElement3(	/* return value: condition code */
	AtElement3 *el,	/* input: contents of orbital elements to be read */
	char *filename,	/* input: the orbital element file name*/
	double mjd0,	/* input: modified Julian Day */
	int kchk)		/* input: if 0, calculate 2nd derivative of mean
								anomaly from A-dot thoretically. */
{
/* Initialized data */
	static double coef1 = .0234375;
	static double coef2 = .3515625;
	static double coef3 = 1.7916667;
	static double coef4 = 1.25;
	static double coef5 = 4.375;
	static double coef6 = 6.642857;
	static double coef7 = 1.7142857;
	static double coef8 = 5.0625;
	static double coef9 = 6.75;
	static double coef10 = 1.928571;
	static double coef11 = 1.666667;
	static double coef12 = .2083333;
	static double coef13 = .8571428;
	static double coef14 = 1.5;
	static double gravc = .0743668;
	static double j2 = .001082626;
	static double days = 86400.;
	static double dycon = 398600.;
	static double aj2 = .001082626;
	static double aj4 = -1.616e-6;
	static double eradi = 6378.137;

/* Local variables */
	double d1, d2, d3, d4;
	double cosi, sini, cosi2, sini2;
	double e2, zzbar, semax0, semax0_2_5, semrec, semrec2, semrec4;
	double pp1, eee, aj22, ppp;
	double con1, con2;

	int istat;

	istat = readfits(el, filename, mjd0);
	if ( 0 != istat ) {
		return istat;
	}

/* CONVERT KOZAI'S MEAN TO BROWSLER' MEAN */

	d1 = sin(el->aincln);
	d2 = sqrt(1.0 - el->eccent * el->eccent);
	d4 = el->semiax / EARTH_RADIUS;
	d3 = 1.0 - j2 * 1.5 * (1.0 - d1 * d1 * 1.5) / (d2 * (d2 * d2)) / (d4 * d4);
	el->semiax /= d3;
	el->adot /= d3;

/* THEORETICAL VALUE OF TIME DERIVIATIVES OF ELEMENTS */

	d1 = el->semiax, d2 = d1;
	pp1 = sqrt(dycon / (d2 * (d1 * d1)));
	ppp = TWO_PI / pp1;
	e2 = el->eccent * el->eccent;
	semrec = (1.0 - e2) * el->semiax / eradi;
	semrec2 = semrec * semrec;
	semrec4 = semrec2 * semrec2;
	eee = sqrt(1.0 - e2);
	cosi = cos(el->aincln);
	sini = sin(el->aincln);
	cosi2 = cosi * cosi;
	sini2 = sini * sini;
	con1 = aj2 / semrec2;
	aj22 = con1 * con1;
	con1 *= 1.5;
	con2 = aj4 / semrec4;

/* zzbar and el->znbar always caluculated. revised 93.04.07 */
/* bug fix in the equation calculating "zzbar". revised 2012.11.22 */
	zzbar = pp1 * (con1 * eee * (1. - sini2 * 1.5) + 1. + coef1 * aj22
		* ((((eee * 25. + 144.) * eee + 105.) * cosi2 - (eee *
		90. + 96.) * eee + 30.) * cosi2 + (eee * 25. + 16.) *
		eee - 15.) * eee - coef2 * con2 * e2 * eee * ((cosi2 * 35.
		- 30.) * cosi2 + 3.));

	el->znbar = zzbar * 60.0;

	if ( el->smodot == 0.0 ) {
		el->smodot = (con1 * zzbar * (2. - sini2 * 2.5) * (con1 * (e2 *
			.5 + 2. - eee * 2. - (coef3 - e2 / 48. - eee * 3.) *
			sini2) + 1.) - coef4 * aj22 * e2 * pp1 * cosi2 * cosi2 -
			coef5 * con2 * pp1 * ((sini2 * 5.25 - coef6) * sini2 + coef7
			+ e2 * ((coef8 * sini2 - coef9) * sini2 + coef10))) * days;
	}

	if ( el->ragdot == 0.0 ) {
			el->ragdot = -cosi * (zzbar * con1 * (con1 * (e2 / 6. + 1.5 -
			eee * 2. - (coef11 - coef12 * e2 - eee * 3.) * sini2) + 1.)
			+ con2 * 4.375 * pp1 * (e2 * 1.5 + 1.) * (coef13 - coef14
			* sini2)) * days;
	}

/*  PROCESS TO COMPUTE ZNBADT FOR ADOT.NE.0.0 */
	if ( el->adot != 0.0 && kchk == 0 ) {
		semax0 = el->semiax / EARTH_RADIUS;
		semax0_2_5 = sqrt(semax0) * semax0 * semax0;
		el->znbadt = -1.5 * gravc / semax0_2_5 * el->adot / EARTH_RADIUS;
	}
	el->perige = el->semiax * (1.0 - el->eccent) - eradi;
	el->apoge = el->semiax * (el->eccent + 1.0) - eradi;

	return 0;
}
