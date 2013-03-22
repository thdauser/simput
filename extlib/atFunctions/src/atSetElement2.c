/*
 * set orbital elements (ELMST2)
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
 *	use AtTimeD instead of AtTime
 *
 * 2012/11/22	C.BALUTA & Y.ISHISAKI	version 3.3
 *      - Bug fix in the equation calculating "zzbar",
 *        the first time derivatitive of the mean anomaly
 *      - The following constants have been updated/corrected:
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
 */

#include "fitsio.h"
#include "atFunctions.h"
#include "atError.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

int atSetElement2(	/* return value: condition code */
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
    static AtTimeD attime = { 0,0,0,0,0,0,0. };

/* Local variables */
    double d__1, d__2, d__3, d__4;
    int jblk;
    double cosi, sini, cosi2, sini2;
    int i;
    double e2, zzbar, semax0, semrec;
    double pp1, eee, aj22;
    double ppp, mjd1, mjdz = 0.0;
    double con1, con2;
    static int ibuff[7] = { 0,0,0,0,0,0,0 };
    static char keyword[FLEN_KEYWORD] = "NAXIS2  ";
    double dbuff[13];

/* Local variables added to read FITS orbital file*/
    int   naxis2 = -13;
    static char comment[FLEN_COMMENT] = " blank blank blank";
    int anyf;
    int   nmove = 1,  hdutype;
    double dnull = 0.;
    int inull = 0, status = 0;
    int colnum = 1, felem = 1, nelements = 1, timeelements = 6;
    fitsfile *inf_ptr = NULL;

/* open */
    ffopen(&inf_ptr, filename, READONLY, &status);
    if (status != 0 ) {	/* if specified file is not a proper FITS */
        fprintf( stderr, "\
Specified file is not a FITS orbital file; trying old style\n");
    	return atSetElement(filename, mjd0, kchk);	/* try old style */
    }

    ffmrhd(inf_ptr, nmove, &hdutype, &status);
    ffgky(inf_ptr, TINT, keyword, &naxis2, comment, &status);

    if (status != 0 ) {
/* YI 25-May-1998 */
		ffclos(inf_ptr, &status);
		return FILE_FORMAT_ERROR;
    }

    jblk = 1;
    while (jblk<=naxis2) {
        ffgcvk(inf_ptr, colnum, jblk, felem, timeelements, inull,
                                        &ibuff[1], &anyf, &status );
/*    Should be able to read all 6 numbers into the short int ibuff array */
		if ( status != 0 ) break;
		attime.yr = ibuff[1];
		attime.mo = ibuff[2];
		attime.dy = ibuff[3];
		attime.hr = ibuff[4];
		attime.mn = ibuff[5];
		attime.sc = ibuff[6];
		atMJulianD(&attime, &mjd1);
		if (jblk==1 || mjd1<=mjd0) {
			for (i=2; i<15; i++) {
				ffgcvd(inf_ptr, i, jblk, felem, nelements, dnull,
					   &dbuff[i-2], &anyf, &status);
			}
			mjdz = mjd1;
			jblk++;
		}
		if (mjd1 > mjd0)  break;
	}

    if (jblk == 1 || status != 0 ) {
/* YI 25-May-1998 */
		ffclos(inf_ptr, &status);
		return FILE_FORMAT_ERROR;
    }
    ffclos(inf_ptr, &status);

    atElement.mjdz = mjdz;
    atMJDate(atElement.mjdz, &(atElement.itz));
    atElement.semiax = dbuff[0];
    atElement.eccent = dbuff[1];
    atElement.aincln = dbuff[2] * DEG2RAD;
    atElement.ragome = dbuff[3] * DEG2RAD;
    atElement.smaome = dbuff[4] * DEG2RAD;
    atElement.omean0 = dbuff[5] * DEG2RAD;
    atElement.adot   = dbuff[6];
    atElement.eccdot = dbuff[7];
    atElement.aindot = dbuff[8] * DEG2RAD;
    atElement.ragdot = dbuff[9] * DEG2RAD;
    atElement.smodot = dbuff[10] * DEG2RAD;

/* CONVERT KOZAI'S MEAN TO BROWSLER' MEAN */

    d__1 = sin(atElement.aincln);
    d__2 = sqrt(1 - atElement.eccent * atElement.eccent);
    d__4 = atElement.semiax / EARTH_RADIUS;
    d__3 = 1 - j2 * 1.5 * (1 - d__1 * d__1 * 1.5)
    			/ (d__2 * (d__2 * d__2)) / (d__4 * d__4);
    atElement.semiax /= d__3;
    atElement.adot /= d__3;

/* THEORETICAL VALUE OF TIME DERIVIATIVES OF ELEMENTS */

    d__1 = atElement.semiax, d__2 = d__1;
    pp1 = sqrt(dycon / (d__2 * (d__1 * d__1)));
    ppp = TWO_PI / pp1;
    e2 = atElement.eccent * atElement.eccent;
    semrec = (1. - e2) * atElement.semiax / eradi;
    eee = sqrt(1. - e2);
    cosi = cos(atElement.aincln);
    sini = sin(atElement.aincln);
    cosi2 = cosi * cosi;
    sini2 = sini * sini;
    con1 = aj2 / semrec / semrec;
    aj22 = con1 * con1;
    con1 *= 1.5;
    con2 = aj4 / semrec / semrec / semrec / semrec;

/* zzbar and atElement.znbar always caluculated. revised 93.04.07 */
/* bug fix in the equation calculating "zzbar". revised 2012.11.22 */
    zzbar = pp1 * (con1 * eee * (1. - sini2 * 1.5) + 1. + coef1 * aj22
	    * ((((eee * 25. + 144.) * eee + 105.) * cosi2 - (eee *
	    90. + 96.) * eee + 30.) * cosi2 + (eee * 25. + 16.) *
	    eee - 15.) * eee - coef2 * con2 * e2 * eee * ((cosi2 * 35.
	    - 30.) * cosi2 + 3.));

    atElement.znbar = zzbar * 60.;

    if (atElement.smodot == 0.) {
	atElement.smodot = (con1 * zzbar * (2. - sini2 * 2.5) * (con1 * (e2 *
		.5 + 2. - eee * 2. - (coef3 - e2 / 48. - eee * 3.) *
		sini2) + 1.) - coef4 * aj22 * e2 * pp1 * cosi2 * cosi2 -
		coef5 * con2 * pp1 * ((sini2 * 5.25 - coef6) * sini2 + coef7
		+ e2 * ((coef8 * sini2 - coef9) * sini2 + coef10))) * days;
    }
    if (atElement.ragdot == 0.) {
	atElement.ragdot = -cosi * (zzbar * con1 * (con1 * (e2 / 6. + 1.5 -
		eee * 2. - (coef11 - coef12 * e2 - eee * 3.) * sini2) + 1.)
		 + con2 * 4.375 * pp1 * (e2 * 1.5 + 1.) * (coef13 - coef14
		* sini2)) * days;
    }
/*  PROCESS TO COMPUTE ZNBADT FOR ADOT.NE.0.0 */
    if (atElement.adot != 0. && kchk == 0) {
	semax0 = atElement.semiax / EARTH_RADIUS;
	atElement.znbadt = gravc * -1.5 / sqrt(semax0) / semax0 / semax0 *
		atElement.adot / EARTH_RADIUS;
    }
    atElement.perige = atElement.semiax * (1. - atElement.eccent) - eradi;
    atElement.apoge = atElement.semiax * (atElement.eccent + 1.) - eradi;
    return 0;
}
