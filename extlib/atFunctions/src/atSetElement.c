/*
 * set orbital elements (ELMST2)
 *              Ver 1.6 (93/04/07)  znbar always calculated.  by N.Kawai
 *
 *        Yoshitaka Ishisaki, Tokyo Metro-U, 25 May 1998
 *        float con1, con2 -> double
 *
 * 2004/07/06	Y.ISHISAKI	version 2.3
 *	add declaration of atElement
 *
 * 2004/09/22	Y.ISHISAKI	version 2.3
 *	use AtTimeD instead of AtTime
 *
 * 2012/11/22	C.BALUTA & Y.ISHISAKI	version 3.3
 *     - Bug fix in the equation calculating "zzbar",
 *       the first time derivatitive of the mean anomaly
 *     - The following constants have been updated/corrected:
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

#include "atFunctions.h"
#include "atError.h"
#include <math.h>
#include <stdio.h>

AtElement atElement;

int atSetElement(	/* return value: condition code */
	char *filename,	/* input: path name of the orbital element file */
	double mjd0,	/* input: modified Julian Day */
	int kchk)	/* input: if 0, calculate 2nd derivative of mean
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
    static int ibuff[7] = { 0,0,0,0,0,0,0 };

/* Local variables */
    double d__1, d__2, d__3, d__4;
    int jblk, code1, code2 = 0;
    double cosi, sini, cosi2, sini2;
    int i;
    double e2, zzbar, semax0, semrec;
    double pp1, eee, aj22;
    double ppp, mjd1, mjdz = 0.0;
    double con1, con2;
    FILE *stream;
    double dbuff[11];


/* open */
    stream = fopen(filename, "r");
    if (stream == NULL ) return OPEN_ERROR;

    jblk = 0;

    while (1) {
		for (i=0; i<7; i++) code1 = fscanf(stream, "%d", &ibuff[i]);
		if (code1!=1) break;
		attime.yr = ibuff[1];
		attime.mo = ibuff[2];
		attime.dy = ibuff[3];
		attime.hr = ibuff[4];
		attime.mn = ibuff[5];
		attime.sc = ibuff[6];
		atMJulianD(&attime, &mjd1);
		if (jblk==0 || mjd1<=mjd0) {
			for (i=0; i<11; i++) code2 = fscanf(stream, "%lf", &dbuff[i]);
			mjdz = mjd1;
			jblk++;
		}
		if (mjd1 > mjd0)  break;
    }

    if (jblk == 0 || code2!=1) {
        fclose(stream);
        return FILE_FORMAT_ERROR;
    }

    fclose(stream);

    atElement.mjdz = mjdz;
    atMJDate(atElement.mjdz, &(atElement.itz));
    atElement.semiax = dbuff[0];
    atElement.eccent = dbuff[1];
    atElement.aincln = dbuff[2] * DEG2RAD;
    atElement.ragome = dbuff[3] * DEG2RAD;
    atElement.smaome = dbuff[4] * DEG2RAD;
    atElement.omean0 = dbuff[5] * DEG2RAD;
    atElement.adot   = dbuff[6];
    atElement.ragdot = dbuff[7] * DEG2RAD;
    atElement.smodot = dbuff[8] * DEG2RAD;

    atElement.eccdot = 0.;	/* change of i and e not supported */
    atElement.aindot = 0.;

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
