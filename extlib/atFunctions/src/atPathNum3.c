/*
 * calc path number (PATNUM)
 *              Ver 1.1 (92/12/07)  Revised by Ay
 *              Ver 1.6 (93/04/07)  return code Revised by N.Kawai
 *
 * 2009/07/24	Y.ISHISAKI	version 3.1
 *		modified to use atElement3 in argument
 *		modified to use AtTimeD instead of AtTime
 */

#include "atFunctions.h"
#include "atError.h"
#include <math.h>
#include <stdio.h>

#define PATNUMLEN 8

int
atPathNum3(
	AtElement3 *el,	/* input: contents of orbital elements */
	double mjd,		/* input: time in MJD */
	char path[11]	/* output: path name */
)
{
    static double mjd0 = 0., mdel0;
    static AtTimeD attime;

    /* Local variables */
    int num, code;
    double e, v, timed, timem, m0, semaxs, mjd1, mut0;

    mjd1 = (double) ((int) (mjd + 5e-9));
    if (mjd1 != mjd0) {
		atMJDateD(mjd1, &attime);
		mjd0 = mjd1;
		timed = mjd1 - el->mjdz;
		timem = timed * 1440.;
		semaxs = el->semiax + el->adot * timed;
		v = -(el->smaome + el->smodot * timed);
		e = atan(sqrt((1 - el->eccent) / (el->eccent + 1)) * tan(v / 2)) * 2;
		m0 = (e - el->eccent * sin(e));
		mut0 = (el->omean0 + (el->znbadt * 0.5 * timed + el->znbar) * timem);
		while ( m0 < 0.0 ) m0 += TWO_PI;
		while ( m0 > TWO_PI ) m0 += -TWO_PI;
		while ( mut0 < 0.0 ) mut0 += TWO_PI;
		while ( mut0 > TWO_PI ) mut0 += -TWO_PI;
		mdel0 = m0 - mut0;
		if ( mdel0 < 0.0 ) mdel0 += TWO_PI;
    }

    timed = mjd - mjd1;
    timem = timed * 1440.0;
    num = (int) ((timem * (el->znbar + el->znbadt * timed)
				  - mdel0) / TWO_PI + 1 + EPS);
    code = (int)sprintf(path, "%02d%02d%02d%02d",
						attime.yr%100, attime.mo, attime.dy, num);
    if ( code == PATNUMLEN || (char *)code == path ) {
		code = NORMAL_END;
    } else {
		code = STRING_FORMAT_ERROR;
    }
    return code;
}
