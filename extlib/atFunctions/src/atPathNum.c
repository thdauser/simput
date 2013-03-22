#include "atFunctions.h"
#include "atError.h"
#include <math.h>
#include <stdio.h>
#define PATNUMLEN 8

/*
 * calc path number (PATNUM)
 *              Ver 1.1 (92/12/07)  Revised by Ay
 *              Ver 1.6 (93/04/07)  return code Revised by N.Kawai
 */
int atPathNum(
	double mjd,	/* input: time in MJD */
	char path[11])	/* output: path name */
{
    static double mjd0 = 0., mdel0;
    static AtTime time;

    /* Local variables */
    int num, code;
    double e, v, timed, timem, m0, semaxs, mjd1, mut0;

    mjd1 = (double) ((int) (mjd + 5e-9));
    if (mjd1 != mjd0) {
		atMJDate(mjd1, &time);
		mjd0 = mjd1;
		timed = mjd1 - atElement.mjdz;
		timem = timed * 1440.;
		semaxs = atElement.semiax + atElement.adot * timed;
		v = -(atElement.smaome + atElement.smodot * timed);
		e = atan(sqrt((1 - atElement.eccent) / (atElement.eccent + 1))
				 * tan(v / 2)) * 2;
		m0 = (e - atElement.eccent * sin(e));
		mut0 = (atElement.omean0 + (atElement.znbadt * .5 * timed +
									atElement.znbar) * timem);
		while(m0 < 0.) m0 += TWO_PI;
		while(m0 > TWO_PI) m0 += -TWO_PI;
		while(mut0 < 0.) mut0 += TWO_PI;
		while(mut0 > TWO_PI) mut0 += -TWO_PI;
		mdel0 = m0 - mut0;
		if (mdel0 < 0.) mdel0 += TWO_PI;
    }

    timed = mjd - mjd1;
    timem = timed * 1440.;
    num = (int) ((timem * (atElement.znbar + atElement.znbadt * timed)
				  - mdel0) / TWO_PI + 1 + EPS);
    code = (int)sprintf(path, "%02d%02d%02d%02d",
						time.yr%100, time.mo, time.dy, num);
    if ( code == PATNUMLEN || (char *)code == path ) {
		code = NORMAL_END;
    } else {
		code = STRING_FORMAT_ERROR;
    }
    return code;
}
