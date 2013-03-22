/************************************************************************
  atReformatAtTime()	reformat AtTime

	2004/01/09 Y.ISHISAKI	version 2.2
		code imported from aste_time.c in astetool-2.5, ignoring extrasec
		extended code for the capability of dealing with wide range of values

************************************************************************/

#include "atFunctions.h"
#include "atError.h"
#include <math.h>

int
atReformatAtTimeD(AtTimeD *attime)
{
	static int mdys[13] = {31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	int yr, mo, dy, hr, mn, sc;
	double ss;

	yr = attime->yr;
	mo = attime->mo;
	dy = attime->dy;
	hr = attime->hr;
	mn = attime->mn;
	sc = attime->sc;
	ss = attime->ss;

	if ( ss < 0.0 ) {
		int minus_sc = (int)(-ss + 1);
		if ( ss == (double)(1 - minus_sc) ) {
			minus_sc--;
		}
		sc = sc - minus_sc;
		ss = ss + minus_sc;
	} else if ( 1.0 <= ss ) {
		int plus_sc = (int)ss;
		sc = sc + plus_sc;
		ss = ss - plus_sc;
	}

	if ( sc < 0 ) {
		int minus_mn = (-sc + 59) / 60;
		mn = mn - minus_mn;
		sc = sc + minus_mn * 60;
	} else if ( 60 <= sc ) {
		int plus_mn = sc / 60;
		mn = mn + plus_mn;
		sc = sc - plus_mn * 60;
	}

	if ( mn < 0 ) {
		int minus_hr = (-mn + 59) / 60;
		hr = hr - minus_hr;
		mn = mn + minus_hr * 60;
	} else if ( 60 <= mn ) {
		int plus_hr = mn / 60;
		hr = hr + plus_hr;
		mn = mn - plus_hr * 60;
	}

	if ( hr < 0 ) {
		int minus_dy = (-hr + 23) / 24;
		dy = dy - minus_dy;
		hr = hr + minus_dy * 24;
	} else if ( 24 <= hr ) {
		int plus_dy = hr / 24;
		dy = dy + plus_dy;
		hr = hr - plus_dy * 24;
	}

	for (;;) {
		if ( mo < 1 ) {
			int minus_yr = (-(mo-1) + 11) / 12;
			yr = yr - minus_yr;
			mo = mo + minus_yr * 12;
		} else if ( 13 <= mo ) {
			int plus_yr = (mo-1) / 12;
			yr = yr + plus_yr;
			mo = mo - plus_yr * 12;
		}

		mdys[2] = 28;
		if ( 0 == yr % 4 && (0 == yr % 400 || 0 != yr % 100) ) {
			mdys[2] = 29;	/* leap year */
		}

		if ( dy < 1 ) {
			mo = mo - 1;
			dy = dy + mdys[mo];
		} else if ( mdys[mo] < dy ) {
			dy = dy - mdys[mo % 12];
			mo = mo + 1;
		} else {
			break;
		}
	}

	attime->yr = yr;
	attime->mo = mo;
	attime->dy = dy;
	attime->hr = hr;
	attime->mn = mn;
	attime->sc = sc;
	attime->ss = ss;

	return NORMAL_END;
}

int
atReformatAtTime(AtTime *attime)
{
	int return_code;
	AtTimeD attimeD;

	atAtTimeToAtTimeD(attime, &attimeD);
	return_code = atReformatAtTimeD(&attimeD);
	atAtTimeDToAtTime(&attimeD, attime);

	return return_code;
}
