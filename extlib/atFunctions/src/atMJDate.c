/***********************************************************************
  atMJDate()	Converting Modified Julian Day to UT (MJDATE)

	Originally coded in Fortran (MJULIA) by K.Mitsuda.
	Translated in C and revised by A.Yoshida for ASTRO-D (ASCA).
	Note: Day of the week is given by ( (mjd+3)%7 )  (0=SUN, 1=MON, etc)
	Be carefull, MJD has only ms accuracy with double precision floating point

	Modified Julian Date = (JD - 2400000.5),
	MJD = 0.0 corresponds 1858/11/17 00:00:00 (Wed)
	JD = 0.0 corresponds 12:00:00 noon 1 Jan -4712 (4713 BC) in Julian Calendar

	There is a gap of ten days between Julian Caldndar and
	Gregorian Calendar, but no discontinuity in Julian dates
	or days of the week: 4 October 1582 (Julian) is a Thursday,
	which begins at JD 2299159.5; and 15 October 1582 (Gregorian) is a Friday,
	which begins at JD 2299160.5. This routine assumes that the
	AtTime structure is always stored in Gregorian Calender
	even before 1582/10/15.

	In the Gregorian calendar currently in use worldwide,
	there is a leap year every year divisible by four except for years
	which are both divisible by 100 and not divisible by 400.
	Therefore, the year 2000 is a leap year, but the years 1700, 1800,
	and 1900 are not.

	Reference: Lieske,J.H., 1979. Astron.Astrophys.,73,282.
	http://aa.usno.navy.mil/data/docs/JulianDate.html

	1992/07/02 A.YOSHIDA	version 1.0

	2004/01/09 Y.ISHISAKI	version 2.2 (almost rewrited)
		changed to use t->ss instead of t->ms
		calculate days in long integer, to avoid floating precision issues

************************************************************************/

#include <stdio.h>
#include <math.h>
#include "atFunctions.h"
#include "atError.h"

static char *pname = "atMJDate";

static char month_norm[365] = {
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,	/* 31 */
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,			/* 28 */
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,	/* 31 */
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,		/* 30 */
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,	/* 31 */
 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,		/* 30 */
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,	/* 31 */
 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,	/* 31 */
 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,		/* 30 */
10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,	/* 31 */
11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,		/* 30 */
12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12	/* 31 */
};

static char month_leap[366] = {
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,	/* 31 */
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,			/* 29 */
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,	/* 31 */
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,		/* 30 */
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,	/* 31 */
 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,		/* 30 */
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,	/* 31 */
 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,	/* 31 */
 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,		/* 30 */
10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,	/* 31 */
11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,		/* 30 */
12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12	/* 31 */
};

static int totdays_norm[13] = {
	0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
};

static int totdays_leap[13] = {
	0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366
};

int		/* Converting Modified Julian Day to UT (MJDATE)   */
atMJDateD(
        double mjd,			/* input:  Modified Julian Day */
        AtTimeD *attime)	/* output: yr,mo,dy,hr,mn,sc,ss */
{
	long days, subdays;
	double frac, ss;
	int y, m, d, hr, mn, sc, leapdays;

	days = (long)floor(mjd);
	frac = mjd - (double)days;

/* add elapsed days for 1858/11/17 since year=0 */
	days += 678941;

/* calculate the elapsed days since year=0
   y	365*y	leapdays
				/4	/100 /400
   -4	-1460	-1   +0   -0
   -3	-1095	-0   +0   -0
   -2	-730	-0   +0   -0
   -1	-365	-0   +0   -0
   0	0		+0   -0   -0
   1	365		+1   -1   +1
   2	730		+1   -1   +1
   3	1095	+1   -1   +1
   4	1460	+1   -1   +1
   5	1825	+2   -1   +1
   100	36500	+25  -1   +1
   101	36865	+26  -2   +1
   400	146000	+100 -4   +1
   401	146365	+101 -5   +2
*/
	if ( 0 <= days ) {
		y = days / 365;
		leapdays = (y+3)/4 - (y+99)/100 + (y+399)/400;
		for (;;) {
			y = (days - leapdays) / 365;
			leapdays = (y+3)/4 - (y+99)/100 + (y+399)/400;
			subdays = days - y*365 - leapdays;
			if ( subdays < 0 ) {	/* this can't be happened */
				fprintf(stderr, "\
%s: calculation error for MJD=%.8f (subdays=%ld)\n", pname, mjd, subdays);
				return NOT_CONVERGED;
			}
			if ( 0 == y % 4 && (0 == y % 400 || 0 != y % 100) ) {
				/* leap year */
				if ( subdays < 366 ) {
					m = month_leap[subdays];
					d = subdays - totdays_leap[m-1] + 1;
					break;
				}
			} else {
				if ( subdays < 365 ) {
					m = month_norm[subdays];
					d = subdays - totdays_norm[m-1] + 1;
					break;
				}
			}
		}
	} else {
		y = (-days + 364) / 365;
		leapdays = y/4 - y/100 + y/400;
		for (;;) {
			y = (-days + 364 - leapdays) / 365;
			leapdays = y/4 - y/100 + y/400;
			subdays = y * 365 + leapdays + days;
			if ( subdays < 0 ) {	/* this can't be happened */
				fprintf(stderr, "\
%s: calculation error for MJD=%.8f (subdays=%ld)\n", pname, mjd, subdays);
				return NOT_CONVERGED;
			}
			if ( 0 == y % 4 && (0 == y % 400 || 0 != y % 100) ) {
				/* leap year */
				if ( subdays < 366 ) {
					m = month_leap[subdays];
					d = subdays - totdays_leap[m-1] + 1;
					break;
				}
			} else {
				if ( subdays < 365 ) {
					m = month_norm[subdays];
					d = subdays - totdays_norm[m-1] + 1;
					break;
				}
			}
		}
		y = -y;
	}

	hr = mn = sc = 0;
	ss = 0.0;

	if ( 0.0 != frac ) {
		frac *= 24;
		hr = (int)frac;
		frac -= hr;
		if ( 0.0 != frac ) {
			frac *= 60;
			mn = (int)frac;
			frac -= mn;
			if ( 0.0 != frac ) {
				frac *= 60;
				sc = (int)frac;
				frac -= sc;
				ss = frac;
			}
		}
	}

	attime->yr = y;
	attime->mo = m;
	attime->dy = d;
	attime->hr = hr;
	attime->mn = mn;
	attime->sc = sc;
	attime->ss = ss;

	return NORMAL_END;
}

int		/* Converting Modified Julian Day to UT (MJDATE)   */
atMJDate(
        double mjd,			/* input:  Modified Julian Day */
        AtTime *attime)		/* output: yr,mo,dy,hr,mn,sc,ss */
{
	int return_code;
	AtTimeD attimeD;

	return_code = atMJDateD(mjd, &attimeD);
	atAtTimeDToAtTime(&attimeD, attime);

	return return_code;
}
