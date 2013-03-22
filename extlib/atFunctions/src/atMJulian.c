/***********************************************************************
  atMJulian()	Converting UT to Modified Julian Day (MJULIA)

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
		treat t->yr as it is, not to add 1900 or 2000
		calculate days in long integer, to avoid floating precision issues

************************************************************************/

#include "atFunctions.h"
#include "atError.h"
#include <math.h>

int		/* Converting UT to Modified Julian Day (MJULIA)     */
atMJulianD(
        AtTimeD *attime,	/* input:  yr,mo,dy,hr,mn,sc,ss */
        double *mjd)		/* output: Modified Julian Day */
{
	static int totdays[13] = {
		0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
	};

	AtTimeD t;
	int y, m, d, leapdays;
	long days;
	double frac;

/* re-format AtTime, first */
	t = *attime;
	atReformatAtTimeD(&t);

/* calculate a fraction in a day */
	frac = ( t.hr + ( t.mn + ( t.sc + t.ss ) /60.0 ) /60.0 ) /24.0;

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
	y = t.yr;
	m = t.mo;
	d = t.dy;
	if ( 0 <= y ) {
		days = y * 365;
		leapdays = (y+3)/4 - (y+99)/100 + (y+399)/400;
		days += leapdays;
	} else {
		y = -y;
		days = y * 365;
		leapdays = y/4 - y/100 + y/400;
		days += leapdays;
		days = -days;
	}

/* add days in a year */
	days += totdays[m-1] + d - 1;
	if ( 0 == y % 4 && (0 == y % 400 || 0 != y % 100) ) {	/* leap year */
		if ( 2 < m ) {
			days++;		/* add +1 for a leap day after March */
		}
	}

/* subtract elapsed days for 1858/11/17 since year=0 */
	days -= 678941;
	*mjd = (double)days + frac;

	return NORMAL_END;
}

int		/* Converting UT to Modified Julian Day (MJULIA)     */
atMJulian(
        AtTime *attime,		/* input:  yr,mo,dy,hr,mn,sc,ms */
        double *mjd)		/* output: Modified Julian Day */
{
	AtTimeD attimeD;

	atAtTimeToAtTimeD(attime, &attimeD);
	return atMJulianD(&attimeD, mjd);
}
