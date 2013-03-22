			/* bug in sprintf corrected.  n.kawai 92/11/20 */
#include <stdio.h>
#include "atFunctions.h"
#include "atError.h"
#include <string.h>

int                             /* Making formated date & time        */
atCTimeD(                       /*               ver 1.0  92/07/01 ay */
	AtTimeD *attime,	/* input: yr,mo,dy,hr,mn,sc,ss */
	char *ctime)		/* output: formated date (char [25]) */
/*  Note: format of *ctime is "yy/mm/dd hh:mm:ss.ssssss" c*24+\0 */
{
	int  rc, us;
	AtTimeD attimeD;

	attimeD = *attime;
/*	atReformatAtTimeD(&attimeD);*/
	us = (int)(attimeD.ss * 1.0e6);
	if ( us > 1000000 ) {
		us = 999999;
	}
	rc = sprintf(ctime, "%02d/%02d/%02d %02d:%02d:%02d.%06d",
				 attimeD.yr%100, attimeD.mo, attimeD.dy,
				 attimeD.hr, attimeD.mn, attimeD.sc, us);
	if ( 24 == rc ) {
		rc = NORMAL_END;
	}
	return (rc);   /* returns # of char being written if ERROR */
}

int						/* Making formated date & time */
atCTimeD2(				/* ver 1.0  1999/05/13 Y.ISHISAKI */
	AtTimeD *attime,	/* input: yr,mo,dy,hr,mn,sc,ss */
	char *ctime)		/* output: formated date (char [27]) */
/*  Note: format of *ctime is "yyyy/mm/dd hh:mm:ss.ssssss" c*26+\0 */
{
	int  rc, us;
	AtTimeD attimeD;

	attimeD = *attime;
/*	atReformatAtTimeD(&attimeD);*/
	us = (int)(attimeD.ss * 1.0e6);
	if ( us > 1000000 ) {
		us = 999999;
	}
	rc = sprintf(ctime, "%04d/%02d/%02d %02d:%02d:%02d.%06d",
				 attimeD.yr%10000, attimeD.mo, attimeD.dy,
				 attimeD.hr, attimeD.mn, attimeD.sc, us);
	if ( 26 == rc ) {
		rc = NORMAL_END;
	}
	return (rc);   /* returns # of char being written if ERROR */
}

int
atCTime(
	AtTime *attime,		/* input: yr,mo,dy,hr,mn,sc,ms */
	char *ctime)		/* output: formated date (char [25]) */
{
	AtTimeD attimeD;

	atAtTimeToAtTimeD(attime, &attimeD);
	return atCTimeD(&attimeD, ctime);
}

int
atCTime2(
	AtTime *attime,		/* input: yr,mo,dy,hr,mn,sc,ms */
	char *ctime)		/* output: formated date (char [25]) */
{
	AtTimeD attimeD;

	atAtTimeToAtTimeD(attime, &attimeD);
	return atCTimeD2(&attimeD, ctime);
}
