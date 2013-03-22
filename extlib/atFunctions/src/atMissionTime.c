/* aste_time.c,v 1.3 1999/06/03 09:11:13 ishisaki Exp */
/************************************************************************
  atMission.c

	AtTimeD <-> mission time (elapsed seconds since specified MJD)
	conversion routines
	[AtTime とミッションタイム (Astro-E time) を変換するルーチン]

	This routine considers the leap seconds, and the leap seconds
	table should be initialized by atMissionTimeInit().
	Details of the leap seconds can be found in web, e.g.,
		http://tycho.usno.navy.mil/leapsec.html
	as,

	In order to keep the cumulative difference in UT1-UTC less than 0.9
	seconds, a leap second is added to the atomic time to decrease the
	difference between the two. This leap second can be either positive or
	negative depending on the Earth's rotation. Since the first leap
	second in 1972, all leap seconds have been positive and there were 22
	leap seconds in the 27 years to January, 1999.  This pattern reflects
	the general slowing trend of the Earth due to tidal braking.

  followings are history of ascatime (by R.Fujimoto).
  [以下、ascatime (by R.Fujimoto) の履歴]

  93/03/06 V1.0	閏秒を考慮するようにした
	 03/16 V1.1	データをメモリに保存するようにした
	 03/28 V1.2	   閏秒が正しく変換されるようにした
	 06/10 V1.3	   閏秒のテーブルを環境変数で指定
	 06/12 V1.4	   leapflagをNOに初期化した

  95/10/14 V2.0	bug fixならびにmjd2asca、asca2mjdという関数を追加。
	 10/15 V2.1	minor bug fix
	 10/16 V2.3	reformAtTimeの引数に閏秒を加える。
			asca2attimeで閏秒が入った瞬間かどうかの判定を改良。
	 10/16 V2.4	asca2mjdtmpを新設。asca2attimeとasca2mjdを修正。
			閏秒前後でmjdが連続になるように定義。

  1998/03/28 M.Ozaki		V0.0
  		time origin changed for ASTRO-E (2000/1/1 0:0:0 UT)
		[ascatime を元に、時間の 0 点だけを修正 (2000/1/1 0:0:0 UT に)]

  1998/04/13 Y.ISHISAKI		V1.0
		astroetime -> astetime

  1998/07/22 Y.ISHISAKI		V1.1
		bug fix (attime2astroe -> attime2aste) in aste2mjdtmp()

  2004/01/10 Y.ISHISAKI		version 2.2
		MAXDAT increased 30 -> 100
		import from astetool-1.25, and almost rewrite

  2004/01/10 Y.ISHISAKI		version 2.3
		add fflush(NULL), before & after leapsec reading message

  2005/08/14 Y.ISHISAKI		version 2.5
		add atMissionTimeResetTable()
		add +1 leapsec at 2006/01/01 00:00:00 in leap_table[], ndata=24
		remember internal_leap_table[] separately
		accept filename="none" in atMissionTimeInit()
		check (flag < 0) for strict mode, (flag & 1) for verbose mode
		change check_leap_table() from FUNCTION into MACRO
		move to 1st extension only when primary in read_fits()
		check fits_open_file() first, for gzip'ed leapsec file

************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atFunctions.h"
#include "fitsio.h"

static char pname[] = "atMissionTime";

#define ENV_TABLE	"LEAPTABLE"
#define DEF_LEAP	"leapsec.fits"

#define check_leap_table() \
	if ( 0 == leap_table_ready ) {\
		char *filename = getenv(ENV_TABLE);\
		if ( NULL == filename ) {\
			filename = DEF_LEAP;\
		}\
		atMissionTimeInit(filename, 1);\
	}\

typedef struct {
	AtTimeD date;
	long mjd;	/* leap sec must be inserted at 00:00:00 */
	int step;	/* this must be -1/0/1 */
} LEAPSEC_TABLE;

#define MAXDAT 100

static LEAPSEC_TABLE leap_table[MAXDAT], internal_leap_table[] = {
	{ { 1972, 01, 01, 00, 00, 00, 0.0 }, 41317, 0 },
	{ { 1972, 07, 01, 00, 00, 00, 0.0 }, 41499, 1 },
	{ { 1973, 01, 01, 00, 00, 00, 0.0 }, 41683, 1 },
	{ { 1974, 01, 01, 00, 00, 00, 0.0 }, 42048, 1 },
	{ { 1975, 01, 01, 00, 00, 00, 0.0 }, 42413, 1 },
	{ { 1976, 01, 01, 00, 00, 00, 0.0 }, 42778, 1 },
	{ { 1977, 01, 01, 00, 00, 00, 0.0 }, 43144, 1 },
	{ { 1978, 01, 01, 00, 00, 00, 0.0 }, 43509, 1 },
	{ { 1979, 01, 01, 00, 00, 00, 0.0 }, 43874, 1 },
	{ { 1980, 01, 01, 00, 00, 00, 0.0 }, 44239, 1 },
	{ { 1981, 07, 01, 00, 00, 00, 0.0 }, 44786, 1 },
	{ { 1982, 07, 01, 00, 00, 00, 0.0 }, 45151, 1 },
	{ { 1983, 07, 01, 00, 00, 00, 0.0 }, 45516, 1 },
	{ { 1985, 07, 01, 00, 00, 00, 0.0 }, 46247, 1 },
	{ { 1988, 01, 01, 00, 00, 00, 0.0 }, 47161, 1 },
	{ { 1990, 01, 01, 00, 00, 00, 0.0 }, 47892, 1 },
	{ { 1991, 01, 01, 00, 00, 00, 0.0 }, 48257, 1 },
	{ { 1992, 07, 01, 00, 00, 00, 0.0 }, 48804, 1 },
	{ { 1993, 07, 01, 00, 00, 00, 0.0 }, 49169, 1 },
	{ { 1994, 07, 01, 00, 00, 00, 0.0 }, 49534, 1 },
	{ { 1996, 01, 01, 00, 00, 00, 0.0 }, 50083, 1 },
	{ { 1997, 07, 01, 00, 00, 00, 0.0 }, 50630, 1 },
	{ { 1999, 01, 01, 00, 00, 00, 0.0 }, 51179, 1 },
	{ { 2006, 01, 01, 00, 00, 00, 0.0 }, 53736, 1 }
};

static int ndata = sizeof(internal_leap_table) / sizeof(*internal_leap_table);
static char leapfile[1024] = "none";
static int leap_table_ready = 0;

static int
read_fits(fitsfile *fits, char *fn, int tbl_size, LEAPSEC_TABLE *tbl, int flag)
{
	long irow, nrows;
	AtTimeD attime;
	int hdunum, itbl, hdutype, anul;
	struct { int mjd, stp; } col;
	double mjd, stp;
	char head[80], *k;

	int istat = 0;

	sprintf(head, "%s: %s:", pname, (flag < 0) ? "ERROR" : "WARNING");

	fits_get_hdu_num(fits, &hdunum);
	if ( 1 == hdunum ) {	/* primary array */
		fits_movrel_hdu(fits, 1, &hdutype, &istat);
	}
	if ( istat || BINARY_TBL != hdutype ) {
		if ( flag < 0 ) {
			fprintf(stderr, "\
%s: BINARY_TBL extension not found in `%s'\n", head, fn);
		} else {
			fprintf(stderr, "\
%s: BINARY_TBL extension not found in `%s',\n\
	internal table is used\n", head, fn);
		}
		return -1;
	}

	fits_get_num_rows(fits, &nrows, &istat);
	if ( istat || nrows <= 0 ) {
		if ( flag < 0 ) {
			fprintf(stderr, "\
%s: fits_get_num_rows() failed in `%s'\n", head, fn);
		} else {
			fprintf(stderr, "\
%s: fits_get_num_rows() failed in `%s',\n\
	internal table is used\n", head, fn);
		}
		return -1;
	}
	if ( tbl_size < nrows ) {
		if ( flag < 0 ) {
			fprintf(stderr, "\
%s: too many lines (%ld) in `%s'\n", head, nrows, fn);
			return -1;
		} else {
			fprintf(stderr, "\
%s: too many lines (%ld) in `%s',\n\
	truncated to %d.\n", head, nrows, fn, tbl_size);
			nrows = tbl_size;
		}
	}

	if (
fits_get_colnum(fits, CASEINSEN, k="MJD", &col.mjd, &istat) ||
fits_get_colnum(fits, CASEINSEN, k="LEAPSECS", &col.stp, &istat) ||
		 0 ) {
		if ( flag < 0 ) {
			fprintf(stderr, "\
%s: %s column not found in `%s'\n", head, k, fn);
		} else {
			fprintf(stderr, "\
%s: %s column not found in `%s',\n\
	internal table is used\n", head, k, fn);
		}
		return -1;
	}

	itbl = 0;
	for (irow = 1; irow <= nrows; irow++) {
		fits_read_col_dbl(fits, col.mjd, irow, 1, 1, 0.0, &mjd, &anul, &istat);
		fits_read_col_dbl(fits, col.stp, irow, 1, 1, 0.0, &stp, &anul, &istat);
		if ( istat ) {
			if ( flag < 0 ) {
				fprintf(stderr, "\
%s: read error in `%s' at irow=%ld\n", head, fn, irow);
				return -1;
			} else {
				fprintf(stderr, "\
%s: read error in `%s' at irow=%ld,\n\
	remaining lines are ignored\n", head, fn, irow);
				break;
			}
		}

		atMJDateD(mjd, &attime);
		if ( attime.hr || attime.mn || attime.sc || 0.0 != attime.ss ) {
			if ( flag < 0 ) {
				fprintf(stderr, "\
%s: MJD is not an integer in '%s' at irow=%ld\n", head, fn, irow);
				return -1;
			} else {
				fprintf(stderr, "\
%s: MJD is not an integer in '%s' at irow=%ld,\n\
	this line is ignored\n", head, fn, irow);
				continue;
			}
		}

		if ( 0.0 != stp && 1.0 != stp && -1.0 == stp ) {
			if ( flag < 0 ) {
				fprintf(stderr, "\
%s: invalid step of leap seconds in '%s' at irow=%ld,\n\
	LEAPSECS=%f is not an integer within +-1\n", head, fn, irow, stp);
				return -1;
			} else {
				fprintf(stderr, "\
%s: invalid step of leap seconds in '%s' at irow=%ld,\n\
	LEAPSECS=%f is not an integer within +-1,\n\
	this line is ignored", head, fn, irow, stp);
				continue;
			}
		}

		tbl[itbl].date = attime;
		tbl[itbl].mjd = (long)mjd;
		tbl[itbl].step = (int)stp;
		itbl++;

	}

	return itbl;
}

static int
read_text_table(FILE *fp, char *fn, int tbl_size, LEAPSEC_TABLE *tbl)
{
	AtTimeD attime;
	double mjd, stp;

	int nrows = 0;

	for (;;) {
		if ( 8 != fscanf(fp, "%d %d %d %d %d %d %lf %lf",
			 &attime.yr, &attime.mo, &attime.dy,
			 &attime.hr, &attime.mn, &attime.sc, &attime.ss, &stp) ) {
			break;
		}
		if ( nrows < tbl_size ) {
			tbl[nrows].date = attime;
			tbl[nrows].step = stp;
			atMJulianD(&attime, &mjd);
			tbl[nrows].mjd = (long)mjd;
			nrows++;
		} else {
			fprintf(stderr, "\
%s: too many lines in `%s', truncated to %d\n", pname, fn, nrows);
			break;
		}
	}

	return nrows;
}

/************************************************************************
  atMissionTimeResetTable
	reset internal table
************************************************************************/
void
atMissionTimeResetTable(void)
{
	leap_table_ready = 0;
	strcpy(leapfile, "none");
	ndata = sizeof(internal_leap_table) / sizeof(*internal_leap_table);
	memcpy(leap_table, internal_leap_table, sizeof(internal_leap_table));
}

/************************************************************************
  atMissionTimeInit
	read leapsec table (or just return file name)

 input:
	char *filename		leapsec filename to read, NULL for query,
						"none" (case sensitive) for internal table

	int flag			+1: print leapsec table contents
						 0: do not print leapsec table contents
						-1: print table, accept only FITS in argument
						-2: do not print table, accept only FITS in argument
  return value:
	leapsec file name, NULL for ERROR
************************************************************************/
char *
atMissionTimeInit(char *filename, int flag)
{
	fitsfile *fits;
	AtTimeD attime;
	int i;
	double stp;

	FILE *fp = NULL;
	int nrows = -1;
	int istat = 0;

/* if filename == NULL, return current leapsec file name */
	if ( NULL == filename || '\0' == *filename ) {
		if ( leap_table_ready ) {
			return leapfile;	/* already initialized */
		} else {
			return NULL;		/* not initialized */
		}
	}

/* if flag < 0 (strict mode), check whether already initialized */
	if ( flag < 0 && leap_table_ready ) {
		fprintf(stderr, "\
%s: ERROR: leapsec table already initialized\n", pname);
		return NULL;
	}

/* if filename == "none", use internal table */
	if ( 0 == strcmp("none", filename) ) {
		fflush(NULL); printf("\
%s: using internal table ...\n", pname);
		goto skip;
	}

/* try open as FITS file first, because CFITSIO accept gzip'ed file */
	fits_open_file(&fits, filename, READONLY, &istat);
	if ( istat ) {
		if ( flag < 0 ) {	/* strict mode */
			fprintf(stderr, "\
%s: ERROR: leapsec file '%s' open failed\n", pname, filename);
			return NULL;
		}

/* try open as text file next */
		fp = fopen(filename, "r");
		if ( NULL == fp ) {
			fprintf(stderr, "\
%s: WARNING: leapsec file `%s' open failed,\n\
	internal table is used\n", pname, filename);
			goto skip;
		}
	}

/* print information message, when verbose mode */
	if ( flag & 1 ) {
		fflush(NULL); printf("\
%s: reading leapsec file '%s' ...\n", pname, filename);
		fflush(NULL);
	}

/* read file content */
	if ( 0 == istat ) {
		nrows = read_fits(fits, filename, MAXDAT, leap_table, flag);
		fits_close_file(fits, &istat);
		if ( istat ) {
			if ( flag < 0 ) {	/* strict mode */
				fprintf(stderr, "\
%s: ERROR: fits_close_file(`%s') failed\n", pname, filename);
				return NULL;
			} else {
				fprintf(stderr, "\
%s: WARNING: fits_close_file(`%s') failed, ignored\n", pname, filename);
			}
		}
	} else if ( NULL != fp ) {
		nrows = read_text_table(fp, filename, MAXDAT, leap_table);
		fclose(fp);
	}

	if ( 0 < nrows ) {
		ndata = nrows;
	} else if ( flag < 0 ) {	/* strict mode */
		return NULL;
	} else {
 skip:	/* entry for internal table */
		filename = "none";
		atMissionTimeResetTable();
	}

/* print table, when verbose mode */
	if ( flag & 1 ) {
		fflush(NULL);
		for (i = 0; i < ndata; i++) {
			attime = leap_table[i].date;
			stp = leap_table[i].step;
			printf("\
%02d: %04d-%02d-%02d %02d:%02d:%02d.%03.0f   %+.1f sec\n",
				   i+1, attime.yr, attime.mo, attime.dy,
				   attime.hr, attime.mn, attime.sc, attime.ss*1000.0, stp);
		}
		printf("\
%s: %d lines were read\n", pname, ndata);
		fflush(NULL);
	}

/* save file name */
	strncpy(leapfile, filename, sizeof(leapfile) - 1);

/* set table ready flag */
	leap_table_ready = 1;

	return leapfile;
}


/* an example of leapsec

Case 1: step = +1

Note: AtTime not always equal MJD

 AtTime
	|					  ／
	|			 o _ _  ／_ _ _ 1993.6.30 23:59:61
	|		   ／|	  ／
	|		 ／  |  ／
	|	   ／	 |／_ _ _ _ _ _ 1993.7.01 00:00:00 (= mjd_leap)
	|	 ／^
	|  ／  |
	|／____|___________________
		   |				 mission time
		  1993.6.30 23:59:60

   MJD
	|					  ／
	|				    ／
	|				  ／
	|			    ／
	|	   o______／_ _ _ _ _ _ 1993.7.01 00:00:00 (= mjd_leap)
	|	 ／^
	|  ／  |
	|／____|___________________
		   |				 mission time
		  1993.6.30 23:59:60

Case 2: step = -1

 AtTime/MJD
	|		  ／
	|	    ／_ _ _ _ _ _ 1993.7.01 00:00:00 (= mjd_leap)
	|	   |
	|	   |
	|	   |
	|	 ／^
	|  ／  |
	|／____|___________________
		   |				 mission time
		  1993.6.30 23:59:59

*/

/************************************************************************
  atAtTimeDToMission
	convert AtTimeD into mission time

 input:
	double mjd_base		reference MJD (normally integer)
	AtTimeD *attime		yy,mo,dy,hr,mn,sc,ss (0<=hr,mn<60, 0<=sc<60+leap_step)
 output:
	double *mission		elapsed seconds since reference MJD
  return value:
	 0: normal end, not in leap second
	+1: in a positive leap second, e.g., sc=60
	-1: in a negative leap second
************************************************************************/
int
atAtTimeDToMission(double mjd_base, AtTimeD *attime, double *mission_time)
{
	static int DAYSEC = 24 * 60 * 60;
	AtTimeD attimeD;
	double ss, mjd;
	int i, sec_in_a_day, trim;

	check_leap_table();

/* consider sec in a day separately */
	attimeD = *attime;
	sec_in_a_day = (attimeD.hr*60 + attimeD.mn)*60 + attimeD.sc;
	ss = attimeD.ss;
	trim = 0;
	attimeD.hr = attimeD.mn = attimeD.sc = 0;
	attimeD.ss = 0.0;
	atMJulianD(&attimeD, &mjd);
	*mission_time = (mjd - mjd_base) * DAYSEC;

	for (i = 0; i < ndata; i++) {
		if ( leap_table[i].mjd <= mjd_base ) {
			continue;		/* ignore leapsec before base MJD */
		}
		if ( leap_table[i].mjd == (long)mjd + 1 ) {
			if ( leap_table[i].step < 0 ) {
				if ( DAYSEC + leap_table[i].step <= sec_in_a_day ) {
					sec_in_a_day = DAYSEC + leap_table[i].step;
					*mission_time += sec_in_a_day + trim;
					return leap_table[i].step;
				}
			} else if ( 0 < leap_table[i].step ) {
				if ( DAYSEC <= sec_in_a_day ) {
					if ( DAYSEC + leap_table[i].step - 1 < sec_in_a_day ) {
						sec_in_a_day = DAYSEC + leap_table[i].step - 1;
					}
					*mission_time += sec_in_a_day + trim + ss;
					return leap_table[i].step;
				}
			}
		} else if ( leap_table[i].mjd <= mjd ) {
			trim += leap_table[i].step;
		}
	}

	*mission_time += sec_in_a_day + trim + ss;
	return 0;
}

/************************************************************************
  atMissionToAtTimeD
	convert mission time into AtTimeD

 input:
	double mjd_base		reference MJD (normally integer)
	double mission_time	elapsed seconds since reference MJD
 output:
	AtTimeD *attime		yy,mo,dy,hr,mn,sc,ss (0<=hr,mn<60, 0<=sc<60+leap_step)
  return value:
	 0: normal end, not in leap second
	+1: in a positive leap second, e.g., sc=60
	[-1: in a negative leap second, this should not occur]
************************************************************************/
int
atMissionToAtTimeD(double mjd_base, double mission_time, AtTimeD *attime)
{
	static int DAYSEC = 24 * 60 * 60;
	int i;
	double mission_frc;
	double mission_leap, mission_sec, mission_day;	/* here must be double,
										because 2^31 sec = only 68 years */
	check_leap_table();

	if ( mission_time <= 0.0 ) {
	/* ignore leapsec before base MJD */
		atMJDateD(mjd_base, attime);
		attime->ss += mission_time;
		atReformatAtTimeD(attime);
		return 0;
	}

	mission_sec = floor(mission_time);
	mission_frc = mission_time - mission_sec;

	for (i = 1; i < ndata; i++) {
		if ( leap_table[i].mjd < mjd_base ) {
			continue;
		}
		atAtTimeDToMission(mjd_base, &leap_table[i].date, &mission_leap);
		if ( mission_time <= mission_leap ) {
			if ( mission_time == mission_leap ) {
				*attime = leap_table[i].date;
				return 0;
			}
			if ( mission_leap - leap_table[i].step <= mission_sec ) {
				double sc = mission_sec - (mission_leap - leap_table[i].step);
				*attime = leap_table[i].date;
				attime->sc -= leap_table[i].step;
				atReformatAtTimeD(attime);
				attime->sc += leap_table[i].step + (int)sc;
				attime->ss = mission_frc;
				return leap_table[i].step;
			}
			break;
		}
	}

	*attime = leap_table[i-1].date;
	atAtTimeDToMission(mjd_base, attime, &mission_leap);
	mission_sec -= mission_leap;
	mission_day = floor(mission_sec / DAYSEC);
	mission_sec -= mission_day * DAYSEC;
	attime->dy += (int)mission_day;
	attime->sc += (int)mission_sec;
	attime->ss = mission_frc;
	atReformatAtTimeD(attime);

	return 0;
}

/************************************************************************
  atMJDToMission
	convert MJD into mission time

 input:
	double mjd_base		reference MJD (normally integer)
	double mjd			Modified Julian Date
 output:
	double *mission_time	elapsed seconds since reference MJD
  return value:
	 0: normal end, not in leap second
	[+1: in a positive leap second, e.g., sc=60, this should not occur]
	-1: in a negative leap second
************************************************************************/
int
atMJDToMission(double mjd_base, double mjd, double *mission_time)
{
	int step;
	AtTimeD attimeD;

	atMJDateD(mjd, &attimeD);
	step = atAtTimeDToMission(mjd_base, &attimeD, mission_time);

	return step;
}

/************************************************************************
  atMissionToMJD
	convert mission time into MJD

 input:
	double mjd_base		reference MJD (normally integer)
	double mission_time	elapsed seconds since reference MJD
 output:
	double *mjd			Modified Julian Date
  return value:
	 0: normal end, not in leap second
	+1: in a positive leap second, MJD != AtTimeD
	[-1: in a negative leap second, this should not occur]
************************************************************************/
int
atMissionToMJD(double mjd_base, double mission_time, double *mjd)
{
	int step;
	AtTimeD attimeD;

	step = atMissionToAtTimeD(mjd_base, mission_time, &attimeD);
	if ( 0 < step ) {
		attimeD.sc = 60;
		attimeD.ss = 0.0;
	}
	atMJulianD(&attimeD, mjd);

	return step;
}

/*	for Emacs
;;; Local Variables: ***
;;; mode:C ***
;;; tab-width:4 ***
;;; c-indent-level:4  ***
;;; End: ***
*/
