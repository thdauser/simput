/************************************************************************
  atAtTimeDToAtTime()	convert AtTimeD to AtTime

	2004/03/10 Y.ISHISAKI	version 2.2
		created for atFunctions-2.2

************************************************************************/

#include "atFunctions.h"

void
atAtTimeDToAtTime(
	AtTimeD *attimeD,	/* input: yr,mo,dy,hr,mn,sc,ss */
	AtTime *attime)		/* input: yr,mo,dy,hr,mn,sc,ms */
{
	attime->yr = attimeD->yr;
	attime->mo = attimeD->mo;
	attime->dy = attimeD->dy;
	attime->hr = attimeD->hr;
	attime->mn = attimeD->mn;
	attime->sc = attimeD->sc;
	attime->ms = attimeD->ss * 1000.0;
}
