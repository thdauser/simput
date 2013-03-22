/************************************************************************
  atAtTimeToAtTimeD()	convert AtTime to AtTimeD

	2004/03/10 Y.ISHISAKI	version 2.2
		created for atFunctions-2.2

************************************************************************/

#include "atFunctions.h"

void
atAtTimeToAtTimeD(
	AtTime *attime,		/* input: yr,mo,dy,hr,mn,sc,ms */
	AtTimeD *attimeD)	/* input: yr,mo,dy,hr,mn,sc,ss */
{
	attimeD->yr = attime->yr;
	attimeD->mo = attime->mo;
	attimeD->dy = attime->dy;
	attimeD->hr = attime->hr;
	attimeD->mn = attime->mn;
	attimeD->sc = attime->sc;
	attimeD->ss = attime->ms / 1000.0;
}
