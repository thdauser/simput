/*============================================================================

  WCSLIB 5.19 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2018, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Direct correspondence concerning WCSLIB to mark@calabretta.id.au

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcsfix_f.c,v 5.19.1.1 2018/07/26 15:41:42 mcalabre Exp mcalabre $
*===========================================================================*/

#include <wcsfix.h>

/* Fortran name mangling. */
#include <wcsconfig_f77.h>
#define wcsfix_  F77_FUNC(wcsfix,  WCSFIX)
#define cdfix_   F77_FUNC(cdfix,   CDFIX)
#define datfix_  F77_FUNC(datfix,  DATFIX)
#define unitfix_ F77_FUNC(unitfix, UNITFIX)
#define celfix_  F77_FUNC(celfix,  CELFIX)
#define spcfix_  F77_FUNC(spcfix,  SPCFIX)
#define cylfix_  F77_FUNC(cylfix,  CYLFIX)

/*--------------------------------------------------------------------------*/

int wcsfix_(int *ctrl, const int naxis[], int *wcs, int stat[])

{
  if (*naxis == 0) naxis = 0x0;
  return wcsfix(*ctrl, naxis, (struct wcsprm *)wcs, stat);
}

/*--------------------------------------------------------------------------*/

int cdfix_(int *wcs)

{
  return cdfix((struct wcsprm *)wcs);
}

/*--------------------------------------------------------------------------*/

int datfix_(int *wcs)

{
  return datfix((struct wcsprm *)wcs);
}

/*--------------------------------------------------------------------------*/

int unitfix_(int *ctrl, int *wcs)

{
  return unitfix(*ctrl, (struct wcsprm *)wcs);
}

/*--------------------------------------------------------------------------*/

int celfix_(int *wcs)

{
  return celfix((struct wcsprm *)wcs);
}

/*--------------------------------------------------------------------------*/

int spcfix_(int *wcs)

{
  return spcfix((struct wcsprm *)wcs);
}

/*--------------------------------------------------------------------------*/

int cylfix_(const int naxis[], int *wcs)

{
  return cylfix(naxis, (struct wcsprm *)wcs);
}
