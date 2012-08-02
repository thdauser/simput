/*============================================================================

  WCSLIB 4.13 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2012, Mark Calabretta

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
  along with WCSLIB.  If not, see <http://www.gnu.org/licenses/>.

  Correspondence concerning WCSLIB may be directed to:
    Internet email: mcalabre@atnf.csiro.au
    Postal address: Dr. Mark Calabretta
                    Australia Telescope National Facility, CSIRO
                    PO Box 76
                    Epping NSW 1710
                    AUSTRALIA

  Author: Mark Calabretta, Australia Telescope National Facility
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id: spx_f.c,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*===========================================================================*/

#include <string.h>

#include <spx.h>

/* Fortran name mangling. */
#include <wcsconfig_f77.h>
#define spxget_ F77_FUNC(spxget, SPXGET)
#define specx_  F77_FUNC(specx,  SPECX)

#define SPX_ERR     200

/*--------------------------------------------------------------------------*/

int spxget_(const int *spx, const int *what, void *value)

{
  int  k;
  int  *ivalp;
  const int *ispxp;
  const struct spxprm *spxp;

  /* Cast pointers. */
  spxp  = (const struct spxprm *)spx;
  ivalp = (int *)value;

  switch (*what) {
  case SPX_ERR:
    /* Copy the contents of the wcserr struct. */
    if (spxp->err) {
      ispxp = (int *)(spxp->err);
      for (k = 0; k < ERRLEN; k++) {
        *(ivalp++) = *(ispxp++);
      }
    } else {
      for (k = 0; k < ERRLEN; k++) {
        *(ivalp++) = 0;
      }
    }
    break;
  default:
    return 1;
  }

  return 0;
}

int spxgti_(const int *spx, const int *what, int *value)
{
  return spxget_(spx, what, value);
}

/*--------------------------------------------------------------------------*/

int specx_(
  const char *type,
  const double *spec,
  const double *restfrq,
  const double *restwav,
  double *specs)

{
  char stype[5];
  strncpy(stype, type, 4);
  stype[4] = '\0';
  return specx(stype, *spec, *restfrq, *restwav, (struct spxprm *)specs);
}

/*--------------------------------------------------------------------------*/

#define SPX_FWRAP(scode, SCODE) \
  int F77_FUNC(scode, SCODE)( \
    const double *rest, \
    const int *n1, \
    const int *s1, \
    const int *s2, \
    const double spec1[], \
    double spec2[], \
    int stat[]) \
  {return scode(*rest, *n1, *s1, *s2, spec1, spec2, stat);}

SPX_FWRAP(freqafrq, FREQAFRQ)
SPX_FWRAP(afrqfreq, AFRQFREQ)

SPX_FWRAP(freqener, FREQENER)
SPX_FWRAP(enerfreq, ENERFREQ)

SPX_FWRAP(freqwavn, FREQWAVN)
SPX_FWRAP(wavnfreq, WAVNFREQ)

SPX_FWRAP(freqvrad, FREQVRAD)
SPX_FWRAP(vradfreq, VRADFREQ)

SPX_FWRAP(freqwave, FREQWAVE)
SPX_FWRAP(wavefreq, WAVEFREQ)

SPX_FWRAP(freqawav, FREQAWAV)
SPX_FWRAP(awavfreq, AWAVFREQ)

SPX_FWRAP(freqvelo, FREQVELO)
SPX_FWRAP(velofreq, VELOFREQ)

SPX_FWRAP(wavevopt, WAVEVOPT)
SPX_FWRAP(voptwave, VOPTWAVE)

SPX_FWRAP(wavezopt, WAVEZOPT)
SPX_FWRAP(zoptwave, ZOPTWAVE)

SPX_FWRAP(waveawav, WAVEAWAV)
SPX_FWRAP(awavwave, AWAVWAVE)

SPX_FWRAP(wavevelo, WAVEVELO)
SPX_FWRAP(velowave, VELOWAVE)

SPX_FWRAP(awavvelo, AWAVVELO)
SPX_FWRAP(veloawav, VELOAWAV)

SPX_FWRAP(velobeta, VELOBETA)
SPX_FWRAP(betavelo, BETAVELO)
