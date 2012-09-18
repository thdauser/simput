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
  $Id: wcsunits_f.c,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*===========================================================================*/

#include <stdio.h>
#include <string.h>

#include <wcsutil.h>
#include <wcsunits.h>

/* Fortran name mangling. */
#include <wcsconfig_f77.h>
#define wcsunitse_ F77_FUNC(wcsunitse, WCSUNITSE)
#define wcsutrne_  F77_FUNC(wcsutrne,  WCSUTRNE)
#define wcsulexe_  F77_FUNC(wcsulexe,  WCSULEXE)

/* Deprecated. */
#define wcsunits_ F77_FUNC(wcsunits, WCSUNITS)
#define wcsutrn_  F77_FUNC(wcsutrn,  WCSUTRN)
#define wcsulex_  F77_FUNC(wcsulex,  WCSULEX)

/*--------------------------------------------------------------------------*/

int wcsunitse_(
  const char have[72],
  const char want[72],
  double *scale,
  double *offset,
  double *power,
  iptr err)

{
  char have_[72], want_[72];

  strncpy(have_, have, 72);
  strncpy(want_, want, 72);
  have_[71] = '\0';
  want_[71] = '\0';

  return wcsunitse(have_, want_, scale, offset, power, (struct wcserr **)err);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int wcsunits_(
  const char have[72],
  const char want[72],
  double *scale,
  double *offset,
  double *power)

{
  return wcsunitse_(have, want, scale, offset, power, 0x0);
}

/*--------------------------------------------------------------------------*/

int wcsutrne_(
  const int *ctrl,
  char unitstr[72],
  iptr err)

{
  int status;
  char unitstr_[72];

  strncpy(unitstr_, unitstr, 72);
  unitstr_[71] = '\0';

  /* This may or may not force the Fortran I/O buffers to be flushed.  If
   * not, try CALL FLUSH(6) before calling WCSUTRNE in the Fortran code. */
  fflush(NULL);

  status = wcsutrne(*ctrl, unitstr_, (struct wcserr **)err);

  wcsutil_blank_fill(72, unitstr_);
  strncpy(unitstr, unitstr_, 72);

  return status;
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int wcsutrn_(
  const int *ctrl,
  char unitstr[72])

{
  return wcsutrne_(ctrl, unitstr, 0x0);
}

/*--------------------------------------------------------------------------*/

int wcsulexe_(
  const char unitstr[72],
  int *func,
  double *scale,
  double units[WCSUNITS_NTYPE],
  iptr err)

{
  char unitstr_[72];

  strncpy(unitstr_, unitstr, 72);
  unitstr_[71] = '\0';

  /* This may or may not force the Fortran I/O buffers to be flushed.  If
   * not, try CALL FLUSH(6) before calling WCSULEXE in the Fortran code. */
  fflush(NULL);

  return wcsulexe(unitstr_, func, scale, units, (struct wcserr **)err);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int wcsulex_(
  const char unitstr[72],
  int *func,
  double *scale,
  double units[WCSUNITS_NTYPE])

{
  return wcsulexe_(unitstr, func, scale, units, 0x0);
}
