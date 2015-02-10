/*============================================================================

  WCSLIB 4.25 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2015, Mark Calabretta

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

  Author: Michael Droetboom, Space Telescope Science Institute,
     and: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: twcs_locale.c,v 4.25.1.2 2015/01/06 01:01:52 mcalabre Exp mcalabre $
*=============================================================================
*
* twcs_locale tests wcslib's handling of locales, such as fr_FR, that use a
* comma as the decimal separator in floating point numbers.  Not part of the
* official test suite.
*
*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <locale.h>

#include "wcs.h"
#include "wcserr.h"
#include "wcshdr.h"
#include "wcsprintf.h"

#define HEADER_SIZE 36000

int main(int argc, char** argv)
{
  struct wcsprm *w;
  char header[HEADER_SIZE];
  size_t real_size;
  FILE *fd;
  int nreject, nwcs;
  int status;
  int nkeyrec;
  char *gen_header;

  wcserr_enable(1);
  wcsprintf_set(stderr);

  setlocale(LC_NUMERIC, "fr_FR");
  wcsprintf("Parsing xmmlss.hdr with locale set to fr_FR.\n");

  fd = fopen("test/xmmlss.hdr", "r");
  real_size = fread(header, 1, HEADER_SIZE, fd);
  fclose(fd);

  if (wcspih(header, real_size / 80, WCSHDR_all, 0, &nreject, &nwcs, &w)) {
    wcserr_prt(w->err, 0x0);
    return 1;
  }

  if (wcsset(w)) {
    wcserr_prt(w->err, 0x0);
    return 1;
  }

  wcsprintf("\nOutput from wcsprt() with this locale\n"
              "-------------------------------------\n");
  wcsprt(w);
  wcsprintf("\n");

  wcsprintf("Output from wcshdo() with the same locale\n"
            "-----------------------------------------\n");
  wcshdo(1, w, &nkeyrec, &gen_header);
  printf("%s", gen_header);

  return 0;
}
