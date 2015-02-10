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

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: tsphdpa.c,v 4.25.1.2 2015/01/06 01:01:52 mcalabre Exp mcalabre $
*=============================================================================
*
* tsphdpa tests sphdpa().
*
*---------------------------------------------------------------------------*/
#include <stdio.h>

#include "sph.h"

int main()

{
  double dist, lat, lat0, lng, lng0, pa;

  printf("\nEnter reference (lng,lat): ");
  scanf("%lf%*[ ,	]%lf", &lng0, &lat0);

  while (1) {
    printf("\nEnter   field   (lng,lat): ");
    scanf("%lf%*[ ,	]%lf", &lng, &lat);

    sphdpa(1, lng0, lat0, &lng, &lat, &dist, &pa);

    printf("(%.4f,%.4f) - (%.4f,%.4f) -> (%.4f,%.4f) (dist,pa)\n",
      lng0, lat0, lng, lat, dist, pa);

    sphpad(1, lng0, lat0, &dist, &pa, &lng, &lat);

    printf("(%.4f,%.4f) + (%.4f,%.4f) -> (%.4f,%.4f) (lng,lat)\n",
      lng0, lat0, dist, pa, lng, lat);
  }

  return 0;
}
