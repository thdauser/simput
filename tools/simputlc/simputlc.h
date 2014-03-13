/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef SIMPUTLC_H
#define SIMPUTLC_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputlc_main
#include "headas_main.c"


struct Parameters {
  /** File name of the SIMPUT file the light curve should be attached
      to. */
  char Simput[SIMPUT_MAXSTR];
  
  /** EXTNAME of the generated HDU. */
  char Extname[SIMPUT_MAXSTR];
  /** EXTVER of the generated HDU. */
  int Extver;

  double MJDREF;

  /** File name of the input ASCII light curve. */
  char LCFile[SIMPUT_MAXSTR];
};


int simputlc_getpar(struct Parameters* const par);


#endif /* SIMPUTLC_H */

