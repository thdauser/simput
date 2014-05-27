/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef SIMPUTSRC_H
#define SIMPUTSRC_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputsrc_main
#include "headas_main.c"


struct Parameters {
  /** File name of the output SIMPUT file. */
  char Simput[SIMPUT_MAXSTR];

  /** ID of the source. */
  int Src_ID;
  /** Name of source. */
  char Src_Name[SIMPUT_MAXSTR];

  /** Source position [deg]. */
  float RA;
  float Dec;

  /** Reference energy band [keV]. */
  float Emin;
  float Emax;
  /** Reference flux [erg/s/cm^2]. */
  float Flux;
  
  char clobber;
};


int simputsrc_getpar(struct Parameters* const par);


#endif /* SIMPUTSRC_H */

