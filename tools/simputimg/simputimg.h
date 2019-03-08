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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef SIMPUTIMG_H
#define SIMPUTIMG_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputimg_main
#include "headas_main.c"


struct Parameters {
  /** File name of the SIMPUT file the image should be attached to. */
  char Simput[SIMPUT_MAXSTR];

  /** EXTNAME of the generated HDU. */
  char Extname[SIMPUT_MAXSTR];
  /** EXTVER of the generated HDU. */
  int Extver;

  /** File name of the input FITS image. */
  char ImageFile[SIMPUT_MAXSTR];
};


int simputimg_getpar(struct Parameters* const par);


#endif /* SIMPUTIMG_H */
