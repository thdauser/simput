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


   Copyright 2015 Thorsten Brand, FAU
*/

#ifndef SIMPUTROTATE_H
#define SIMPUTROTATE_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputrotate_main
#include "headas_main.c"


typedef struct{
  char incat[SIMPUT_MAXSTR];
  char outcat[SIMPUT_MAXSTR];
  float c1_ra;
  float c1_dec;
  float c2_ra;
  float c2_dec;
} Parameters ;

int simputrotate_getpar(Parameters* const par);


#endif /* SIMPUTROTATE_H */

