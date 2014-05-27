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

#ifndef SIMPUTPSD_H
#define SIMPUTPSD_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputpsd_main
#include "headas_main.c"


struct Parameters {
  /** File name of the SIMPUT file the PSD should be attached to. */
  char Simput[SIMPUT_MAXSTR];
  
  /** EXTNAME of the generated HDU. */
  char Extname[SIMPUT_MAXSTR];
  /** EXTVER of the generated HDU. */
  int Extver;

  /** PSD general parameters */
  long PSDnpt;
  float PSDfmin;
  float PSDfmax;

  /** PSD: Zero-frequency Lorentzian parameters */
  float LFQ;
  float LFrms;

  /** PSD: Horizontal branch Lorentzian parameters */
  float HBOf;
  float HBOQ;
  float HBOrms;

  /** PSD: Quasi-periodic Lorentzian parameters (1-3) */
  float Q1f;
  float Q1Q;
  float Q1rms;

  float Q2f;
  float Q2Q;
  float Q2rms;

  float Q3f;
  float Q3Q;
  float Q3rms;

  /** File name of the input ASCII PSD. */
  char PSDFile[SIMPUT_MAXSTR];
};


int simputpsd_getpar(struct Parameters* const par);


#endif /* SIMPUTPSD_H */

