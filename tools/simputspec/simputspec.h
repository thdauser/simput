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

#ifndef SIMPUTSPEC_H
#define SIMPUTSPEC_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"
#include "arf.h"
#include "rmf.h"

#define TOOLSUB simputspec_main
#include "headas_main.c"


struct Parameters {
  /** File name of the output SIMPUT file. */
  char Simput[SIMPUT_MAXSTR];

  /** EXTNAME of the generated HDU. */
  char Extname[SIMPUT_MAXSTR];
  /** EXTVER of the generated HDU. */
  int Extver;

  /** Lower and upper boundary of the generated spectrum [keV]. */
  float Elow;
  float Eup;
  float Estep;

  /** Power law. */
  float plPhoIndex;
  float plFlux;

  /** Black body temperature [keV]. */
  float bbkT;
  float bbFlux;

  /** Line dispersion [keV]. */
  float flSigma;
  float flFlux;

  float rflSpin;
  float rflFlux;

  /** Absorption column [10^22 atoms/cm^2] */
  float NH;

  /** Reference energy band [keV]. */
  float Emin;
  float Emax;

  /** File name of the input ISIS parameter file containing a spectral
      model. */
  char ISISFile[SIMPUT_MAXSTR];
  /** File name for optional preperation script (f. e. to load additional
   models). */
  char ISISPrep[SIMPUT_MAXSTR];
  char ISISPostCmd[SIMPUT_MAXSTR];

  /** File name of the Xspec spectral model. */
  char XSPECPrep[SIMPUT_MAXSTR];
  char XSPECFile[SIMPUT_MAXSTR];
  char XSPECPostCmd[SIMPUT_MAXSTR];

  /** File name of the input ASCII spectrum. */
  char ASCIIFile[SIMPUT_MAXSTR];

  /** File name of the input PHA spectrum. */
  char PHAFile[SIMPUT_MAXSTR];

  /** number of bins of the spectrum **/
  int nbins;
  /** should the bins be spaced logarithmically? **/
  int logegrid;


};


int simputspec_getpar(struct Parameters* const par);


#endif /* SIMPUTSPEC_H */

