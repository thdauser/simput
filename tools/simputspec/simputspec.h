#ifndef SIMPUTSPEC_H
#define SIMPUTSPEC_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputspec_main
#include "headas_main.c"


struct Parameters {
  /** File name of the output SIMPUT file. */
  char Simput[SIMPUT_MAXSTR];

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

  /** File name of the input ASCII spectrum. */
  char XSPECFile[SIMPUT_MAXSTR];
};


int simputspec_getpar(struct Parameters* const par);


#endif /* SIMPUTSPEC_H */

