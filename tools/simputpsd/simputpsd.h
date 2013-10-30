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

