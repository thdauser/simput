#ifndef SIMPUTFILE_H
#define SIMPUTFILE_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputfile_main
#include "headas_main.c"


struct Parameters {
  /** File name of the output SIMPUT file. */
  char Simput[SIMPUT_MAXSTR];

  /** Name of the X-ray source. */
  char Src_Name[SIMPUT_MAXSTR];

  /** Source position [deg]. */
  float RA;
  float Dec;

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

  /** File name of the input PHA spectrum. */
  char PHAFile[SIMPUT_MAXSTR];

  /** File name of the input ASCII light curve. */
  char LCFile[SIMPUT_MAXSTR];

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

  int chatter;
  char clobber;
  char history;
};


int simputfile_getpar(struct Parameters* const par);


#endif /* SIMPUTFILE_H */

