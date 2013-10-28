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

