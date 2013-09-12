#ifndef SIMPUTPSD_H
#define SIMPUTPSD_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputpsd_main
#include "headas_main.c"


struct Parameters {
  /** File name of the SIMPUT file, where the PSD should be attached
      to. */
  char Simput[SIMPUT_MAXSTR];
  
  /** File name of the input ASCII PSD. */
  char PSDFile[SIMPUT_MAXSTR];
};


int simputpsd_getpar(struct Parameters* const par);


#endif /* SIMPUTPSD_H */

