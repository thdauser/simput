#ifndef SIMPUTLC_H
#define SIMPUTLC_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputlc_main
#include "headas_main.c"


struct Parameters {
  /** File name of the SIMPUT file the light curve should be attached
      to. */
  char Simput[SIMPUT_MAXSTR];
  
  /** File name of the input ASCII light curve. */
  char LCFile[SIMPUT_MAXSTR];
};


int simputlc_getpar(struct Parameters* const par);


#endif /* SIMPUTLC_H */

