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

