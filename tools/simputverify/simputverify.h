#ifndef SIMPUTVERIFY_H
#define SIMPUTVERIFY_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputverify_main
#include "headas_main.c"


struct Parameters {
  char Simput[SIMPUT_MAXSTR];
};


int simputverify_getpar(struct Parameters* const par);


#endif /* SIMPUTVERIFY_H */

