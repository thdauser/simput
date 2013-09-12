#ifndef SIMPUTMERGE_H
#define SIMPUTMERGE_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"

#define TOOLSUB simputmerge_main
#include "headas_main.c"


struct Parameters {
  char Infile1[SIMPUT_MAXSTR];
  char Infile2[SIMPUT_MAXSTR];
  char Outfile[SIMPUT_MAXSTR];
  char FetchExtensions;
  
  char clobber;
};


int simputmerge_getpar(struct Parameters* const par);


#endif /* SIMPUTMERGE_H */

