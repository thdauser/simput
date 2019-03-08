/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include <stdio.h>
#include <stdlib.h>
#include <pil.h>
#include <err.h>
#include <strings.h>
#include <math.h>

#include "labnh.h"
#include "posstring.h"

// coordinate systems
#define J2000 0
#define GAL 1


int main (int argc,char **argv) {
  int err;
  if ((err=PILInit(argc,argv)) != PIL_OK) {
    errx(1,"Problem in PIL initialization.\n");
  }

  //
  // Verbosity and Sampling
  //
  int verbosity,sampling;
  if ((err=PILGetInt("verbosity",&verbosity))!=PIL_OK) {
    errx(1,"Problem with getting verbosity parameter");
  }
  if ((err=PILGetInt("sampling",&sampling))!=PIL_OK) {
    errx(1,"Problem with getting sampling parameter");
  }
  nh_verbosity(verbosity);
  nh_sampling(sampling);

  //
  // Determine coordinate system used
  //
  char coorsys[PIL_LINESIZE];
  if ((err=PILGetAsString("coordsys",coorsys))!=PIL_OK) {
    errx(1,"Problem with getting coordsys parameter");
  }

  int coord=-1;
  if (strncasecmp(coorsys,"J2000",PIL_LINESIZE)==0) {
    coord=J2000;
  }
  if (strncasecmp(coorsys,"gal",PIL_LINESIZE)==0) {
    coord=GAL;
  }

  if (coord==-1) {
    errx(1,"Illegal coordinate system given");
  }

  double column;

  if (coord==J2000) {
    // ... work in equatorial coordinates
    char rastr[PIL_LINESIZE], decstr[PIL_LINESIZE];
    double ra,dec;

    if ((err=PILGetAsString("ra",rastr))!=PIL_OK) {
      errx(1,"Problem with getting ra parameter");
    }
    ra=parse_posstring(rastr,HOURPOS);

    if ((err=PILGetAsString("dec",decstr))!=PIL_OK) {
      errx(1,"Problem with getting dec parameter");
    }
    dec=parse_posstring(decstr,LATPOS);

    printf("Calculating N_H for position\nRA : %09.5lf deg  DEC: %+08.5f deg\n",
	   ra,dec);
    char *rrst=posstring(ra,HOURPOS);
    char *ddst=posstring(dec,LATPOS);
    if (verbosity>0) {
      printf("     %s        %s\n",rrst,ddst);
    }
    free(rrst);
    free(ddst);

    column=nh_equ(ra,dec);

  }

  if (coord==GAL) {
    // ... work in Galactic coordinates
    char liistr[PIL_LINESIZE], biistr[PIL_LINESIZE];
    double lii,bii;
    if ((err=PILGetAsString("lii",liistr))!=PIL_OK) {
      errx(1,"Problem with getting lii parameter");
    }
    lii=parse_posstring(liistr,LONPOS);

    if ((err=PILGetAsString("bii",biistr))!=PIL_OK) {
      errx(1,"Problem with getting bii parameter");
    }
    bii=parse_posstring(biistr,LATPOS);

    if (verbosity>0) {
      printf("Calculating N_H for position\nlii: %9.5lf deg  bii: %9.5lf deg\n",
	     lii,bii);
    }

    column=nh_gal(lii,bii);

  }

  if((err=PILPutReal("nh",column))!=PIL_OK) {
    errx(1,"Problem writing nh parameter");
  }

  if (verbosity>0) {
    printf("N_H = %e cm^2\n",column);
  }

  PILFlushParameters();
  PILClose(PIL_OK);
  exit(0);
}
