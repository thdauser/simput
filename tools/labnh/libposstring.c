/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include "posstring.h"

#include <math.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>

double parse_posstring(const char *str, const int type) {
  // parse a string giving a position
  // The following syntax is supported:
  //    sDD.DDDD     - e.g.,  -1.234 -- ALWAYS DEGREES
  //    sDD:MM:SS.SS - e.g., +23:15:24.235 see type!
  //    sDD MM SS.SS - e.g.,  23 14 24.234
  //
  // if type is 0  (= DEGPOS), then this string is something in degrees.
  //   no further renormalization is performed
  // if type is 1  (= LATPOS), then this string is a declination and
  //   allowed to be in [-90,+90]
  // if type is 2  (= LONPOS), then this string is a longitude, and it
  //   will be renormalized to [0,360[
  // if type is 3 ( = HOURPOS), then this string is a right ascension
  //   in case the position is given in dd:mm:ss.ss notation, an
  //   automatic conversion to degrees is performed
  //

  // match a double first
  double value;
  int doubstr=1;

  char *endptr;

  value=strtod(str,&endptr);

  // if endptr is not \0, we are not matched fully, i.e.,
  // we aren't a simple double, so get the digits
  if (*endptr!='\0') {
    doubstr=0;

    // determine the sign. Note: any whitespace before
    // sign is already removed by the PIL, so this
    // approach should always work
    double sign=+1;
    if (str[0]=='-' || str[0]=='+') {
      if (str[0]=='-') { sign=-1.;}
    }

    int dd;
    unsigned int mm;
    double ss;

    char buffer1[3],buffer2[3];
    int ret=sscanf(str," %4i%1[: ]%2u%1[: ]%lf",&dd,buffer1,&mm,buffer2,&ss);

    if (ret!=5) {
      if (ret==EOF) {
	err(1,"Error:");
      }
      errx(1,"Must give string of format DD:MM:SS.SSS or DD MM SS.SSS\n");
    }

    value=sign*((ss/60. + mm)/60.+fabs(dd));
  }

  // renormalization and sanity checks
  if (type==LATPOS) {
    if (value<-90. || value>+90.) {
      errx(1,"Degree value must be between -90 and +90 deg (inclusive)\n");
    }
  }

  if (type==LONPOS) {
    // convert a longitude to [0.,360.[
    value=fmod(value,360.);
    if (value<0) {value=value+360.;}
  }

  if (type==HOURPOS) {
    if (doubstr==0) {
      if (value<0. || value>=24.) {
	errx(1,"Hour value must be >=0h and <24h\n");
      }
      value=15.*value;
    } else {
      value=fmod(value,360.);
      if (value<0) {value=value+360.;}
    }

  }

  return value;
}

char *posstring(double position, const int type) {
  // format a position
  // type is as parse_posstring above, the angle position
  // is assumed to be IN DEGREES

  char *result;

  if (type==HOURPOS) { position=position/15.;}

  int sign = (position<0.) ? (int) '-' : (int) '+';
  position=fabs(position);
  int dd=(int) position;
  position=(position-dd)*60.;
  int mm=(int) position;
  position=(position-mm)*60.;

  // round ss to 3 digits for RA, 2 digits otherwise
  double ss;
  if (type==HOURPOS) {
    ss=((long)(position*1000.+0.5))/1000.;
  } else {
    ss=((long)(position*100.+0.5))/100.;
  }

  if (ss>=60.) {
    ss=0.;
    mm=mm+1;
    if (mm>=60) {
      mm=0;
      dd=dd+1;
      if (type==HOURPOS) {
	if (dd==24) dd=0;
      } else {
	if (dd==360) dd=0;
      }
    }
  }

  if ( type==HOURPOS) {
    result=malloc(13);
    if (result==NULL) {
      errx(1,"Out of memory");
    }
    snprintf(result,13,"%02i:%02i:%06.3f",dd,mm,ss);
  } else {
    result=malloc(14);
    if (result==NULL) {
      errx(1,"Out of memory");
    }
    if (type==LATPOS) {
      snprintf(result,13,"%03i:%02i:%05.2f",dd,mm,ss);
    } else {
      snprintf(result,14,"%04i:%02i:%05.2f",dd,mm,ss);
    }
    result[0]=sign;
  }
  return result;
}
