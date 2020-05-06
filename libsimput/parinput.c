/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Thomas Dauser, Joern Wilms FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "parinput.h"
#include "common.h"
#include "simput.h"
#include "ape/ape_trad.h"

#define PAR_FLOAT 1
#define PAR_DOUBLE 2
#define PAR_INT 3
#define PAR_LONG 4
#define PAR_BOOL 5


//
// dispatcher routine for simple parameters
//
void query_simput_parameter(char* name, const int type, void *retval, int* status ){

  if (*status!=EXIT_SUCCESS) return;

  switch (type) {
  case PAR_FLOAT:
      *status=ape_trad_query_float(name,(float *) retval);
      break;
  case PAR_DOUBLE:
      *status=ape_trad_query_double(name,(double *) retval);
      break;
  case PAR_INT:
      *status = ape_trad_query_int(name,(int *) retval);
      break;
  case PAR_LONG:
      *status = ape_trad_query_long(name,(long *) retval);
      break;
  case PAR_BOOL:
      *status = ape_trad_query_bool(name,(char *) retval);
      break;
  }

  if (*status!=EXIT_SUCCESS) {
    fprintf(stderr,"Failed to read parameter %s from the command line. \n",name);
  }
}

// return 1 (boolean true) if *file is "none" (in any of its capitalized versions),
// return 0 (boolean false) otherwise.
int is_empty_file_name(char *file){
  if (strncasecmp(file, "none", 5) == 0) {
    return 1;
  } else {
    return 0;
  }
}

// Read a file name parameter and return it in a pointer that can be freed.
// If name is "none" (in any of its capitalized versions), *field is set to NULL.
void query_simput_parameter_file_name(char *name, char **field, int *status){
  if (*status!=EXIT_SUCCESS) return;

  char *buf=NULL;
  *status=ape_trad_query_file_name(name, &buf);
  if (*status!=EXIT_SUCCESS) return;

  if (is_empty_file_name(buf)) {
    *field = NULL;
  } else {
    *field=strdup(buf);
    free(buf);
  }
}

// read a file name parameter into a preallocated buffer
void query_simput_parameter_file_name_buffer(char *name, char * const field, int buflen, int *status){
  if (*status!=EXIT_SUCCESS) return;
  char *buf=NULL;
  *status=ape_trad_query_file_name(name, &buf);
  if (*status!=EXIT_SUCCESS) return;

  if (is_empty_file_name(buf)) {
    field[0]='\0';
  } else {
    strncpy(field,buf,buflen-1);
    field[buflen-1]='\0';
  }
}

// read a string parameter and return it in a pointer that can be freed
void query_simput_parameter_string(char *name, char **field, int *status){
  if (*status!=EXIT_SUCCESS) return;
  char *buf=NULL;
  *status=ape_trad_query_string(name, &buf);
  if (*status!=EXIT_SUCCESS) return;

  *field=strdup(buf);
  free(buf);
}

// read a string parameter into a preallocated buffer
void query_simput_parameter_string_buffer(char *name, char * const field, int buflen, int *status){
  if (*status!=EXIT_SUCCESS) return;
  char *buf=NULL;
  *status=ape_trad_query_string(name, &buf);
  strncpy(field,buf,buflen-1);
  field[buflen-1]='\0';
  free(buf);
}


// for most datatypes we just call the dispatcher routine
void query_simput_parameter_int(char *name, int *field, int *status){
  query_simput_parameter(name, PAR_INT,field,status );
}

void query_simput_parameter_long(char *name, long *field, int *status){
  query_simput_parameter(name, PAR_LONG, field, status );
}

void query_simput_parameter_float(char *name, float *field, int *status){
  query_simput_parameter(name, PAR_FLOAT, field, status );
}

void query_simput_parameter_double(char *name, double *field, int *status){
  query_simput_parameter(name, PAR_DOUBLE, field, status );
}

void query_simput_parameter_bool(char *name, int *field, int *status){
  char dummy;
  query_simput_parameter(name, PAR_BOOL, &dummy, status );
  *field=(int) dummy;
}
