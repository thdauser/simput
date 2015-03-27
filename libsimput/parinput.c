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
*/


#include "parinput.h"
#include "common.h"
#include "simput.h"
#include "ape/ape_trad.h"

#define PAR_MAXBUFSIZE 1024

void* query_simput_parameter(char* name, const int type, int* status ){

  switch (type) {
  case PAR_STRING:
  case PAR_FILE:
    {
      // JW: great danger of buffer overruns, 
      // but I do not know what to do about this...
      //      char *buffer=malloc(PAR_MAXBUFSIZE*sizeof(char));  
      //CHECK_NULL_RET(buffer, *status,"memory allocation failed",NULL);
      char *buffer=NULL;
      if (type==PAR_STRING) {
	*status=ape_trad_query_string(name, &buffer);
      } else {
	*status=ape_trad_query_file_name(name, &buffer);
      }
      char *sbuf;
      sbuf=strndup(buffer,PAR_MAXBUFSIZE-1);
      free(buffer);
      return sbuf;
    }
  case PAR_FLOAT:
    {
      float *fbuf=malloc(sizeof(float));
      CHECK_NULL_RET(fbuf, *status,"memory allocation failed",NULL);
      *status=ape_trad_query_float(name,fbuf);
      return fbuf;
    }
  case PAR_DOUBLE:
    {
      double *dbuf=malloc(sizeof(double));
      CHECK_NULL_RET(dbuf, *status,"memory allocation failed",NULL);
      *status=ape_trad_query_double(name,dbuf);
      return dbuf;
    }
  case PAR_INT:
    {
      int *ibuf=malloc(sizeof(int));
      CHECK_NULL_RET(ibuf, *status,"memory allocation failed",NULL);
      *status = ape_trad_query_int(name,ibuf);
      return ibuf;
    }
  case PAR_BOOL:
    {
      char *cbuf=malloc(sizeof(char));
      CHECK_NULL_RET(cbuf, *status,"memory allocation failed",NULL);
      *status = ape_trad_query_bool(name,cbuf);
      return cbuf;
    }
  }

  *status = EXIT_FAILURE;
  printf("Failed to read parameter %s from the command line.",name);
  return NULL;
}

static void simput_getpar_empty_file_name(char *file_name){
	if ((0==strcmp(file_name, "none"))|| (0==strcmp(file_name, "NONE"))) {
		strcpy(file_name, "");
	}
}

void query_simput_parameter_file_name(char *name, char **field, int *status){
  char *sbuf = (char *) query_simput_parameter(name, PAR_FILE, status );
  CHECK_NULL_VOID(sbuf, *status,"reading simput parameter failed"); 
  if (*status == EXIT_SUCCESS){
    simput_getpar_empty_file_name(sbuf);
    *field=strdup(sbuf);
  }
  free(sbuf);
}

void query_simput_parameter_string(char *name, char **field, int *status){
  char *sbuf = (char *) query_simput_parameter(name, PAR_STRING, status );
  if (*status == EXIT_SUCCESS) {
    *field=strdup(sbuf);
  }
  free(sbuf);
}

void query_simput_parameter_int(char *name, int *field, int *status){
  int *ibuf = (int *) query_simput_parameter(name, PAR_INT,status );
  if (*status == EXIT_SUCCESS){
    *field =  *ibuf;
  }
  free(ibuf);
}

void query_simput_parameter_long(char *name, long *field, int *status){
  long *lbuf = (long *) query_simput_parameter(name, PAR_LONG, status );
  if (*status == EXIT_SUCCESS){
    *field =  *lbuf;
  }
  free(lbuf);
}

void query_simput_parameter_float(char *name, float *field, int *status){
  float *fbuf = (float *) query_simput_parameter(name, PAR_FLOAT, status );
  if (*status == EXIT_SUCCESS){
    *field =  *fbuf;
  }
  free(fbuf);
}

void query_simput_parameter_double(char *name, double *field, int *status){
  double *fbuf = (double *) query_simput_parameter(name, PAR_DOUBLE, status );
  if (*status == EXIT_SUCCESS) {
    *field =  *fbuf;
  }
  free(fbuf);
}

void query_simput_parameter_bool(char *name, char *field, int *status){
  char *bbuf = (char *) query_simput_parameter(name, PAR_BOOL, status );
  if (*status == EXIT_SUCCESS) {
    *field =  *bbuf;
  }
  free(bbuf);
}
