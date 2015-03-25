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


   Copyright 2007-2014 Thomas Dauser FAU
*/


#include "parinput.h"

static void* query_simput_parameter(char* name, char* type, int* status ){

	if (0==strcmp(type,"string") || 0==strcmp(type,"file") ){
		char *sbuf=NULL; // malloc here as well?
		if (0==strcmp(type,"string"))
			*status=ape_trad_query_string(name, &sbuf);
		else
			*status=ape_trad_query_file_name(name, &sbuf);
		return (char*) sbuf;
	} else if (0==strcmp(type,"float")) {
		float *fbuf=malloc(sizeof(float));
		CHECK_NULL_RET(fbuf, *status,"memory allocation failed",NULL);
		*status=ape_trad_query_float(name,fbuf);
		return (float *)fbuf;
	} else if (0==strcmp(type,"double")) {
		float *fbuf=malloc(sizeof(double));
		CHECK_NULL_RET(fbuf, *status,"memory allocation failed",NULL);
		*status=ape_trad_query_float(name,fbuf);
		return (double *)fbuf;
	} else if (0==strcmp(type,"int") ) {
		int *ibuf=malloc(sizeof(int));
		CHECK_NULL_RET(ibuf, *status,"memory allocation failed",NULL);
		*status = ape_trad_query_int(name,ibuf);
		return (int*) ibuf;
	} else if (0==strcmp(type,"bool") ) {
		char *cbuf=malloc(sizeof(char));
		CHECK_NULL_RET(cbuf, *status,"memory allocation failed",NULL);
		*status = ape_trad_query_bool(name,cbuf);
		return (char*) cbuf;
	} else {
		*status = EXIT_FAILURE;
		printf("Failed to read parameter %s from the command line.",name);
		return NULL;
	}

}

static void simput_getpar_empty_file_name(char *file_name){
	if ((0==strcmp(file_name, "none"))|| (0==strcmp(file_name, "NONE"))) {
		strcpy(file_name, "");

	}
}

void query_simput_parameter_file_name(char *name, char *field, int *status){
	char *sbuf = (char *) query_simput_parameter(name, "file", status );
	CHECK_NULL_VOID(sbuf, *status,"reading simput parameter failed"); // need to free a NULL pointer?
	if (*status == EXIT_SUCCESS){
		strcpy(field,sbuf);
		simput_getpar_empty_file_name(sbuf);
	}
	free(sbuf);
}

void query_simput_parameter_string(char *name, char *field, int *status){
	char *sbuf = (char *) query_simput_parameter(name, "string", status );
	if (*status == EXIT_SUCCESS) {
		strcpy(field,sbuf);
	}
	free(sbuf);
}

void query_simput_parameter_int(char *name, int *field, int *status){
	int *ibuf = (int *) query_simput_parameter(name, "int", status );
	if (*status == EXIT_SUCCESS){
		*field =  *ibuf;
	}
	free(ibuf); // do we need this?
}

void query_simput_parameter_long(char *name, long *field, int *status){
	int *ibuf = (long *) query_simput_parameter(name, "long", status );
	if (*status == EXIT_SUCCESS){
		*field =  *ibuf;
	}
	free(ibuf); // do we need this?
}

void query_simput_parameter_float(char *name, float *field, int *status){
	// do we a free here need this?
	float *fbuf = (float *) query_simput_parameter(name, "float", status );
	if (*status == EXIT_SUCCESS){
		*field =  *fbuf;
	}
	free(fbuf);
}

void query_simput_parameter_double(char *name, double *field, int *status){
	// do we a free here need this?
	double *fbuf = (double *) query_simput_parameter(name, "double", status );
	if (*status == EXIT_SUCCESS)
		*field =  *fbuf;
	free(fbuf);
}

void query_simput_parameter_bool(char *name, char *field, int *status){
	char *bbuf = (char *) query_simput_parameter(name, "bool", status );
	if (*status == EXIT_SUCCESS) {
		*field =  *bbuf;
	}
	free(bbuf); // do we need this?
}
