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


   Copyright 2007-2014 Thomas Dauser, FAU
*/

#ifndef PARINPUT_H
#define PARINPUT_H


// read parameter name as a file name. Field is a char *. The routine
// allocates memory for the parameter. field can be free'd if the parameter is
// not needed anymore
void query_simput_parameter_file_name(char *name, char **field, int *status);

// read parameter name as a string. Field is a char *. The routine
// allocates memory for the parameter. field can be free'd if the parameter is
// not needed anymore
void query_simput_parameter_string(char *name, char **field, int *status);

// read parameter name as a file name. Field is a pointer to a pre-allocated buffer
// of length buflen which will be overwritten with the null-terminated file name is written. If the parameter
// is longer than buflen, only buflen-1 characters are written. field is guaranteed
// to be null-terminated.
// files named "none" (capitals are also allowed) are returned as strings of length 0
void query_simput_parameter_file_name_buffer(char *name, char * const field, int buflen, int *status);

// read parameter name as a string. Field is a pointer to a pre-allocated buffer
// of length buflen which will be overwritten with the null-terminated string. If the parameter
// is longer than buflen, only buflen-1 characters are written. field is guaranteed
// to be null-terminated.
void query_simput_parameter_string_buffer(char *name, char * const field, int buflen, int *status);

// read parameter name as an int, return in field
void query_simput_parameter_int(char *name, int *field, int *status);

// read parameter name as a long, return in field
void query_simput_parameter_long(char *name, long *field, int *status);

// read parameter name as a float, return in field
void query_simput_parameter_float(char *name, float *field, int *status);

// read parameter name as a double, return in field
void query_simput_parameter_double(char *name, double *field, int *status);

// read parameter name as a boolean, return in field
void query_simput_parameter_bool(char *name, int *field, int *status);

#endif /* PARINPUT_H */
