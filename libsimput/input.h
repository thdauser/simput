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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef INPUT_H
#define INPUT_H (1)

#include "simput.h"
#include "common.h"
#include "ape/ape_trad.h"

void query_simput_parameter_file_name(char *name, char *field, int *status);
void query_simput_parameter_string(char *name, char *field, int *status);
void query_simput_parameter_int(char *name, int *field, int *status);
void query_simput_parameter_float(char *name, float *field, int *status);
void query_simput_parameter_bool(char *name, char *field, int *status);

#endif /* INPUT_H */
