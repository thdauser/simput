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

#ifndef _POSSTRING_H
#define _POSSTRING_H

#define DEGPOS  0
#define LATPOS  1
#define LONPOS  2
#define HOURPOS 3


double parse_posstring(const char *str, const int type);
char *posstring(double position, const int type);

#endif
