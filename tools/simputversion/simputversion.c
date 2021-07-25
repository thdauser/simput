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


   Copyright 2007-2018 Remeis-Sternwarte
*/

#include <simputconfig.h>

#include <stdio.h>
#include <stdlib.h>

int main()
{
  printf("SIMPUT version %s\n",PACKAGE_VERSION);
  printf("Compiled %s, %s\n",__DATE__,__TIME__);
  exit(0);
}
