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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef ARF_H
#define ARF_H 1

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns an empty ARF data structure with
    NULL-initialized pointers. */
struct ARF* getARF(int* const status);

/** get a ARF from arrays
    low_energy, high_energy and eff_area all needto have NumberEnergyBins elements. */
struct ARF* getARFfromarrays(long NumberEnergyBins, float low_energy[], float high_energy[], float eff_area[], char* telescope, int* const status);

/** Load an ARF from a response file. */
struct ARF* loadARF(char* filename, int* const status);

/** Destructor for the ARF data structure. Warning: As there is no
    internal destructor for the ARF data structure in the HEASP
    library, the memory allocated by the function ReadARF is released
    by a self-implemented function, which is not guaranteed to work
    properly. */
void freeARF(struct ARF* const arf);


#endif /* ARF_H */
