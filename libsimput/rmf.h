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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef RMF_H
#define RMF_H 1

#include "arf.h"
#include <unistd.h>

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


/** Constructor. Returns an empty RMF data structure with
    pointers initialized with NULL. */
struct RMF* getRMF(int* const status);

/** Load an RMF matrix and the corresponding EBOUNDS from a response
    file. */
struct RMF* loadRMF(char* const filename, int* const status);

/** Load an RMF matrix and the corresponding EBOUNDS from a response
    file and check for normalization. */
struct RMF* loadNormalizedRMF(char* const filename, int* const status);

/** Load an RSP matrix and the corresponding EBOUNDS from a response
    file and split the RSP into an ARF and an RMF. */
void loadArfRmfFromRsp(char* const filename,
		       struct ARF** arf,
		       struct RMF** rmf,
		       int* const status);

/** Check the validity of a FITS response file, using ftchkrmf
    from the HEASoft FTOOLS. */
void checkRMFValidity(char* const filename, int* const status);

/** Check if the RMF file contains matrix rows with a sum of more than 1.
    In that case the RSP probably also contains the mirror ARF, what should
    not be the case for this simulation. Row sums with a value of less than
    1 should actually also not be used, but can be handled by the simulation. */
void checkRMFNormalization(const struct RMF* const rmf, int* const status);

/** Destructor for the RMF data structure. Warning: As there is no
    internal destructor for the RMF data structure in the HEASP
    library, the memory allocated by the function ReadRMFMatrix is not
    realeased. */
void freeRMF(struct RMF* const rmf);

/** Returns a randomly selected channel for a photon of the given
    input energy. The channel number can start at 0, 1, or any other
    positive value defined in the RMF. The return value '-1' means
    that no appropriate channel could be determined. */
void returnRMFChannel(struct RMF *rmf,
		      const float energy,
		      long* const channel);

/** Load the EBOUNDS extension from an RMF or RSP file. */
void loadEbounds(struct RMF* rmf, char* const filename, int* const status);

/** Determines the PHA channel corresponding to a given energy
    according to the EBOUNDS table of the detector response. The
    routine performs a binary search to obtain the PHA channel the
    specified energy lies within. The energy has to be given in the
    same units as the EBOUNDS are, which usually is [keV]. Note that
    the routine is NOT doing an RMF randomization of the measured
    channel. In case the channel cannot be determined, because there
    is no RMF available or the desired energy lies outside the covered
    range, a negative value is returned. */
long getEBOUNDSChannel(const float energy, const struct RMF* const rmf);

/** Determine the signal corresponding to a particular PHA channel
    according to the EBOUNDS table. The input channel must have the
    same offset as in the EBOUNDS table. I.e. if the first channel in
    the EBOUNDS has the number 1, the numbering starts at 1. If the
    first channel has the number 0, the numbering starts at 0.  The
    returned energy is given in the same units as the EBOUNDS,
    (usually [keV]). */
void getEBOUNDSEnergyLoHi(const long channel,
			  const struct RMF* const rmf,
			  float* const lo,
			  float* const hi,
			  int* const status);


#endif /* RMF_H */
