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

#include "arf.h"
#include "common.h"


struct ARF* getARF(int* const status)
{
  struct ARF* arf=(struct ARF*)malloc(sizeof(struct ARF));
  CHECK_NULL_RET(arf, *status, "memory allocation for ARF failed", arf);

  // Initialize.
  arf->NumberEnergyBins=0;
  arf->LowEnergy       =NULL;
  arf->HighEnergy      =NULL;
  arf->EffArea         =NULL;
  arf->ARFVersion[0]='\0';
  arf->Telescope[0]='\0';
  arf->Instrument[0]='\0';
  arf->Detector[0]='\0';
  arf->Filter[0]='\0';
  arf->ARFExtensionName[0]='\0';
  
  return(arf);
}


struct ARF* getARFfromarrays(long NumberEnergyBins, float low_energy[], float high_energy[], float eff_area[], char* telescope, int* const status)
{
  struct ARF* arf=getARF(status);
  CHECK_STATUS_RET(*status, arf);

  arf->NumberEnergyBins=NumberEnergyBins;

  arf->LowEnergy = (float*)malloc(arf->NumberEnergyBins*sizeof(float));
  CHECK_NULL_RET(arf->LowEnergy, *status, "memory allocation for ARF failed", arf);
  arf->LowEnergy = low_energy;

  arf->HighEnergy = (float*)malloc(arf->NumberEnergyBins*sizeof(float));
  CHECK_NULL_RET(arf->HighEnergy, *status, "memory allocation for ARF failed", arf);
  arf->HighEnergy = high_energy;

  arf->EffArea = (float*)malloc(arf->NumberEnergyBins*sizeof(float));
  CHECK_NULL_RET(arf->EffArea, *status, "memory allocation for ARF failed", arf);
  arf->EffArea = eff_area;

  strcpy(arf->ARFVersion,"1.0");          /* SPECRESP extension format version */
  strcpy(arf->Telescope, telescope);
  strcpy(arf->Instrument, "any");
  strcpy(arf->Detector, "any");
  strcpy(arf->Filter, "none");
  strcpy(arf->ARFExtensionName, "none");
  return arf;
}


struct ARF* loadARF(char* filename, int* const status) 
{
  struct ARF* arf=getARF(status);
  CHECK_STATUS_RET(*status, arf);

  // Load the ARF from the FITS file using the HEAdas ARF access
  // routines (part of libhdsp).
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, arf);
  
  // Read the 'SPECRESP'.
  *status=ReadARF(fptr, 0, arf);
  CHECK_STATUS_RET(*status, arf);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, arf);

  // Print some information:
  headas_chat(5, "ARF loaded with %ld energy bins\n",
	      arf->NumberEnergyBins);

  int ii;
  for (ii=0; ii<arf->NumberEnergyBins; ii++){
	  if(! (arf->EffArea[ii] >=0)){
		  headas_chat(5, " *** Warning: found value %e in ARF for bin %i, setting it to 0.0\n",
				  arf->EffArea,ii);
		  arf->EffArea[ii] = 0.0;
	  }
  }

  return(arf);
}


void freeARF(struct ARF* const arf) 
{
  if (NULL!=arf) {
    if (NULL!=arf->LowEnergy) {
      free(arf->LowEnergy);
    }
    if (NULL!=arf->HighEnergy) {
      free(arf->HighEnergy);
    }
    if (NULL!=arf->EffArea) {
      free(arf->EffArea);
    }
    free(arf);
  }
}

