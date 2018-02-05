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


   Copyright 2007-2015 Philippe Peille, IRAP
*/

#ifndef SIMPUTMULTICELL_H
#define SIMPUTMULTICELL_H

#include "multispec.h"

#define TOOLSUB simputmulticell_main
#include "headas_main.c"

struct Parameters {
  /** File name of the output SIMPUT file. */
  char *Simput;

  /* Lower and upper boundary of the generated spectrum [keV].*/
  float Elow;
  float Eup;
  float Estep;

  /* Reference energy band [keV].*/
  float Emin;
  float Emax;

  /** File name of the input ISIS parameter file containing a spectral
       model	. */
  char *ISISFile;


  /* File name for optional preperation script (f. e. to load additional
    models). */
  char *ISISPrep;

  /* File name of the input Xspec spectral model. */
  char *XSPECFile;

   /* File name of the additional Xspec script. */
   char *XSPECPrep;


  // Param to handle the param info to read from the input file
  char *ParamFile;
  char *ParamNames;
  char *ParamInputNames;
  char *ParamsNumValues;
  char *ParamsLogScale;
  char *InputType;
  int FluxIndex;
  int DirectMatch;
  long StartIndex;

  // Usual tool parameters
  int chatter;
  int clobber;
  int history;

  int nbins;
  int logegrid;
};

// Parameter provider function
// this function is called whenever the new cell parameters have to be read in
typedef int (*param_provider) (void *paraminfo, double* para_array, double* ra, double* dec, double* flux, int* status);
// Function to get the min and max parameter values arrays out of the input data
typedef void (*minmax_provider) (void *paraminfo, double* min_array, double* max_array, int* status);
// Function to release memory of the input data
typedef void (*paraminfo_destructor) (void **paraminfo, int* status);


/** Struct to read in the spectral parameter of each cell */
typedef struct{
	void* paraminfo;
	param_provider get_next_cell;
	minmax_provider get_minmax;
	paraminfo_destructor free_param_info;
	int npara;
}ParamReader;

/** Struct for reading in a FITS image as a cells distribution */
typedef struct{
	SimputImg* image;

	/** Number of the current row in the FITS image */
	long row;

	/** Number of parameters sampled */
	int num_par;

	/** Indexes of those parameters in the image */
	int* index_array;

	/** Index of the ra, dec, and flux column */
	int flux_index;

}ImageCellDistrib;

/** Struct for reading in a FITS table as a cells distribution */
typedef struct{
	/** Pointer to fits file*/
	fitsfile* fptr;

	/** Number of the current row in the FITS table */
	long row;

	/** Total number of rows in the FITS table */
	long nrows;

	/** Number of parameters sampled */
	int num_par;

	/** Column indexes of the parameters to read */
	int* col_indexes;

	/** Index of the flux column */
	int col_ra_index,col_dec_index,col_flux_index;

}TableCellDistrib;

void simputmulticell_getpar(struct Parameters* const par,int* status);
void freeParStrings(struct Parameters* par);

#endif /* SIMPUTMULTICELL_H */
