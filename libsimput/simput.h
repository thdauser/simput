#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;
  
  /** Column number of the essential columns. */
  int csrc_id, csrc_name, cra, cdec, cflux, ce_min, ce_max;
  int cspectrum, clightcur, cimage;

} simput_srcctlg_file;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

// TODO Change this -> Instead of open and close routine just provide one routine to add a source to a source catalog (and create a new source catalog if there is none in the FITS file).

/** Open a SIMPUT source catalog FITS file. Open either an existing
    SIMPUT file or create a new and empty one and open it. The file
    access mode can be either READWRITE or READONLY. */
int simput_open_srcctlg(simput_srcctlg_file** const srcctlg,
			const char* const filename,
			const int mode,
			int* const status);

/** Close an open SIMPUT source catalog FITS file. */
int simput_close_srcctlg(simput_srcctlg_file** const srcctlg,
			 int* const status);

/** Add a new line with source data to a source catalog. */
int simput_add_src(simput_srcctlg_file* const srcctlg,
		   long src_id,
		   char* const src_name,
		   float ra,
		   float dec,
		   float flux,
		   float e_min,
		   float e_max,
		   char* const spectrum,
		   char* const lightcur,
		   char* const image,
		   int* const status);

/** Store a flux density spectrum in the given file. If the file
    already exists, append a new table. */
int simput_store_spectrum(const char* const filename,
			  const long nbins,
			  float* const e_min,
			  float* const e_max,
			  float* const flux,
			  float phase,
			  int* hdunum,
			  int* const status);

/** Store a flux density spectrum in the given file. If the file
    already exists, append a new table. */
int simput_store_lightcur(const char* const filename,
			  const long nbins,
			  double* const time,
			  float* const phase,
			  float* const flux,
			  float* const pol_frac,
			  float* const pol_dir,
			  float e_min,
			  float e_max,
			  int* hdunum,
			  int* const status);


#endif /* SIMPUT_H */
