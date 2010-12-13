#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Add a new line with source data to a source catalog file. If the
    file or the source catalog extension does not exists yet, create a
    new file with the appropriate extension. The ID (src_id) must be
    unique. If a source with the same ID as given by the parameter
    already exists in the source catalog, an error is returned. */
void simput_add_src(const char* const filename,
		    long src_id,
		    char* src_name,
		    float ra,
		    float dec,
		    float flux,
		    float e_min,
		    float e_max,
		    char* spectrum,
		    char* lightcur,
		    char* image,
		    int* const status);

/** Store a flux density spectrum in the given file. If the file
    already exists, append a new table. */
void simput_store_spectrum(const char* const filename,
			   char* const extname,
			   const long nbins,
			   float* const e_min,
			   float* const e_max,
			   float* const flux,
			   float phase,
			   int* extver,
			   int* const status);

/** Store a flux density spectrum in the given file. If the file
    already exists, append a new table. */
void simput_store_lightcur(const char* const filename,
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

/** Insert a reference to a spectrum into a specific source
    description in the given source catalog. */
void simput_add_spectrum(const char* const srcctlg_filename,
			 const long src_id,
			 const char* const spec_filename,
			 int* const status);

#endif /* SIMPUT_H */
