#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "simput.h"


#define N_SRC_CAT_COLUMNS  (10)
#define N_SPEC_COLUMNS     (3)
#define N_LIGHTCUR_COLUMNS (5)

#define SIMPUT_ERROR(msg) fprintf(stderr, "Error in %s: %s!\n", __func__, msg)
#define CHECK_STATUS(a) if (EXIT_SUCCESS!=a) return(a)
#define CHECK_NULL(a,status,msg) if (NULL==a) { SIMPUT_ERROR(msg); status=EXIT_FAILURE; return(status); }


/** Open a SIMPUT source catalog FITS file. Open either an existing
    SIMPUT file or create a new and empty one and open it. The file
    access mode can be either READWRITE or READONLY. */
int simput_open_srcctlg(simput_srcctlg_file** const srcctlg,
			const char* const filename,
			const int mode,
			int* const status)
{
  *srcctlg = (simput_srcctlg_file*)malloc(sizeof(simput_srcctlg_file));
  CHECK_NULL(*srcctlg, *status, "memory allocation failed");

  // Initialize pointer with NULL.
  (*srcctlg)->fptr = NULL;

  // Return value of CFITSIO routines.
  int ret;

  // Check if the requested file already exists.
  int exists;
  ret = fits_file_exists(filename, &exists, status);
  CHECK_STATUS(ret);

  // If it doesn't exist, create a new empty file containing
  // the basic table structure (source catalog table).
  if (1!=exists) { 
    // Maybe have to change the 1 in the comparison above.
    // There might also be other possible values.

    // Create and open a new empty FITS file.
    ret = fits_create_file(&(*srcctlg)->fptr, filename, status);
    CHECK_STATUS(ret);

    // Create header keywords.

    // Create the table structure for the source catalog.
    char *ttype[] = 
      { "SRC_ID", "SRC_NAME", "RA", "DEC", "FLUX", "E_MIN", "E_MAX",
	"SPECTRUM", "LIGHTCUR", "IMAGE" };
    char *tform[] = 
      { "J", "80A", "E", "E", "E", "E", "E", "21A", "21A", "21A" };
    char *tunit[] = 
      { "", "", "deg", "deg", "erg/s/cm^2", "keV", "keV", "", "", "" };
    ret = fits_create_tbl((*srcctlg)->fptr, BINARY_TBL, 0, N_SRC_CAT_COLUMNS, 
			  ttype, tform, tunit, "SRC_CAT", status);
    CHECK_STATUS(ret);

  } else {
    // The file does already exist.

    // Try to open the file via CFITSIO with the requested mode.
    ret = fits_open_file(&(*srcctlg)->fptr, filename, mode, status);
    CHECK_STATUS(ret);
    
    // Note: we do not want to apply a full consistency check here
    // because we also want to be able to access a file, which is
    // not complete yet.
  }
  // END of check if file exists or not.

  // Move the internal HDU pointer of the fitsfile data structure
  // to the right extension containing the source catalog.
  ret = fits_movnam_hdu((*srcctlg)->fptr, BINARY_TBL, "SRC_CAT", 0, status);
  CHECK_STATUS(ret);

  // Determine the column numbers of the essential columns.
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "SRC_ID", &(*srcctlg)->csrc_id, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "SRC_NAME", &(*srcctlg)->csrc_name, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "RA", &(*srcctlg)->cra, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "DEC", &(*srcctlg)->cdec, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "FLUX", &(*srcctlg)->cflux, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "E_MIN", &(*srcctlg)->ce_min, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "E_MAX", &(*srcctlg)->ce_max, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "SPECTRUM", &(*srcctlg)->cspectrum, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "LIGHTCUR", &(*srcctlg)->clightcur, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum((*srcctlg)->fptr, CASEINSEN, "IMAGE", &(*srcctlg)->cimage, status);
  CHECK_STATUS(ret);

  return(*status);
}


/** Close an open SIMPUT source catalog FITS file. */
int simput_close_srcctlg(simput_srcctlg_file** const srcctlg,
			 int* const status) 
{
  if (NULL!=*srcctlg) {
    int ret = fits_close_file((*srcctlg)->fptr, status);
    CHECK_STATUS(ret);
  }

  free(*srcctlg);
  *srcctlg=NULL;

  return(*status);
}


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
		   int* const status) 
{
  int ret;

  // Check if the simputfile pointer points to valid data.
  CHECK_NULL(srcctlg, *status, "simput_srcctlg_file not initialized");
  CHECK_NULL(srcctlg->fptr, *status, "no open FITS file");

  // Determine the current number of lines in the source catalog.
  long nrows;
  ret = fits_get_num_rows(srcctlg->fptr, &nrows, status);
  CHECK_STATUS(ret);

  // Add a new line at the end of the catalog.
  ret = fits_insert_rows(srcctlg->fptr, nrows++, 1, status);
  CHECK_STATUS(ret);

  // Insert the source data to the table.
  ret = fits_write_col(srcctlg->fptr, TLONG, srcctlg->csrc_id, 
		       nrows, 1, 1, &src_id, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TSTRING, srcctlg->csrc_name, 
		       nrows, 1, 1, src_name, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TFLOAT, srcctlg->cra, 
		       nrows, 1, 1, &ra, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TFLOAT, srcctlg->cdec, 
		       nrows, 1, 1, &dec, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TFLOAT, srcctlg->cflux, 
		       nrows, 1, 1, &flux, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TFLOAT, srcctlg->ce_min, 
		       nrows, 1, 1, &e_min, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TFLOAT, srcctlg->ce_max, 
		       nrows, 1, 1, &e_max, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TSTRING, srcctlg->cspectrum, 
		       nrows, 1, 1, spectrum, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TSTRING, srcctlg->clightcur, 
		       nrows, 1, 1, lightcur, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(srcctlg->fptr, TSTRING, srcctlg->cimage, 
		       nrows, 1, 1, image, status);
  CHECK_STATUS(ret);

  return(*status);
}




/** Store a flux density spectrum in the given file. If the file
    already exists, append a new table. */
int simput_store_spectrum(const char* const filename,
			  const long nbins,
			  float* const e_min,
			  float* const e_max,
			  float* const flux,
			  float phase,
			  int* const status)
{
  fitsfile* fptr=NULL;

  // Check if the file already exists.
  int exists;
  int ret = fits_file_exists(filename, &exists, status);
  CHECK_STATUS(ret);

  // If no, create a new file.
  if (1!=exists) { 
    // Create and open a new empty FITS file.
    ret = fits_create_file(&fptr, filename, status);
    CHECK_STATUS(ret);
  } else {
    // The file does already exist, so just open it.
    ret = fits_open_file(&fptr, filename, READWRITE, status);
    CHECK_STATUS(ret);
  }
  // END of check, whether the file already exists or not.

  // Create a new table for the spectrum.
  char *ttype[] = { "E_MIN", "E_MAX", "FLUX" };
  char *tform[] = { "E", "E", "E" };
  char *tunit[] = { "keV", "keV", "erg/s/cm^2/keV" };
  ret = fits_create_tbl(fptr, BINARY_TBL, 0, N_SPEC_COLUMNS, 
			ttype, tform, tunit, "SPECTRUM", status);
  CHECK_STATUS(ret);

  // Insert header keywords.
  ret = fits_update_key(fptr, TFLOAT, "PHASE", &phase, 
			"phase for which the spectrum is valid", status);
  CHECK_STATUS(ret);

  // Determine the column numbers of the essential columns.
  int ce_min, ce_max, cflux;
  ret = fits_get_colnum(fptr, CASEINSEN, "E_MIN", &ce_min, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum(fptr, CASEINSEN, "E_MAX", &ce_max, status);
  CHECK_STATUS(ret);
  ret = fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
  CHECK_STATUS(ret);

  // Store the spectrum in the table.
  ret = fits_insert_rows(fptr, 1, nbins, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(fptr, TFLOAT, ce_min, 
		       1, 1, nbins, e_min, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(fptr, TFLOAT, ce_max, 
		       1, 1, nbins, e_max, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(fptr, TFLOAT, cflux, 
		       1, 1, nbins, flux, status);
  CHECK_STATUS(ret);

  // Close the file.
  ret = fits_close_file(fptr, status);
  CHECK_STATUS(ret);

  return(*status);
}



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
			  int* const status)
{
  fitsfile* fptr=NULL;

  // Check if the file already exists.
  int exists;
  int ret = fits_file_exists(filename, &exists, status);
  CHECK_STATUS(ret);

  // If no, create a new file.
  if (1!=exists) { 
    // Create and open a new empty FITS file.
    ret = fits_create_file(&fptr, filename, status);
    CHECK_STATUS(ret);
  } else {
    // The file does already exist, so just open it.
    ret = fits_open_file(&fptr, filename, READWRITE, status);
    CHECK_STATUS(ret);
  }
  // END of check, whether the file already exists or not.

  // Create a new table for the spectrum.
  char *ttype[N_LIGHTCUR_COLUMNS];
  char *tform[N_LIGHTCUR_COLUMNS];
  char *tunit[N_LIGHTCUR_COLUMNS];
  int ncolumns=0;
  if (NULL!=time) {
    ttype[ncolumns] = "TIME";
    tform[ncolumns] = "D";
    tunit[ncolumns] = "s";
    ncolumns++;
  }
  if (NULL!=phase) {
    ttype[ncolumns] = "PHASE";
    tform[ncolumns] = "E";
    tunit[ncolumns] = "";
    ncolumns++;
  }
  ttype[ncolumns] = "FLUX";
  tform[ncolumns] = "E";
  tunit[ncolumns] = "erg/s/cm^2";
  ncolumns++;
  if (NULL!=pol_frac) {
    ttype[ncolumns] = "POL_FRAC";
    tform[ncolumns] = "E";
    tunit[ncolumns] = "";
    ncolumns++;
  }
  if (NULL!=phase) {
    ttype[ncolumns] = "POL_DIR";
    tform[ncolumns] = "E";
    tunit[ncolumns] = "radians";
    ncolumns++;
  }
   ret = fits_create_tbl(fptr, BINARY_TBL, 0, ncolumns, 
			ttype, tform, tunit, "LIGHTCUR", status);
  CHECK_STATUS(ret);

  // Insert header keywords.
  ret = fits_update_key(fptr, TFLOAT, "E_MIN", &e_min, 
			"lower value of the reference energy band", status);
  CHECK_STATUS(ret);
  ret = fits_update_key(fptr, TFLOAT, "E_MAX", &e_max, 
			"upper value of the reference energy band", status);
  CHECK_STATUS(ret);

  // Determine the column numbers of the essential columns and store 
  // the light curve in the table.
  ret = fits_insert_rows(fptr, 1, nbins, status);
  CHECK_STATUS(ret);
  int column;
  if (NULL!=time) {
    ret = fits_get_colnum(fptr, CASEINSEN, "TIME", &column, status);
    CHECK_STATUS(ret);
    ret = fits_write_col(fptr, TDOUBLE, column, 1, 1, nbins, time, status);
    CHECK_STATUS(ret);
  }
  if (NULL!=phase) {
    ret = fits_get_colnum(fptr, CASEINSEN, "PHASE", &column, status);
    CHECK_STATUS(ret);
    ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, phase, status);
    CHECK_STATUS(ret);
  }
  ret = fits_get_colnum(fptr, CASEINSEN, "FLUX", &column, status);
  CHECK_STATUS(ret);
  ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, flux, status);
  CHECK_STATUS(ret);
  if (NULL!=pol_frac) {
    ret = fits_get_colnum(fptr, CASEINSEN, "POL_FRAC", &column, status);
    CHECK_STATUS(ret);
    ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, pol_frac, status);
    CHECK_STATUS(ret);
  }
  if (NULL!=pol_dir) {
    ret = fits_get_colnum(fptr, CASEINSEN, "POL_DIR", &column, status);
    CHECK_STATUS(ret);
    ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, pol_dir, status);
    CHECK_STATUS(ret);
  }

  // Close the file.
  ret = fits_close_file(fptr, status);
  CHECK_STATUS(ret);

  return(*status);
}
