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

#define MAXMSG (1024)

#define N_SRC_CAT_COLUMNS  (10)
#define N_SPEC_COLUMNS     (3)
#define N_LIGHTCUR_COLUMNS (5)

/////////////////////////////////////////////////////////////////
// Structures.
/////////////////////////////////////////////////////////////////

typedef struct {

	fitsfile* fptr;

	/** Column numbers. */
	int csrc_id, csrc_name, cra, cdec, cflux, ce_min, ce_max;
	int cspectrum, clightcur, cimage;

} SIMPUT_SrcCtlg;

/////////////////////////////////////////////////////////////////
// Macros.
/////////////////////////////////////////////////////////////////

#define SIMPUT_ERROR(msg) \
  fprintf(stderr, "Error in %s: %s!\n", __func__, msg)

#define CHECK_STATUS_RET(a,b) \
  if (EXIT_SUCCESS!=a) return(b)

#define CHECK_STATUS(a) \
  CHECK_STATUS_RET(a,a)

#define CHECK_NULL(a,status,msg,ret) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return(ret);\
  }

/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Determine the line number of a particular source in a source
 catalog. The source is identified by its source ID. If no source
 with this ID is contained in the catalog, the return value is
 zero. */
static long simput_get_src_linenum(fitsfile* const fptr, const long src_id,
		int* const status);

/** Create an extension identifier string that can be used to uniquely
 identify a particular FITS extension like, e.g., a binary table
 containing a source spectrum. */
static void simput_ext_id(char* const id, const char* const filename,
		const char* const extname, const int extver);

/////////////////////////////////////////////////////////////////
// Function Definitions.
/////////////////////////////////////////////////////////////////

int simput_add_src(const char* const filename, long src_id, char* src_name,
		float ra, float dec, float flux, float e_min, float e_max,
		char* spectrum, char* lightcur, char* image, int* const status) {
	fitsfile* fptr = NULL;

	// Check if the file already exists.
	int exists;
	int ret = fits_file_exists(filename, &exists, status);
	CHECK_STATUS(ret);

	// If no, create a new file.
	if (1 != exists) {
		// Create and open a new empty FITS file.
		ret = fits_create_file(&fptr, filename, status);
		CHECK_STATUS(ret);
	} else {
		// The file does already exist, so just open it.
		ret = fits_open_file(&fptr, filename, READWRITE, status);
		CHECK_STATUS(ret);
	}
	// END of check, whether the file already exists or not.

	// Check if a source catalog extension exists.
	// Try to move the internal HDU pointer of the fitsfile data structure
	// to the right extension containing the source catalog.
	fits_write_errmark();
	int temp_status = EXIT_SUCCESS; // We have to use another status variable here.
	fits_movnam_hdu(fptr, BINARY_TBL, "SRC_CAT", 0, &temp_status);
	fits_clear_errmark();

	if (BAD_HDU_NUM == temp_status) {
		// Create the table structure for the source catalog.
		char *ttype[] = { "SRC_ID", "SRC_NAME", "RA", "DEC", "FLUX", "E_MIN",
				"E_MAX", "SPECTRUM", "LIGHTCUR", "IMAGE" };
		char *tform[] = { "J", "20A", "E", "E", "E", "E", "E", "21A", "21A",
				"21A" };
		char *tunit[] = { "", "", "deg", "deg", "erg/s/cm^2", "keV", "keV", "",
				"", "" };
		ret = fits_create_tbl(fptr, BINARY_TBL, 0, N_SRC_CAT_COLUMNS, ttype,
				tform, tunit, "SRC_CAT", status);
		CHECK_STATUS(ret);

		// Insert header keywords.
		// ...

	}

	// Note: we do not want to apply a full consistency check here
	// because we also want to be able to access a file, which is
	// not complete yet.

	// Determine the column numbers.
	int csrc_id, csrc_name, cra, cdec, cflux, ce_min, ce_max;
	int cspectrum, clightcur, cimage;
	ret = fits_get_colnum(fptr, CASEINSEN, "SRC_ID", &csrc_id, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "SRC_NAME", &csrc_name, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "RA", &cra, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "DEC", &cdec, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "E_MIN", &ce_min, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "E_MAX", &ce_max, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "SPECTRUM", &cspectrum, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "LIGHTCUR", &clightcur, status);
	CHECK_STATUS(ret);
	ret = fits_get_colnum(fptr, CASEINSEN, "IMAGE", &cimage, status);
	CHECK_STATUS(ret);

	// Store the data:

	// Check if a source of this ID is already contained in the catalog.
	long row = simput_get_src_linenum(fptr, src_id, status);
	CHECK_STATUS(ret);

	// If the source is already contained in the catalog, create an
	// error message.
	if (row > 0) {
		// ERRMSG
		*status = EXIT_FAILURE;
		CHECK_STATUS(*status);
	}

	// Determine the current number of lines in the source catalog.
	long nrows;
	ret = fits_get_num_rows(fptr, &nrows, status);
	CHECK_STATUS(ret);

	// Insert the source data at the end of the table.
	nrows++;
	ret = fits_write_col(fptr, TLONG, csrc_id, nrows, 1, 1, &src_id, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TSTRING, csrc_name, nrows, 1, 1, &src_name,
			status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, cra, nrows, 1, 1, &ra, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, cdec, nrows, 1, 1, &dec, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, cflux, nrows, 1, 1, &flux, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, ce_min, nrows, 1, 1, &e_min, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, ce_max, nrows, 1, 1, &e_max, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TSTRING, cspectrum, nrows, 1, 1, &spectrum,
			status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TSTRING, clightcur, nrows, 1, 1, &lightcur,
			status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TSTRING, cimage, nrows, 1, 1, &image, status);
	CHECK_STATUS(ret);

	// Close the file.
	ret = fits_close_file(fptr, status);
	CHECK_STATUS(ret);

	return (*status);
}

static long simput_get_src_linenum(fitsfile* const fptr, const long src_id,
		int* const status) {
	// Determine the column of the source ID.
	int csrc_id;
	int ret = fits_get_colnum(fptr, CASEINSEN, "SRC_ID", &csrc_id, status);
	CHECK_STATUS(ret);

	// Determine the number of lines in the source catalog.
	long nrows;
	ret = fits_get_num_rows(fptr, &nrows, status);
	CHECK_STATUS(ret);

	long row, id, nulval = 0;
	int anynul = 0;
	for (row = 1; row <= nrows; row++) {
		ret = fits_read_col(fptr, TLONG, csrc_id, row, 1, 1, &nulval, &id,
				&anynul, status);
		CHECK_STATUS(ret);

		if (src_id == id)
			return (row);
	}

	return (0);
}

int simput_store_spectrum(const char* const filename, char* const extname,
		const long nbins, float* const e_min, float* const e_max,
		float* const flux, float phase, int* extver, int* const status) {
	fitsfile* fptr = NULL;

	// Check if the file already exists.
	int exists;
	int ret = fits_file_exists(filename, &exists, status);
	CHECK_STATUS(ret);

	// If no, create a new file.
	if (1 != exists) {
		// Create and open a new empty FITS file.
		ret = fits_create_file(&fptr, filename, status);
		CHECK_STATUS(ret);
	} else {
		// The file does already exist, so just open it.
		ret = fits_open_file(&fptr, filename, READWRITE, status);
		CHECK_STATUS(ret);
	}
	// END of check, whether the file already exists or not.

	// Check if a particular name is specified for the new extension.
	char extname2[MAXMSG];
	if (NULL == extname) {
		// If no explicit exname is given, use the default value "SPECTRUM".
		strcpy(extname2, "SPECTRUM");
	}

	// Check if the FITS file already contains an extension with
	// the same EXTNAME. Try to move the internal HDU pointer of the
	// fitsfile data structure to such an extension.
	*extver = 1;
	fits_write_errmark();
	int temp_status;
	do {
		temp_status = EXIT_SUCCESS; // We have to use another status variable here.
		fits_movnam_hdu(fptr, BINARY_TBL, extname2, *extver, &temp_status);
	} while (BAD_HDU_NUM == temp_status);
	fits_clear_errmark();

	// Create a new table for the spectrum with a unique
	// EXTNAME and EXTVER combination.
	char *ttype[] = { "E_MIN", "E_MAX", "FLUX" };
	char *tform[] = { "E", "E", "E" };
	char *tunit[] = { "keV", "keV", "erg/s/cm^2/keV" };
	ret = fits_create_tbl(fptr, BINARY_TBL, 0, N_SPEC_COLUMNS, ttype, tform,
			tunit, extname2, status);
	CHECK_STATUS(ret);
	ret = fits_update_key(fptr, TINT, "EXTVER", extver, "extension identifier",
			status);
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
	ret = fits_write_col(fptr, TFLOAT, ce_min, 1, 1, nbins, e_min, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, ce_max, 1, 1, nbins, e_max, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, cflux, 1, 1, nbins, flux, status);
	CHECK_STATUS(ret);

	// Close the file.
	ret = fits_close_file(fptr, status);
	CHECK_STATUS(ret);

	return (*status);
}

int simput_store_lightcur(const char* const filename, const long nbins,
		double* const time, float* const phase, float* const flux,
		float* const pol_frac, float* const pol_dir, float e_min, float e_max,
		int* hdunum, int* const status) {
	fitsfile* fptr = NULL;

	// Check if the file already exists.
	int exists;
	int ret = fits_file_exists(filename, &exists, status);
	CHECK_STATUS(ret);

	// If no, create a new file.
	if (1 != exists) {
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
	int ncolumns = 0;
	if (NULL != time) {
		ttype[ncolumns] = "TIME";
		tform[ncolumns] = "D";
		tunit[ncolumns] = "s";
		ncolumns++;
	}
	if (NULL != phase) {
		ttype[ncolumns] = "PHASE";
		tform[ncolumns] = "E";
		tunit[ncolumns] = "";
		ncolumns++;
	}
	ttype[ncolumns] = "FLUX";
	tform[ncolumns] = "E";
	tunit[ncolumns] = "erg/s/cm^2";
	ncolumns++;
	if (NULL != pol_frac) {
		ttype[ncolumns] = "POL_FRAC";
		tform[ncolumns] = "E";
		tunit[ncolumns] = "";
		ncolumns++;
	}
	if (NULL != phase) {
		ttype[ncolumns] = "POL_DIR";
		tform[ncolumns] = "E";
		tunit[ncolumns] = "radians";
		ncolumns++;
	}
	ret = fits_create_tbl(fptr, BINARY_TBL, 0, ncolumns, ttype, tform, tunit,
			"LIGHTCUR", status);
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
	int column;
	if (NULL != time) {
		ret = fits_get_colnum(fptr, CASEINSEN, "TIME", &column, status);
		CHECK_STATUS(ret);
		ret = fits_write_col(fptr, TDOUBLE, column, 1, 1, nbins, time, status);
		CHECK_STATUS(ret);
	}
	if (NULL != phase) {
		ret = fits_get_colnum(fptr, CASEINSEN, "PHASE", &column, status);
		CHECK_STATUS(ret);
		ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, phase, status);
		CHECK_STATUS(ret);
	}
	ret = fits_get_colnum(fptr, CASEINSEN, "FLUX", &column, status);
	CHECK_STATUS(ret);
	ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, flux, status);
	CHECK_STATUS(ret);
	if (NULL != pol_frac) {
		ret = fits_get_colnum(fptr, CASEINSEN, "POL_FRAC", &column, status);
		CHECK_STATUS(ret);
		ret = fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, pol_frac,
				status);
		CHECK_STATUS(ret);
	}
	if (NULL != pol_dir) {
		ret = fits_get_colnum(fptr, CASEINSEN, "POL_DIR", &column, status);
		CHECK_STATUS(ret);
		ret
				= fits_write_col(fptr, TFLOAT, column, 1, 1, nbins, pol_dir,
						status);
		CHECK_STATUS(ret);
	}

	// Store the HDU number of the new extension.
	*hdunum = fits_get_hdu_num(fptr, hdunum);

	// Close the file.
	ret = fits_close_file(fptr, status);
	CHECK_STATUS(ret);

	return (*status);
}

static void simput_ext_id(char* const id, const char* const filename,
		const char* const extname, const int extver) {
	sprintf(id, "%s[%s,%d]", filename, extname, extver);

	/*
	 strcpy(id, filename);
	 strcat(id, "[");
	 strcat(id, extname);
	 strcat(id, ",");
	 strcat(id, sprintf("%d", extver));
	 strcat(id, "]");
	 */
}

static SIMPUT_SrcCtlg* simput_open_srctbl(const char* const filename,
		int* const status) {
	SIMPUT_SrcCtlg* srcctlg = (SIMPUT_SrcCtlg*) malloc(sizeof(SIMPUT_SrcCtlg));
	CHECK_NULL(srcctlg,*status,"memory allocation for data structure failed",srcctlg);

	// Check if the file already exists.
	int exists;
	int ret = fits_file_exists(filename, &exists, status);
	CHECK_STATUS_RET(ret,srcctlg);

	// If no, create a new file.
	if (1 != exists) {
		// Create and open a new empty FITS file.
		ret = fits_create_file(&srcctlg->fptr, filename, status);
		CHECK_STATUS_RET(ret,srcctlg);
	} else {
		// The file does already exist, so just open it.
		ret = fits_open_file(&srcctlg->fptr, filename, READWRITE, status);
		CHECK_STATUS_RET(ret,srcctlg);
	}
	// END of check, whether the file already exists or not.

	// Check if a source catalog extension exists.
	// Try to move the internal HDU pointer of the fitsfile data structure
	// to the right extension containing the source catalog.
	fits_write_errmark();
	int temp_status = EXIT_SUCCESS; // We have to use another status variable here.
	fits_movnam_hdu(srcctlg->fptr, BINARY_TBL, "SRC_CAT", 0, &temp_status);
	fits_clear_errmark();

	if (BAD_HDU_NUM == temp_status) {
		// Create the table structure for the source catalog.
		char *ttype[] = { "SRC_ID", "SRC_NAME", "RA", "DEC", "FLUX", "E_MIN",
				"E_MAX", "SPECTRUM", "LIGHTCUR", "IMAGE" };
		char *tform[] = { "J", "20A", "E", "E", "E", "E", "E", "21A", "21A",
				"21A" };
		char *tunit[] = { "", "", "deg", "deg", "erg/s/cm^2", "keV", "keV", "",
				"", "" };
		ret = fits_create_tbl(srcctlg->fptr, BINARY_TBL, 0, N_SRC_CAT_COLUMNS, ttype,
				tform, tunit, "SRC_CAT", status);
		CHECK_STATUS_RET(ret,srcctlg);

		// Insert header keywords.
		// ...

	}

	// Note: we do not want to apply a full consistency check here
	// because we also want to be able to access a file, which is
	// not complete yet.

	// Determine the column numbers.
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "SRC_ID", &srcctlg->csrc_id, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "SRC_NAME", &srcctlg->csrc_name, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "RA", &srcctlg->cra, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "DEC", &srcctlg->cdec, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "FLUX", &srcctlg->cflux, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "E_MIN", &srcctlg->ce_min, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "E_MAX", &srcctlg->ce_max, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "SPECTRUM", &srcctlg->cspectrum, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "LIGHTCUR", &srcctlg->clightcur, status);
	CHECK_STATUS_RET(ret,srcctlg);
	ret = fits_get_colnum(srcctlg->fptr, CASEINSEN, "IMAGE", &srcctlg->cimage, status);
	CHECK_STATUS_RET(ret,srcctlg);

	return(srcctlg);
}
