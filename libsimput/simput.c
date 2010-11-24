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


#define N_SRC_CAT_COLUMNS (10)

#define CHECKSTATUS(a) if (EXIT_SUCCESS!=a) return(a)


/** Open a SIMPUT FITS file. Open either an existing SIMPUT file or
    create a new and empty one and open it. The file access mode can
    be either READWRITE or READONLY. */
int simput_open_srcctlg(simputfile** fptr,
			const char* const filename,
			const int mode,
			int* const status)
{
  // Initialize pointer with NULL.
  *fptr = NULL;

  // Return value of CFITSIO routines.
  int ret;

  // Check if the requested file already exists.
  int exists;
  ret = fits_file_exists(filename, &exists, status);
  CHECKSTATUS(ret);

  // If it doesn't exist, create a new empty file containing
  // the basic table structure (source catalog table).
  if (1!=exists) { 
    // Maybe have to change the 1 in the comparison above.
    // There might also be other possible values.

    // Create and open a new empty FITS file.
    ret = fits_create_file(fptr, filename, status);
    CHECKSTATUS(ret);

    // Create header keywords.

    // Create the table structure for the source catalog.
    char *ttype[] = 
      { "SRC_ID", "SRC_NAME", "RA", "DEC", "FLUX", "E_MIN", "E_MAX",
	"SPECTRUM", "LIGHTCUR", "IMAGE" };
    char *tform[] = 
      { "J", "21A", "E", "E", "E", "E", "E", "21A", "21A", "21A" };
    char *tunit[] = 
      { "", "", "deg", "deg", "erg/s/cm^2", "keV", "keV", "", "", "" };
    ret = fits_create_tbl(*fptr, BINARY_TBL, 0, N_SRC_CAT_COLUMNS, 
			  ttype, tform, tunit, "SRC_CAT", status);
    CHECKSTATUS(ret);

  } else {
    // The file does already exist.

    // Try to open the file via CFITSIO with the requested mode.
    ret = fits_open_file(fptr, filename, mode, status);
    CHECKSTATUS(ret);
    
    // Note: we do not want to apply a full consistency check here
    // because we also want to be able to access a file, which is
    // not complete yet.
  }
  // END of check if file exists or not.

  // TODO Collect some basic information about the file like
  // the HDU number of the source catalog and store it 
  // in an appropriate data structure.

  return(EXIT_SUCCESS);
}


/** Close an open SIMPUT FITS file. */
int simput_close_srcctlg(simputfile** fptr,
			 int* const status) 
{
  int ret = fits_close_file(*fptr, status);
  CHECKSTATUS(ret);

  *fptr = NULL;

  return(EXIT_SUCCESS);
}


int simput_add_src(simputfile* fptr) 
{

}

