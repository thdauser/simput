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

#include "common.h"




// Pointer to the spectrum cache
// This is initialized when the first spectrum is read
static SimputSpecExtCache *SpecCache = NULL;
static char **NamePtr = NULL;

static void read_unit(fitsfile* const fptr, const int column,
		      char* unit, int* const status)
{
  // Read the header keyword.
  char keyword[SIMPUT_MAXSTR], comment[SIMPUT_MAXSTR];
  sprintf(keyword, "TUNIT%d", column);
  fits_read_key(fptr, TSTRING, keyword, unit, comment, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not read FITS keyword '%s'", keyword);
    SIMPUT_ERROR(msg);
    return;
  }
}


static float unit_conversion_rad(const char* const unit)
{
  if (0==strcmp(unit, "rad")) {
    return(1.);
  } else if (0==strcmp(unit, "deg")) {
    return(M_PI/180.);
  } else if (0==strcmp(unit, "arcmin")) {
    return(M_PI/180./60.);
  } else if (0==strcmp(unit, "arcsec")) {
    return(M_PI/180./3600.);
  } else {
    // Unknown units.
    return(0.);
  }
}


static float unit_conversion_keV(const char* const unit)
{
  if (0==strcmp(unit, "keV")) {
    return(1.);
  } else if (0==strcmp(unit, "eV")) {
    return(0.001);
  } else {
    // Unknown units.
    return(0.);
  }
}


static float unit_conversion_ergpspcm2(const char* const unit)
{
  if (0==strcmp(unit, "erg/s/cm**2")) {
    return(1.);
  } else {
    // Unknown units.
    return(0.);
  }
}


static float unit_conversion_phpspcm2pkeV(const char* const unit)
{
  if (0==strcmp(unit, "photon/s/cm**2/keV")) {
    return(1.);
  } else if (0==strcmp(unit, "photons/s/cm**2/keV")) {
    return(1.);
  } else {
    // Unknown units.
    return(0.);
  }
}


static float unit_conversion_s(const char* const unit)
{
  if ((0==strcmp(unit, "s")) ||
      (0==strcmp(unit, "Hz^-1")) || (0==strcmp(unit, "Hz**-1")) ||
      (0==strcmp(unit, "1/Hz"))){
    return(1.);
  } else if (0==strcmp(unit, "min")) {
    return(60.);
  } else if (0==strcmp(unit, "h")) {
    return(3600.);
  } else if (0==strcmp(unit, "d")) {
    return(24.*3600.);
  } else if (0==strcmp(unit, "yr")) {
    return(365.25*24.*3600.);
  } else {
    // Unknown units.
    return(0.);
  }
}


static float unit_conversion_Hz(const char* const unit)
{
  if ((0==strcmp(unit, "Hz")) ||
      (0==strcmp(unit, "s^-1")) || (0==strcmp(unit, "s**-1")) ){
    return(1.);
  } else {
    // Unknown units.
    return(0.);
  }
}


SimputCtlg* openSimputCtlg(const char* const filename,
			   const int mode,
			   const int maxstrlen_src_name,
			   const int maxstrlen_spectrum,
			   const int maxstrlen_image,
			   const int maxstrlen_timing,
			   int* const status)
{
  SimputCtlg* cat=newSimputCtlg(status);
  CHECK_STATUS_RET(*status, cat);

  char **tform=NULL;

  do { // Error handling loop.

    // Store the filename and filepath of the FITS file containing
    // the source catalog.
    char cfilename[SIMPUT_MAXSTR];
    char rootname[SIMPUT_MAXSTR];
    // Make a local copy of the filename variable in order to avoid
    // compiler warnings due to discarded const qualifier at the
    // subsequent function call.
    strcpy(cfilename, filename);
    fits_parse_rootname(cfilename, rootname, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not determine rootname of '%s'", cfilename);
      SIMPUT_ERROR(msg);
      break;
    }
    CHECK_STATUS_BREAK(*status);

    // Split rootname into the file path and the file name.
    char* lastslash=strrchr(rootname, '/');
    if (NULL==lastslash) {
      cat->filepath=(char*)malloc(sizeof(char));
      CHECK_NULL_BREAK(cat->filepath, *status,
		       "memory allocation for filepath failed");
      cat->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
      CHECK_NULL_BREAK(cat->filename, *status,
		       "memory allocation for filename failed");
      cat->filepath[0]='\0';
      strcpy(cat->filename, rootname);
    } else {
      lastslash++;
      cat->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
      CHECK_NULL_BREAK(cat->filename, *status,
		       "memory allocation for filename failed");
      strcpy(cat->filename, lastslash);

      *lastslash='\0';
      cat->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
      CHECK_NULL_BREAK(cat->filepath, *status,
		       "memory allocation for filepath failed");
      strcpy(cat->filepath, rootname);
    }
    // END of storing the filename and filepath.

    // Check if the file already exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "check failed whether file '%s' exists", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (1==exists) {
      // The file already exists => open it.
      fits_open_file(&cat->fptr, filename, mode, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not open file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

    } else if (READWRITE==mode) {
      // The file does not exist, but it shall be created.
      fits_create_file(&cat->fptr, filename, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

    } else {
      // The file should be opened for read access, but it does not exist.
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "file '%s' does not exist", filename);
      SIMPUT_ERROR(msg);
      *status=EXIT_FAILURE;
      break;
    }

    // Try to move to the source catalog extension.
    int status2=EXIT_SUCCESS;
    fits_write_errmark();
    fits_movnam_hdu(cat->fptr, BINARY_TBL, "SRC_CAT", 0, &status2);
    fits_clear_errmark();
    if (BAD_HDU_NUM==status2) {
      if (READWRITE==mode) {
	// The file does not contain a source catalog => create one.
	char *ttype[]={"SRC_ID", "SRC_NAME", "RA", "DEC", "IMGROTA", "IMGSCAL",
		       "E_MIN", "E_MAX", "FLUX", "SPECTRUM", "IMAGE", "TIMING"};
	char *tunit[]={"", "", "deg", "deg", "deg", "",
		       "keV", "keV", "erg/s/cm**2", "", "", ""};
	tform=(char**)malloc(12*sizeof(char*));
	CHECK_NULL_BREAK(tform, *status,
			 "memory allocation for table parameters failed");
	int ii;
	for (ii=0; ii<12; ii++) {
	  tform[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	  CHECK_NULL_BREAK(tform[ii], *status,
			   "memory allocation for table parameters failed");
	}
	CHECK_STATUS_BREAK(*status);
	strcpy(tform[0], "J");
	if (0==maxstrlen_src_name) {
	  strcpy(tform[1], "1PA");
	} else {
	  sprintf(tform[1], "%dA", maxstrlen_src_name);
	}
	strcpy(tform[2], "D");
	strcpy(tform[3], "D");
	strcpy(tform[4], "E");
	strcpy(tform[5], "E");
	strcpy(tform[6], "E");
	strcpy(tform[7], "E");
	strcpy(tform[8], "E");
	if (0==maxstrlen_spectrum) {
	  strcpy(tform[9], "1PA");
	} else {
	  sprintf(tform[9], "%dA", maxstrlen_spectrum);
	}
	if (0==maxstrlen_image) {
	  strcpy(tform[10], "1PA");
	} else {
	  sprintf(tform[10], "%dA", maxstrlen_image);
	}
	if (0==maxstrlen_timing) {
	  strcpy(tform[11], "1PA");
	} else {
	  sprintf(tform[11], "%dA", maxstrlen_timing);
	}
	fits_create_tbl(cat->fptr, BINARY_TBL, 0, 12, ttype, tform, tunit,
			"SRC_CAT", status);
	if (EXIT_SUCCESS!=*status) {
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "could not create binary table for source catalog "
		  "in file '%s'", filename);
	  SIMPUT_ERROR(msg);
	  break;
	}

	// Write the necessary header keywords.
	fits_write_key(cat->fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
	fits_write_key(cat->fptr, TSTRING, "HDUCLAS1", "SRC_CAT", "", status);
	fits_write_key(cat->fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
	fits_write_key(cat->fptr, TSTRING, "RADESYS", "FK5", "", status);
	float equinox=2000.0;
	fits_write_key(cat->fptr, TFLOAT, "EQUINOX", &equinox, "", status);
	if (EXIT_SUCCESS!=*status) {
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "failed writing FITS keywords in file '%s'", filename);
	  SIMPUT_ERROR(msg);
	  break;
	}

      } else {
	// The file is opened for read access, but does not contain
	// a source catalog extension.
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "file '%s' does not contain a source catalog", filename);
	SIMPUT_ERROR(msg);
	*status=EXIT_FAILURE;
	break;
      }
    }

    // Get the column names.
    // Required columns:
    fits_get_colnum(cat->fptr, CASEINSEN, "SRC_ID", &cat->csrc_id, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'SRC_ID' in source catalog");
      break;
    }

    fits_get_colnum(cat->fptr, CASEINSEN, "RA", &cat->cra, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'RA' in source catalog");
      break;
    }

    fits_get_colnum(cat->fptr, CASEINSEN, "DEC", &cat->cdec, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'DEC' in source catalog");
      break;
    }

    fits_get_colnum(cat->fptr, CASEINSEN, "E_MIN", &cat->ce_min, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'E_MIN' in source catalog");
      break;
    }

    fits_get_colnum(cat->fptr, CASEINSEN, "E_MAX", &cat->ce_max, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'E_MAX' in source catalog");
      break;
    }

    fits_get_colnum(cat->fptr, CASEINSEN, "FLUX", &cat->cflux, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'FLUX' in source catalog");
      break;
    }

    fits_get_colnum(cat->fptr, CASEINSEN, "SPECTRUM", &cat->cspectrum, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not find column 'SPECTRUM' in source catalog");
      break;
    }

    // Optional columns:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(cat->fptr, CASEINSEN, "SRC_NAME", &cat->csrc_name, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(cat->fptr, CASEINSEN, "IMGROTA", &cat->cimgrota, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(cat->fptr, CASEINSEN, "IMGSCAL", &cat->cimgscal, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(cat->fptr, CASEINSEN, "IMAGE", &cat->cimage, &opt_status);
    opt_status=EXIT_SUCCESS;

    fits_get_colnum(cat->fptr, CASEINSEN, "TIMING", &cat->ctiming, &opt_status);
    // For compatibility with SIMPUT version 1.0.0.
    if (EXIT_SUCCESS!=opt_status) {
      opt_status=EXIT_SUCCESS;
      fits_get_colnum(cat->fptr, CASEINSEN, "LIGHTCUR", &cat->ctiming, &opt_status);
    }
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Take care of the units. Determine conversion factors.
    char ura[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->cra, ura, status);
    CHECK_STATUS_BREAK(*status);
    cat->fra = unit_conversion_rad(ura);
    if (0.==cat->fra) {
      SIMPUT_ERROR("unknown units of 'RA' column");
      *status=EXIT_FAILURE;
      break;
    }

    char udec[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->cdec, udec, status);
    CHECK_STATUS_BREAK(*status);
    cat->fdec = unit_conversion_rad(udec);
    if (0.==cat->fdec) {
      SIMPUT_ERROR("unknown units of 'DEC' column");
      *status=EXIT_FAILURE;
      break;
    }

    if (cat->cimgrota>0) {
      char uimgrota[SIMPUT_MAXSTR];
      read_unit(cat->fptr, cat->cimgrota, uimgrota, status);
      CHECK_STATUS_BREAK(*status);
      cat->fimgrota = unit_conversion_rad(uimgrota);
      if (0.==cat->fimgrota) {
	SIMPUT_ERROR("unknown units of 'IMGROTA' column");
	*status=EXIT_FAILURE;
	break;
      }
    }

    char ue_min[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->ce_min, ue_min, status);
    CHECK_STATUS_BREAK(*status);
    cat->fe_min = unit_conversion_keV(ue_min);
    if (0.==cat->fe_min) {
      SIMPUT_ERROR("unknown units of 'E_MIN' column");
      *status=EXIT_FAILURE;
      break;
    }

    char ue_max[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->ce_max, ue_max, status);
    CHECK_STATUS_BREAK(*status);
    cat->fe_max = unit_conversion_keV(ue_max);
    if (0.==cat->fe_max) {
      SIMPUT_ERROR("unknown units of 'E_MAX' column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->cflux, uflux, status);
    CHECK_STATUS_BREAK(*status);
    cat->fflux = unit_conversion_ergpspcm2(uflux);
    if (0.==cat->fflux) {
      SIMPUT_ERROR("unknown units of 'FLUX' column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries.
    fits_get_num_rows(cat->fptr, &cat->nentries, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine number of entries in source catalog");
      break;
    }

  } while(0); // END of error handling loop.

  // Release memory.
  if (NULL!=tform) {
    int ii;
    for (ii=0; ii<12; ii++) {
      if (NULL!=tform[ii]) {
	free(tform[ii]);
      }
    }
  }

  return(cat);
}


SimputSrc* loadSimputSrc(SimputCtlg* const cat,
			 const long row,
			 int* const status)
{
  SimputSrc* src=NULL;

  // String buffers.
  char* src_name[1]={NULL};
  char* spectrum[1]={NULL};
  char* image[1]   ={NULL};
  char* timing[1]  ={NULL};

  // Check if the specified row number is within the valid range.
  if ((row<=0) || (row>cat->nentries)) {
    SIMPUT_ERROR("invalid row number");
    return(NULL);
  }

  do { // Beginning of error handling loop.

    // Allocate memory for string buffers.
    src_name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(src_name[0], *status,
		     "memory allocation for string buffer failed");
    spectrum[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(spectrum[0], *status,
		     "memory allocation for string buffer failed");
    image[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(image[0], *status,
		     "memory allocation for string buffer failed");
    timing[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(timing[0], *status,
		     "memory allocation for string buffer failed");

    long src_id=0;
    double ra=0., dec=0.;
    float imgrota=0., imgscal=1.;
    float e_min=0., e_max=0., flux=0.;

    // Read the data from the table.
    int anynul=0;
    fits_read_col(cat->fptr, TLONG, cat->csrc_id, row, 1, 1,
		  &src_id, &src_id, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading source ID from source catalog");
      break;
    }

    if (cat->csrc_name>0) {
      fits_read_col_str(cat->fptr, cat->csrc_name, row, 1, 1, "",
			src_name, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading source name from source catalog");
	break;
      }
    } else {
      src_name[0][0]='\0';
    }

    fits_read_col(cat->fptr, TDOUBLE, cat->cra, row, 1, 1, &ra, &ra, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading right ascension from source catalog");
      break;
    }
    ra *=cat->fra; // Convert to [rad].

    fits_read_col(cat->fptr, TDOUBLE, cat->cdec, row, 1, 1, &dec, &dec, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading declination from source catalog");
      break;
    }
    dec*=cat->fdec; // Convert to [rad].

    if (cat->cimgrota>0) {
      fits_read_col(cat->fptr, TFLOAT, cat->cimgrota, row, 1, 1,
		    &imgrota, &imgrota, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading image rotation angle from source catalog");
	break;
      }
      imgrota*=cat->fimgrota; // Convert to [rad].
    }
    if (cat->cimgscal>0) {
      fits_read_col(cat->fptr, TFLOAT, cat->cimgscal, row, 1, 1,
		    &imgscal, &imgscal, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
        SIMPUT_ERROR("failed reading image scaling from source catalog");
        break;
      }

      // Check if imgscal value is valid
      if (imgscal == 0) {
        char buffer[512];
        sprintf(buffer, "IMGSCAL of source %ld in %s must be non-zero",
          src_id, cat->filename);
        SIMPUT_ERROR(buffer);
        *status=EXIT_FAILURE;
        break;
      }
      if (imgscal < 1.e-8) {
        char buffer[512];
        sprintf(buffer, "IMGSCAL value of source %ld in %s is very small (IMGSCAL = %e)",
          src_id, cat->filename, imgscal);
        SIMPUT_WARNING(buffer);
      }
    }

    fits_read_col(cat->fptr, TFLOAT, cat->ce_min, row, 1, 1,
		  &e_min, &e_min, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading 'E_MIN' from source catalog");
      break;
    }
    e_min*=cat->fe_min; // Convert to [keV].

    fits_read_col(cat->fptr, TFLOAT, cat->ce_max, row, 1, 1,
		  &e_max, &e_max, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading 'E_MAX' from source catalog");
      break;
    }
    e_max*=cat->fe_max; // Convert to [keV].

    fits_read_col(cat->fptr, TFLOAT, cat->cflux, row, 1, 1,
		  &flux, &flux, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading reference flux from source catalog");
      break;
    }
    flux*=cat->fflux; // Convert to [erg/s/cm**2].
    if (flux > 1.e20) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "flux (%e erg/s/cm**2) exceeds maximum value, "
	      "therefore reset to 0", flux);
      SIMPUT_WARNING(msg);
      flux=0.;
    }

    fits_read_col(cat->fptr, TSTRING, cat->cspectrum, row, 1, 1,
		  "", spectrum, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading spectrum reference from source catalog");
      break;
    }

    if (cat->cimage>0) {
      fits_read_col(cat->fptr, TSTRING, cat->cimage, row, 1, 1,
		    "", image, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading image reference from source catalog");
	break;
      }
    } else {
      image[0][0]='\0';
    }

    if (cat->ctiming>0) {
      fits_read_col(cat->fptr, TSTRING, cat->ctiming, row, 1, 1,
		    "", timing, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading timing reference from source catalog");
	break;
      }
    } else {
      timing[0][0]='\0';
    }

    // Create a new SimputSource data structure.
    src=newSimputSrcV(src_id, src_name[0], ra, dec, imgrota, imgscal,
		      e_min, e_max, flux, spectrum[0], image[0],
		      timing[0], status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  // Release memory.
  if (NULL!=src_name[0]) free(src_name[0]);
  if (NULL!=spectrum[0]) free(spectrum[0]);
  if (NULL!=image[0])    free(image[0]);
  if (NULL!=timing[0])   free(timing[0]);

  return(src);
}


void appendSimputSrc(SimputCtlg* const cat,
		     SimputSrc* const src,
		     int* const status)
{
  // Write the data.
  fits_insert_rows(cat->fptr, cat->nentries, 1, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed appending new row to source catalog");
    return;
  }
  cat->nentries++;

  fits_write_col(cat->fptr, TLONG, cat->csrc_id, cat->nentries, 1, 1,
		 &src->src_id, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing source ID to catalog");
    return;
  }

  fits_write_col(cat->fptr, TSTRING, cat->csrc_name, cat->nentries, 1, 1,
		 &src->src_name, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing source name to catalog");
    return;
  }

  double ra = src->ra*180./M_PI;
  fits_write_col(cat->fptr, TDOUBLE, cat->cra, cat->nentries, 1, 1,
		 &ra, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing right ascension to source catalog");
    return;
  }

  double dec = src->dec*180./M_PI;
  fits_write_col(cat->fptr, TDOUBLE, cat->cdec, cat->nentries, 1, 1,
		 &dec, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing declination to source catalog");
    return;
  }

  float imgrota = src->imgrota*180./M_PI;
  fits_write_col(cat->fptr, TFLOAT, cat->cimgrota, cat->nentries, 1, 1,
		 &imgrota, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing image rotation angle to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TFLOAT, cat->cimgscal, cat->nentries, 1, 1,
		 &src->imgscal, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing image scaling to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TFLOAT, cat->ce_min, cat->nentries, 1, 1,
		 &src->e_min, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing 'E_MIN' to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TFLOAT, cat->ce_max, cat->nentries, 1, 1,
		 &src->e_max, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing 'E_MAX' to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TFLOAT, cat->cflux, cat->nentries, 1, 1,
		 &src->eflux, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing reference flux to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TSTRING, cat->cspectrum, cat->nentries, 1, 1,
		 &src->spectrum, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing spectrum reference to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TSTRING, cat->cimage, cat->nentries, 1, 1,
		 &src->image, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing image reference to source catalog");
    return;
  }

  fits_write_col(cat->fptr, TSTRING, cat->ctiming, cat->nentries, 1, 1,
		 &src->timing, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed writing timing reference to source catalog");
    return;
  }
}


void appendSimputSrcBlock(SimputCtlg* const cat,
			  SimputSrc** const src,
			  const long nsources,
			  int* const status)
{
  // Insert new rows.
  fits_insert_rows(cat->fptr, cat->nentries, nsources, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("failed appending new rows to source catalog");
    return;
  }

  long first=cat->nentries+1;
  cat->nentries+=nsources;

  // Buffers.
  long* src_id=NULL;
  double* ra=NULL;
  double* dec=NULL;
  float* imgrota=NULL;
  float* imgscal=NULL;
  float* e_min=NULL;
  float* e_max=NULL;
  float* eflux=NULL;

  // Beginning of error handling loop.
  do {

    // Allocate memory for buffers.
    src_id=(long*)malloc(nsources*sizeof(long));
    CHECK_NULL_BREAK(src_id, *status, "memory allocation failed");
    ra =(double*)malloc(nsources*sizeof(double));
    CHECK_NULL_BREAK(ra, *status, "memory allocation failed");
    dec=(double*)malloc(nsources*sizeof(double));
    CHECK_NULL_BREAK(dec, *status, "memory allocation failed");
    imgrota=(float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(imgrota, *status, "memory allocation failed");
    imgscal=(float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(imgscal, *status, "memory allocation failed");
    e_min  =(float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(e_min, *status, "memory allocation failed");
    e_max  =(float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(e_max, *status, "memory allocation failed");
    eflux  =(float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(eflux, *status, "memory allocation failed");


    // Write the data.
    long ii;
    for (ii=0; ii<nsources; ii++) {

      // Copy values to buffers.
      src_id[ii] =src[ii]->src_id;
      ra[ii]     =src[ii]->ra *180./M_PI;
      dec[ii]    =src[ii]->dec*180./M_PI;
      imgrota[ii]=src[ii]->imgrota;
      imgscal[ii]=src[ii]->imgscal;
      e_min[ii]  =src[ii]->e_min;
      e_max[ii]  =src[ii]->e_max;
      eflux[ii]  =src[ii]->eflux;

      // Write strings.
      fits_write_col(cat->fptr, TSTRING, cat->csrc_name, first+ii, 1, 1,
		     &(src[ii]->src_name), status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing source name to catalog");
	break;
      }

      fits_write_col(cat->fptr, TSTRING, cat->cspectrum, first+ii, 1, 1,
		     &(src[ii]->spectrum), status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing spectrum reference to source catalog");
	break;
      }

      fits_write_col(cat->fptr, TSTRING, cat->cimage, first+ii, 1, 1,
		     &(src[ii]->image), status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing image reference to source catalog");
	break;
      }

      fits_write_col(cat->fptr, TSTRING, cat->ctiming, first+ii, 1, 1,
		     &(src[ii]->timing), status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing timing reference to source catalog");
	break;
      }
    }
    CHECK_STATUS_BREAK(*status);
    // END of loop over all sources.

    // Write the buffers to the file.
    fits_write_col(cat->fptr, TLONG, cat->csrc_id, first, 1, nsources,
		   src_id, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing source IDs to catalog");
      break;
    }

    fits_write_col(cat->fptr, TDOUBLE, cat->cra, first, 1, nsources,
		   ra, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing right ascensions to catalog");
      break;
    }

    fits_write_col(cat->fptr, TDOUBLE, cat->cdec, first, 1, nsources,
		   dec, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing declinations to catalog");
      break;
    }

    fits_write_col(cat->fptr, TFLOAT, cat->cimgrota, first, 1, nsources,
		   imgrota, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing image rotation angles to catalog");
      break;
    }

    fits_write_col(cat->fptr, TFLOAT, cat->cimgscal, first, 1, nsources,
		   imgscal, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing image scalings to catalog");
      break;
    }

    fits_write_col(cat->fptr, TFLOAT, cat->ce_min, first, 1, nsources,
		   e_min, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing 'E_MIN' values to catalog");
      break;
    }

    fits_write_col(cat->fptr, TFLOAT, cat->ce_max, first, 1, nsources,
		   e_max, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing 'E_MAX' values to catalog");
      break;
    }

    fits_write_col(cat->fptr, TFLOAT, cat->cflux, first, 1, nsources,
		   eflux, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing reference fluxes to catalog");
      break;
    }

  } while(0); // END of error handling loop.

  // Free memory.
  if (NULL!=src_id) free(src_id);
  if (NULL!=ra)     free(ra);
  if (NULL!=dec)    free(dec);
  if (NULL!=imgrota) free(imgrota);
  if (NULL!=imgscal) free(imgscal);
  if (NULL!=e_min)   free(e_min);
  if (NULL!=e_max)   free(e_max);
  if (NULL!=eflux)   free(eflux);
}


// WARNING: THIS FUNCTION IS NOT WORKING YET AS IT SHOULD BE !!!
uniqueSimputident* get_simput_ident(char* filename, int type, int *status){

	uniqueSimputident *ident = (uniqueSimputident*) malloc (sizeof(uniqueSimputident));
	CHECK_NULL_RET(ident,*status,"error: memory allocation failed",NULL);

	// Open the specified FITS file. The filename must uniquely identify
	// the spectrum contained in a binary table via the extended filename
	// syntax. It must even specify the row, in which the spectrum is
	// contained. Therefore we do not have to care about the HDU or row
	// number.
	fitsfile* fptr=NULL;

	if (type == SIMPUT_IMG_TYPE){
		fits_open_image(&fptr, filename, READONLY, status);
	} else {
		fits_open_file(&fptr, filename, READONLY, status);
	}
	if (EXIT_SUCCESS!=*status) {
		char msg[SIMPUT_MAXSTR];
		sprintf(msg, "could not open FITS table in file '%s'", filename);
		SIMPUT_ERROR(msg);
	}

	CHECK_NULL_RET(fptr,*status,"error: empty pointer to fitsfile",NULL);
	CHECK_NULL_RET(fptr->Fptr,*status,"error: empty pointer to fits structure",NULL);

	ident->filename = strdup((fptr->Fptr)->filename);
	ident->io_pos = (fptr->Fptr)->io_pos+(fptr->Fptr)->datastart;

	if (NULL!=fptr) fits_close_file(fptr, status);
	CHECK_STATUS_RET(*status, NULL);

	return ident;
}


SimputMIdpSpec* loadSimputMIdpSpec(const char* const filename,
				   int* const status)
{
  // String buffer.
  char* name[1]={NULL};

  SimputMIdpSpec* spec;

  if ( SpecCache == NULL )
  {
    headas_chat(5, "Initializing spectrum cache\n");
    initSpecCache();
  }

  if ( SpecCache )
  {
    char *basename, *extname;
    int extver;
    char msg[SIMPUT_MAXSTR];
    char *expr;
    long ind;
    long row;

    do {

      expr = scanSpecFileName((char *) filename, &basename, &extname, &extver, status);
      if ( *status )
      {
        sprintf(msg, "Error scanning filename, can not use caching, falling back to default");
        SIMPUT_WARNING(msg);
	*status = EXIT_SUCCESS;
	break ;
      }

      headas_chat(5, "Checking whether %s is already cached\n", filename);
      ind = specIsCached(basename, extname, extver);

      if ( ind == -1 )
      {
	headas_chat(5, "No, have to open, looking for where to open\n");
	ind = getNextSpecCache();
	headas_chat(5, "Opening to index %ld\n", ind);
	if ( SpecCache->filename[ind] != NULL )
	{
	  headas_chat(5, "Index %ld was in use, clearing\n", ind);
	  destroyNthSpecCache(SpecCache, ind);
	}
	openNthSpecCache(basename, extname, extver, ind, status);
      } else {
	headas_chat(5, "Yes, not opening again\n");
      }

      row = getSpecRow(expr, ind);

      if (row>=0) {
	 spec = readCacheSpec(ind, row, (char *) filename, status);
      } else {
	 // not found
	 spec=NULL;
	 *status=EXIT_FAILURE;
      }

      if ( spec == NULL || *status != EXIT_SUCCESS )
      {
	headas_chat(5, "*** Warning in loadSimputMIdpSpec: Not possible to use spectrum cache! ***\n");
	*status = EXIT_SUCCESS;
	break;
      }

      free(basename);
      free(extname);

      return spec;
    } while (0) ;
    free(basename);
    free(extname);
  }

  headas_chat(5, "Falling back to default reading method\n");

  spec=newSimputMIdpSpec(status);
  CHECK_STATUS_RET(*status, spec);


  // Open the specified FITS file. The filename must uniquely identify
  // the spectrum contained in a binary table via the extended filename
  // syntax. It must even specify the row, in which the spectrum is
  // contained. Therefore we do not have to care about the HDU or row
  // number.
  fitsfile* fptr=NULL;
  fits_open_table(&fptr, filename, READONLY, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not open FITS table in file '%s'", filename);
    SIMPUT_ERROR(msg);
    return(spec);
  }

  do { // Error handling loop.

    // Get the column names.
    int cenergy=0, cfluxdensity=0, cname=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'ENERGY' in spectrum '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "FLUXDENSITY", &cfluxdensity, status);
    if (EXIT_SUCCESS!=*status) {
      // For compatibility with SIMPUT version 1.0.0.
      *status=EXIT_SUCCESS;
      fits_get_colnum(fptr, CASEINSEN, "FLUX", &cfluxdensity, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not find column 'FLUXDENSITY' in spectrum '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    fits_clear_errmark();

    // Optional columnes:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "NAME", &cname, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Determine the unit conversion factors.
    char uenergy[SIMPUT_MAXSTR];
    read_unit(fptr, cenergy, uenergy, status);
    CHECK_STATUS_BREAK(*status);
    float fenergy=unit_conversion_keV(uenergy);
    if (0.==fenergy) {
      SIMPUT_ERROR("unknown units of 'ENERGY' column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(fptr, cfluxdensity, uflux, status);
    CHECK_STATUS_BREAK(*status);
    float fflux=unit_conversion_phpspcm2pkeV(uflux);
    if (0.==fflux) {
      SIMPUT_ERROR("unknown units of 'FLUXDENSITY' column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries in the 2 vector columns.
    int typecode;
    long nenergy, nfluxdensity, width;
    fits_get_coltype(fptr, cenergy, &typecode, &nenergy, &width, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine type of column 'ENERGY'");
      break;
    }

    fits_get_coltype(fptr, cfluxdensity, &typecode, &nfluxdensity, &width, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine type of column 'FLUXDENSITY'");
      break;
    }

    // If the columns are of variable-length data type, the returned repeat
    // value is 1. In that case we have to use another routine to get the
    // number of elements in a particular row.
    if (1==nenergy) {
      long offset;
      fits_read_descript(fptr, cenergy, 1, &nenergy, &offset, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("could not determine type of column 'ENERGY'");
	break;
      }

      fits_read_descript(fptr, cfluxdensity, 1, &nfluxdensity, &offset, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("could not determine type of column 'FLUXDENSITY'");
	break;
      }
    }

    // The number of energy bins and of flux entries must be identical.
    if (nenergy!=nfluxdensity) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "number of energy and flux entries in spectrum '%s' is "
	      "not equivalent", filename);
      SIMPUT_ERROR(msg);
      *status=EXIT_FAILURE;
      break;
    }
    spec->nentries=nenergy;

    // Informational output.
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "spectrum '%s' contains %ld data points",
	    filename, spec->nentries);
    SIMPUT_INFO(msg);

    // Allocate memory for the arrays.
    spec->energy=(float*)malloc(spec->nentries*sizeof(float));
    CHECK_NULL_BREAK(spec->energy, *status,
		     "memory allocation for spectrum failed");
    spec->fluxdensity=(float*)malloc(spec->nentries*sizeof(float));
    CHECK_NULL_BREAK(spec->fluxdensity, *status,
		     "memory allocation for spectrum failed");

    // Allocate memory for string buffer.
    name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(name[0], *status,
		     "memory allocation for string buffer failed");

    // Read the data from the table.
    int anynul=0;
    fits_read_col(fptr, TFLOAT, cenergy, 1, 1, spec->nentries,
		  NULL, spec->energy, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading energy values from spectrum");
      break;
    }

    fits_read_col(fptr, TFLOAT, cfluxdensity, 1, 1, spec->nentries,
		  NULL, spec->fluxdensity, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading flux values from spectrum");
      break;
    }

    if (cname>0) {
    	fits_read_col(fptr, TSTRING, cname, 1, 1, 1, "", name, &anynul, status);
    	if (EXIT_SUCCESS!=*status) {
    		SIMPUT_ERROR("failed reading designator of spectrum");
    		break;
    	}
    } else {
    	name[0][0]='\0';
    }

    // Multiply with unit scaling factor.
    long ii;
    for (ii=0; ii<spec->nentries; ii++) {
      spec->energy[ii]*=fenergy;
      spec->fluxdensity[ii]*=fflux;
    }

    // Copy the name (ID) of the spectrum from the string buffer
    // to the data structure.
    spec->name = (char*)malloc((strlen(name[0]))*sizeof(char)+1);
    CHECK_NULL_BREAK(spec->name, *status,
		     "memory allocation for name string failed");
    strcpy(spec->name, name[0]);

    // Store the file reference to the spectrum for later comparisons.
    spec->fileref=
      (char*)malloc((strlen(filename)+strlen(name[0])+1)*sizeof(char));
    CHECK_NULL_BREAK(spec->fileref, *status,
		     "memory allocation for file reference failed");
    sprintf(spec->fileref, "%s", filename);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=name[0]) free(name[0]);

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, spec);

  return(spec);
}


void loadCacheAllSimputMIdpSpec(SimputCtlg* const cat,
				const char* const filename,
				int* const status)
{
  fitsfile* fptr=NULL;
  char* name[1]={NULL};
  SimputMIdpSpec** sb=NULL; // Buffer for reading in the spectra.
  long nrows;

  do { // Error handling loop.

    // Check if the filename refers to a binary table extension
    // containing mission-independent spectra.
    int exttype=getSimputExtType(cat, filename, status);
    CHECK_STATUS_BREAK(*status);

    // Only mission-independent spectra can be pre-loaded into
    // the cache.
    if (EXTTYPE_MIDPSPEC!=exttype) {
      break;
    }

    // This routine can only be used if the catalog internal buffer
    // is empty.
    if (NULL!=cat->midpspecbuff) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "for pre-loading all spectra from '%s' "
	      "the buffer must be empty", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Open the specified FITS file. The filename must uniquely identify
    // the extension containing the spectra via the extended filename
    // syntax.
    fits_open_table(&fptr, filename, READONLY, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not open FITS table in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Determine the number of rows in the table.
    fits_get_num_rows(fptr, &nrows, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine number of spectra");
      break;
    }

    // Allocate memory for buffering the spectra
    sb=(SimputMIdpSpec**)malloc(nrows*sizeof(SimputMIdpSpec*));
    CHECK_NULL_BREAK(sb, *status,
		     "memory allocation for buffer of spectra failed");

    // Determine the column numbers.
    int cenergy=0, cfluxdensity=0, cname=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'ENERGY' in spectrum '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "FLUXDENSITY", &cfluxdensity, status);
    if (EXIT_SUCCESS!=*status) {
      // For compatibility with SIMPUT version 1.0.0.
      *status=EXIT_SUCCESS;
      fits_get_colnum(fptr, CASEINSEN, "FLUX", &cfluxdensity, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not find column 'FLUXDENSITY' in spectrum '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    fits_clear_errmark();

    // Optional columns:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "NAME", &cname, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Determine the unit conversion factors.
    char uenergy[SIMPUT_MAXSTR];
    read_unit(fptr, cenergy, uenergy, status);
    CHECK_STATUS_BREAK(*status);
    float fenergy=unit_conversion_keV(uenergy);
    if (0.==fenergy) {
      SIMPUT_ERROR("unknown units of 'ENERGY' column");
      *status=EXIT_FAILURE;
      break;
    }

    char ufluxdensity[SIMPUT_MAXSTR];
    read_unit(fptr, cfluxdensity, ufluxdensity, status);
    CHECK_STATUS_BREAK(*status);
    float ffluxdensity=unit_conversion_phpspcm2pkeV(ufluxdensity);
    if (0.==ffluxdensity) {
      SIMPUT_ERROR("unknown units of 'FLUXDENSITY' column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries in the 2 vector columns.
    int typecode;
    long nenergy, nfluxdensity, width;
    fits_get_coltype(fptr, cenergy, &typecode, &nenergy, &width, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine type of column 'ENERGY'");
      break;
    }

    fits_get_coltype(fptr, cfluxdensity, &typecode, &nfluxdensity, &width, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine type of column 'FLUXDENSITY'");
      break;
    }

    // If the columns are of variable-length data type, the returned repeat
    // value is 1. In that case we have to use another routine to get the
    // number of elements in a particular row.
    if ((1==nenergy)&&(1==nfluxdensity)) {
      long offset;
      fits_read_descript(fptr, cenergy, 1, &nenergy, &offset, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("could not determine type of column 'ENERGY'");
	break;
      }

      fits_read_descript(fptr, cfluxdensity, 1, &nfluxdensity, &offset, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("could not determine type of column 'FLUXDENSITY'");
	break;
      }
    }

    // The number of energy bins and of flux entries must be identical.
    if (nenergy!=nfluxdensity) {
      SIMPUT_ERROR("number of energy and flux density entries in spectrum is "
		   "not equivalent");
      *status=EXIT_FAILURE;
      break;
    }

    // Allocate memory for string buffer.
    name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(name[0], *status,
		     "memory allocation for string buffer failed");

    // Load the spectra.
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "load %ld spectra with %ld data points each",
	    nrows, nenergy);
    SIMPUT_INFO(msg);

    long jj;
    for (jj=0; jj<nrows; jj++) {
      // Allocate memory for a new spectrum.
      SimputMIdpSpec* spec=newSimputMIdpSpec(status);
      CHECK_STATUS_BREAK(*status);

      spec->nentries=nenergy;

      // Allocate memory for the arrays.
      spec->energy=(float*)malloc(spec->nentries*sizeof(float));
      CHECK_NULL_BREAK(spec->energy, *status,
		       "memory allocation for spectrum failed");
      spec->fluxdensity=(float*)malloc(spec->nentries*sizeof(float));
      CHECK_NULL_BREAK(spec->fluxdensity, *status,
		       "memory allocation for spectrum failed");

      // Read the data from the table.
      int anynul=0;
      fits_read_col(fptr, TFLOAT, cenergy, jj+1, 1, spec->nentries,
		    NULL, spec->energy, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading energy values from spectrum");
	break;
      }

      fits_read_col(fptr, TFLOAT, cfluxdensity, jj+1, 1, spec->nentries,
		    NULL, spec->fluxdensity, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading flux values from spectrum");
	break;
      }

      if (cname>0) {
	fits_read_col(fptr, TSTRING, cname, jj+1, 1, 1, "", name,
		      &anynul, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading designator of spectrum");
	  break;
	}
      } else {
	name[0][0]='\0';
      }

      // Multiply with unit scaling factor.
      long ii;
      for (ii=0; ii<spec->nentries; ii++) {
	spec->energy[ii]*=fenergy;
	spec->fluxdensity[ii]*=ffluxdensity;
      }

      // Copy the name (ID) of the spectrum from the string buffer
      // to the data structure.
      spec->name=(char*)malloc((strlen(name[0])+1)*sizeof(char));
      CHECK_NULL_BREAK(spec->name, *status,
		      "memory allocation for name string failed");
      strcpy(spec->name, name[0]);

      // Store the file reference to the spectrum for later comparisons.
      spec->fileref=
	(char*)malloc((strlen(filename)+strlen(name[0])+11)*sizeof(char));
      CHECK_NULL_BREAK(spec->fileref, *status,
		       "memory allocation for file reference failed");
      sprintf(spec->fileref, "%s[NAME=='%s']", filename, name[0]);

      // Add the spectrum to the buffer.
      sb[jj]=spec;
    }
    CHECK_STATUS_BREAK(*status);
    // END of reading all spectra.

    // Insert the spectra into the binary tree buffer of the
    // SimputCtlg data structure.
    buildSimputMIdpSpecBuffer(&(cat->midpspecbuff), sb, nrows, 0, status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=name[0]) free(name[0]);
  if (NULL!=sb) free(sb);
  // Note: Do NOT release the memory of the individual spectra contained
  // in the buffer! They are now part of the catalog-internal binary tree
  // buffer.

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


void saveSimputMIdpSpec(SimputMIdpSpec* const spec,
			const char* const filename,
			char* const extname,
			int extver,
			int* const status)
{
  fitsfile* fptr=NULL;

  // String buffer.
  char* name[1]={NULL};

  char **ttype=NULL;
  char **tform=NULL;
  char **tunit=NULL;

  do { // Error handling loop.

    // Check if the EXTNAME has been specified.
    if (NULL==extname) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }
    if (0==strlen(extname)) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }

    // Check if the specified file exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "check failed whether file '%s' exists", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not open file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

    } else {
      // If no, create a new file.
      fits_create_file(&fptr, filename, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    // END of check, whether the specified file exists.


    // Try to move to the specified extension.
    long nrows=0;
    int status2=EXIT_SUCCESS;
    fits_write_errmark();
    fits_movnam_hdu(fptr, BINARY_TBL, extname, extver, &status2);
    fits_clear_errmark();
    if (BAD_HDU_NUM==status2) {
      // If that does not work, create a new binary table.
      // Allocate memory for the format strings.
      ttype=(char**)malloc(3*sizeof(char*));
      tform=(char**)malloc(3*sizeof(char*));
      tunit=(char**)malloc(3*sizeof(char*));
      CHECK_NULL_BREAK(ttype, *status, "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tform, *status, "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tunit, *status, "memory allocation for string buffer failed");
      int ii;
      for (ii=0; ii<3; ii++) {
	ttype[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	tform[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	tunit[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	CHECK_NULL_BREAK(ttype[ii], *status,
			 "memory allocation for string buffer failed");
	CHECK_NULL_BREAK(tform[ii], *status,
			 "memory allocation for string buffer failed");
	CHECK_NULL_BREAK(tunit[ii], *status,
			 "memory allocation for string buffer failed");
      }
      CHECK_STATUS_BREAK(*status);

      // Set up the table format.
      strcpy(ttype[0], "ENERGY");
      sprintf(tform[0], "1PE");
      strcpy(tunit[0], "keV");

      strcpy(ttype[1], "FLUXDENSITY");
      sprintf(tform[1], "1PE");
      strcpy(tunit[1], "photon/s/cm**2/keV");

      strcpy(ttype[2], "NAME");
      strcpy(tform[2], "48A");
      strcpy(tunit[2], "");

      // Create the table.
      fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, extname, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create binary table for spectrum in file '%s'",
		filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // Write header keywords.
      fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS1", "SPECTRUM", "", status);
      fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
      fits_write_key(fptr, TINT, "EXTVER", &extver, "", status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "failed writing FITS keywords in file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // The new table contains now data up to now.
      nrows=0;

    } else {
      // The extension already exists.
      // Determine the number of contained rows.
      fits_get_num_rows(fptr, &nrows, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("could not determine number of entries in spectrum");
	break;
      }
    }
    // END of check, whether the specified extension exists.


    // Determine the column numbers.
    int cenergy=0, cfluxdensity=0, cname=0;
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'ENERGY' in spectrum '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    fits_get_colnum(fptr, CASEINSEN, "FLUXDENSITY", &cfluxdensity, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'FLUXDENSITY' in spectrum '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Optional name column:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "NAME", &cname, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // If data structure contains a name, check if is unique.
    if (NULL!=spec->name) {
      if (strlen(spec->name)>0) {
	// Check if the NAME string is too long.
	if (strlen(spec->name)>48) {
	  SIMPUT_ERROR("'NAME' of spectrum contains more than 48 characters");
	  *status=EXIT_FAILURE;
	  break;
	}

	// Check if the NAME column is present.
	if (0==cname) {
	  SIMPUT_ERROR("spectrum extension does not contain a 'NAME' column");
	  *status=EXIT_FAILURE;
	  break;
	}

	// Allocate memory for string buffer.
	name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	CHECK_NULL_BREAK(name[0], *status,
			 "memory allocation for string buffer failed");
	long row;
	for(row=0; row<nrows; row++) {
	  int anynul=0;
	  fits_read_col(fptr, TSTRING, cname, row+1, 1, 1, "", name, &anynul, status);
	  if (EXIT_SUCCESS!=*status) {
	    SIMPUT_ERROR("failed reading 'NAME' from spectrum");
	    break;
	  }
	  if (0==strcmp(name[0], spec->name)) {
	    SIMPUT_ERROR("'NAME' in spectrum data structure is not unique");
	    *status=EXIT_FAILURE;
	    break;
	  }
	}
	CHECK_STATUS_BREAK(*status);
      }
    }

    // Create a new row in the table and store the data of the spectrum in it.
    fits_insert_rows(fptr, nrows++, 1, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed appending new row to spectrum extension");
      return;
    }

    fits_write_col(fptr, TFLOAT, cenergy, nrows, 1, spec->nentries,
		   spec->energy, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing energy values to spectrum");
      break;
    }

    fits_write_col(fptr, TFLOAT, cfluxdensity, nrows, 1, spec->nentries,
		   spec->fluxdensity, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing flux values to spectrum");
      break;
    }

    if ((cname>0) && (NULL!=spec->name)) {
      fits_write_col(fptr, TSTRING, cname, nrows, 1, 1, &spec->name, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing 'NAME' value to spectrum");
	break;
      }
    }
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=name[0]) free(name[0]);

  if (NULL!=ttype) {
    int ii;
    for (ii=0; ii<3; ii++) {
      if (NULL!=ttype[ii]) free(ttype[ii]);
    }
    free(ttype);
  }
  if (NULL!=tform) {
    int ii;
    for (ii=0; ii<3; ii++) {
      if (NULL!=tform[ii]) free(tform[ii]);
    }
    free(tform);
  }
  if (NULL!=tunit) {
    int ii;
    for (ii=0; ii<3; ii++) {
      if (NULL!=tunit[ii]) free(tunit[ii]);
    }
    free(tunit);
  }

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


void saveSimputMIdpSpecBlock(SimputMIdpSpec** const spec,
			     const long nspec,
			     const char* const filename,
			     char* const extname,
			     int extver,
			     int* const status)
{
  fitsfile* fptr=NULL;

  // String buffer.
  char **ttype=NULL;
  char **tform=NULL;
  char **tunit=NULL;

  do { // Error handling loop.

    // Check if the EXTNAME has been specified.
    if (NULL==extname) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }
    if (0==strlen(extname)) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }

    // Check if the specified file exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "check failed whether file '%s' exists", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not open file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

    } else {
      // If no, create a new file.
      fits_create_file(&fptr, filename, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    // END of check, whether the specified file exists.


    // Try to move to the specified extension.
    long nrows=0;
    int status2=EXIT_SUCCESS;
    fits_write_errmark();
    fits_movnam_hdu(fptr, BINARY_TBL, extname, extver, &status2);
    fits_clear_errmark();
    if (BAD_HDU_NUM==status2) {
      // If that does not work, create a new binary table.
      // Allocate memory for the format strings.
      ttype=(char**)malloc(3*sizeof(char*));
      tform=(char**)malloc(3*sizeof(char*));
      tunit=(char**)malloc(3*sizeof(char*));
      CHECK_NULL_BREAK(ttype, *status, "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tform, *status, "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tunit, *status, "memory allocation for string buffer failed");
      int ii;
      for (ii=0; ii<3; ii++) {
	ttype[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	tform[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	tunit[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
	CHECK_NULL_BREAK(ttype[ii], *status,
			 "memory allocation for string buffer failed");
	CHECK_NULL_BREAK(tform[ii], *status,
			 "memory allocation for string buffer failed");
	CHECK_NULL_BREAK(tunit[ii], *status,
			 "memory allocation for string buffer failed");
      }
      CHECK_STATUS_BREAK(*status);

      // Set up the table format.
      strcpy(ttype[0], "ENERGY");
      sprintf(tform[0], "1PE");
      strcpy(tunit[0], "keV");

      strcpy(ttype[1], "FLUXDENSITY");
      sprintf(tform[1], "1PE");
      strcpy(tunit[1], "photon/s/cm**2/keV");

      strcpy(ttype[2], "NAME");
      strcpy(tform[2], "48A");
      strcpy(tunit[2], "");

      // Create the table.
      fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, extname, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create binary table for spectrum in file '%s'",
		filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // Write header keywords.
      fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS1", "SPECTRUM", "", status);
      fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
      fits_write_key(fptr, TINT, "EXTVER", &extver, "", status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "failed writing FITS keywords in file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // The new table contains now data up to now.
      nrows=0;

    } else {
      // The extension already exists.
      // Determine the number of contained rows.
      fits_get_num_rows(fptr, &nrows, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("could not determine number of entries in spectrum");
	break;
      }
    }
    // END of check, whether the specified extension exists.


    // Determine the column numbers.
    int cenergy=0, cfluxdensity=0, cname=0;
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'ENERGY' in spectrum '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    fits_get_colnum(fptr, CASEINSEN, "FLUXDENSITY", &cfluxdensity, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'FLUXDENSITY' in spectrum '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Optional name column:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "NAME", &cname, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Create new rows in the table and store the data of the spectra in them.
    fits_insert_rows(fptr, nrows, nspec, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed appending new rows to spectrum extension");
      return;
    }

    long jj;
    for (jj=0; jj<nspec; jj++) {
      nrows++;
      fits_write_col(fptr, TFLOAT, cenergy, nrows, 1, spec[jj]->nentries,
		     spec[jj]->energy, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing energy values to spectrum extension");
	break;
      }

      fits_write_col(fptr, TFLOAT, cfluxdensity, nrows, 1, spec[jj]->nentries,
		     spec[jj]->fluxdensity, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing flux density values to spectrum extension");
	break;
      }

      if ((cname>0) && (NULL!=spec[jj]->name)) {
	fits_write_col(fptr, TSTRING, cname, nrows, 1, 1, &spec[jj]->name, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed writing 'NAME' value to spectrum");
	  break;
	}
      }
      CHECK_STATUS_BREAK(*status);
    }
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=ttype) {
    int ii;
    for (ii=0; ii<3; ii++) {
      if (NULL!=ttype[ii]) free(ttype[ii]);
    }
    free(ttype);
  }
  if (NULL!=tform) {
    int ii;
    for (ii=0; ii<3; ii++) {
      if (NULL!=tform[ii]) free(tform[ii]);
    }
    free(tform);
  }
  if (NULL!=tunit) {
    int ii;
    for (ii=0; ii<3; ii++) {
      if (NULL!=tunit[ii]) free(tunit[ii]);
    }
    free(tunit);
  }

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


SimputLC* loadSimputLC(const char* const filename, int* const status)
{
  // String buffers.
  char* spectrum[1]={NULL};
  char* image[1]={NULL};

  SimputLC* lc=NULL;
  fitsfile* fptr=NULL;

  do { // Error handling loop.

    // Open the specified FITS file. The filename must uniquely identify
    // the light curve contained in a binary table via the extended filename
    // syntax. Therefore we do not have to care about the HDU number.
    fits_open_table(&fptr, filename, READONLY, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not open FITS table in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Get an empty SimputLC data structure.
    lc=newSimputLC(status);
    CHECK_STATUS_BREAK(*status);

    // Get the column names.
    int ctime=0, cphase=0, cflux=0, cspectrum=0, cimage=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'FLUX' in light curve '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Optional columnes:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "TIME", &ctime, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(fptr, CASEINSEN, "PHASE", &cphase, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(fptr, CASEINSEN, "SPECTRUM", &cspectrum, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(fptr, CASEINSEN, "IMAGE", &cimage, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Check, whether there is either a TIME or a PHASE column (but not both).
    if ((0==ctime)&&(0==cphase)) {
      SIMPUT_ERROR("table extension contains neither 'TIME' nor 'PHASE' column");
      *status=EXIT_FAILURE;
      return(lc);
    } else if ((ctime>0)&&(cphase>0)) {
      SIMPUT_ERROR("table extension contains both 'TIME' and 'PHASE' column");
      *status=EXIT_FAILURE;
      return(lc);
    }

    // We check that the phase and time columns are of the correct type:
    int typecode=-1;
    long repeat=-1;
    long width=-1;
    // Check the phase:
    if (cphase){
      fits_get_coltype(fptr, cphase, &typecode, &repeat, &width, status);
      if (repeat!=1){
	SIMPUT_ERROR("'PHASE' column should contain scalar values, not arrays");
	*status=EXIT_FAILURE;
	return(lc);
      }
      *status=EXIT_SUCCESS;
    }
    // Check the time:
    else { // We have already check that there is either phase or time column
      fits_get_coltype(fptr, ctime, &typecode, &repeat, &width, status);
      if (repeat!=1){
	SIMPUT_ERROR("'TIME' column should contain scalar values, not arrays");
	*status=EXIT_FAILURE;
	return(lc);
      }
      *status=EXIT_SUCCESS;
    }

    // Determine the unit conversion factors.
    float ftime=0.;
    if (ctime>0) {
      char utime[SIMPUT_MAXSTR];
      read_unit(fptr, ctime, utime, status);
      CHECK_STATUS_BREAK(*status);
      ftime=unit_conversion_s(utime);
      if (0.==ftime) {
	SIMPUT_ERROR("unknown units of 'TIME' column");
	*status=EXIT_FAILURE;
	break;
      }
    }
    // END of determine unit conversion factors.


    // Read the header keywords.
    char comment[SIMPUT_MAXSTR];
    // Required keywords.
    fits_read_key(fptr, TDOUBLE, "MJDREF",   &lc->mjdref,   comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'MJDREF' from light curve '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      break;
    }

    fits_read_key(fptr, TDOUBLE, "TIMEZERO", &lc->timezero, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'TIMEZERO' from light curve '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (cphase>0) {
      // Only for periodic light curves.
      fits_read_key(fptr, TDOUBLE, "PHASE0", &lc->phase0, comment, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not read FITS keyword 'PHASE0' from light curve '%s'",
		filename);
	SIMPUT_ERROR(msg);
	break;
      }

      fits_read_key(fptr, TDOUBLE, "PERIOD", &lc->period, comment, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not read FITS keyword 'PERIOD' from light curve '%s'",
		filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // Verify that the specified period is positive.
      if (lc->period<=0.) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "period of light curve '%s' does not have a positive value",
		filename);
	SIMPUT_ERROR(msg);
	*status=EXIT_FAILURE;
	break;
      }

      int status2=EXIT_SUCCESS;
      fits_write_errmark();
      fits_read_key(fptr, TDOUBLE, "DPERIOD", &lc->dperiod, comment, &status2);
      fits_clear_errmark();
      if (EXIT_SUCCESS!=status2) {
	// The keyword 'DPERIOD' does not exist. Therefore we can set
	// the respective attribute to 0.
	lc->dperiod=0.;
      } else {
	// Check if dperiod is too small such that the applied formula
	// that is used for handling the period change, is guaranteed
	// to work properly.
	if (fabs(lc->dperiod)<0) {
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "value of DPERIOD of the light curve '%s' is negative (%e)",
		  filename, lc->dperiod);
	  SIMPUT_ERROR(msg);
	  *status=EXIT_FAILURE;
	  break;
	}
      }

    } else { // END of if (cphase>0) (periodic light curves)
      lc->phase0 =0.;
      lc->period =0.;
      lc->dperiod=0.;
    }

    // Optional keywords.
    opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_read_key(fptr, TFLOAT, "FLUXSCAL", &lc->fluxscal, comment, &opt_status);
    if (EXIT_SUCCESS!=opt_status) {
      // FLUXSCAL is not given in the FITS header. We therefore assume
      // that it has a value of 1.
      lc->fluxscal=1.;
      opt_status=EXIT_SUCCESS;
    }
    fits_clear_errmark();

    // Make sure that FLUXSCAL is positive.
    if (lc->fluxscal<=0.0) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "light curve '%s' has zero or negative value of keyword 'FLUXSCAL'",
	      filename);
      SIMPUT_ERROR(msg);
      *status=EXIT_FAILURE;
      break;
    }

    // Determine the number of rows in the table.
    fits_get_num_rows(fptr, &lc->nentries, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine number of entries in light curve");
      break;
    }

    // We correct the number of entries in the stored light curve.
    // Two correcions may be needed:
    // 1. Usually we will have light curves were phase 1.0 is not explicitly given,
    //    so an extra entry has to be added to store the first phase again at the end
    //    of the array so to be able to interpolate along the whole light curve.
    //    If it is not the case, we won't add any extra entry by this reason.
    // 2. We can also have a light curve were phase 0.0 is not given, so an extra entry
    //    has to be added at the begining so to fulfill the begining of the lightcurve.
    int nfitsentries=lc->nentries;
    int anynul=0;

    // Flags indicating wether the first phase is 0.0 and the last one is 1.0
    int first_phase_no0=0;
    int last_phase_no1=0;

    if (cphase>0) {
      // Correction 1.) Unless phase 1.0 is explicitly given, we add an extra row
      // where we will store the first phase again.
      // We first check that the last phase is not phase 1.0:
      double *last_phase=(double*)malloc(sizeof(double));
      fits_read_col(fptr, TDOUBLE, cphase, nfitsentries, 1, 1,
                    0, last_phase, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
        SIMPUT_ERROR("failed reading phase values from light curve");
        break;
      }
      if (*last_phase!=1.0){
        // phase 1.0 is not given, so we add a new entry
        lc->nentries++;
        last_phase_no1=1;
      }

      // Correction 2.) There is an extra correction, corresponding to the case
      // where phase0.0 spectrum is not given. In this situation, the
      // last phase given in the fits file should be taken into account
      // in first place to generate photons in the first part of the lc.
      double *first_phase=(double*)malloc(sizeof(double));
      fits_read_col(fptr, TDOUBLE, cphase, 1, 1, 1,
                    0, first_phase, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
        SIMPUT_ERROR("failed reading phase values from light curve");
        break;
      }
      if (*first_phase!=0.0){
	// phase 0.0 is not given, so we add a new entry
	lc->nentries++;
	first_phase_no0=1;
      }
    }
    // Now lc->nentries does not stand for the number of entries in the
    // fits table any more, but the number of entries in the loaded lightcurve

    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "light curve '%s' contains %ld data points\n",
	    filename, lc->nentries);
    SIMPUT_INFO(msg);

    // Allocate memory for the arrays.

    if (ctime>0) {
      lc->time =(double*)malloc(lc->nentries*sizeof(double));
      CHECK_NULL_BREAK(lc->time, *status,
                       "memory allocation for light curve failed");
    }
    if (cphase>0) {
      lc->phase=(double*)malloc(lc->nentries*sizeof(double));
      CHECK_NULL_BREAK(lc->phase, *status,
                       "memory allocation for light curve failed");
    }
    lc->flux   =(float*)malloc(lc->nentries*sizeof(float));
    CHECK_NULL_BREAK(lc->flux, *status,
                     "memory allocation for light curve failed");
    if (cspectrum>0) {
      lc->spectrum=(char**)malloc(lc->nentries*sizeof(char*));
      CHECK_NULL_BREAK(lc->spectrum, *status,
                       "memory allocation for light curve failed");
      // String buffer.

      spectrum[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      CHECK_NULL_BREAK(spectrum[0], *status,
                       "memory allocation for string buffer failed");
    }
    if (cimage>0) {
      lc->spectrum=(char**)malloc(lc->nentries*sizeof(char*));
      CHECK_NULL_BREAK(lc->spectrum, *status,
                       "memory allocation for light curve failed");
      // String buffer.

      image[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      CHECK_NULL_BREAK(image[0], *status,
                       "memory allocation for string buffer failed");
    }

    // Read the data from the table.
    // TIME
    if (ctime>0) {
      fits_read_col(fptr, TDOUBLE, ctime, 1, 1, nfitsentries,
                    0, lc->time, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
        SIMPUT_ERROR("failed reading time values from light curve");
        break;
      }
      // Multiply with unit scaling factor.
      long row;
      for (row=0; row<nfitsentries; row++) {
        lc->time[row]*=ftime;
      }
    }

    // PHASE
    if (cphase>0) {
      // We take the phases in the fits file. They are allocated at the
      // begining of the array lc->phase.
      fits_read_col(fptr, TDOUBLE, cphase, 1, 1, nfitsentries,
                    0, lc->phase, &anynul, status);
      // In the case we have no phase 1.0, we add the first one at the end:
      if (last_phase_no1){
	lc->phase[lc->nentries-1]=(lc->phase[0]+1.);
      }
      // In case that phase0.0 is not given, we need an extra correction:
      if (first_phase_no0){
	// We shift the phases one cell:
	for (int ii=0;ii<nfitsentries;ii++){
	  lc->phase[lc->nentries-2-ii]=lc->phase[lc->nentries-3-ii];
	}
	// We add the last phase given in the fits at the begining:
	lc->phase[0]=(lc->phase[lc->nentries-2]-1.);
      }

      if (EXIT_SUCCESS!=*status) {
        SIMPUT_ERROR("failed reading phase values from light curve");
        break;
	}
    }

    // FLUX
    // We take the phases in the fits file:
    fits_read_col(fptr, TFLOAT, cflux, 1, 1, nfitsentries,
                  0, lc->flux, &anynul, status);
    if (cphase>0) {
      // In the case we have no phase 1.0, we add the first one at the end:
      if(last_phase_no1){
	lc->flux[lc->nentries-1]=lc->flux[0];
      }
      //In case that phase0.0 is not given, we need an extra correction:
      if (first_phase_no0){
        // We shift the fluxes one cell:
        for (int ii=0;ii<nfitsentries;ii++){
          lc->flux[lc->nentries-2-ii]=lc->flux[lc->nentries-3-ii];
        }
        // We add the last flux given in the fits at the begining:
        lc->flux[0]=(lc->flux[lc->nentries-2]);
      }
    }
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading flux values from light curve");
      break;
    }

    // SPECTRUM
    if (cspectrum>0) {
      long row;
      // We take the phases in the fits file:
      for (row=0; row<nfitsentries; row++) {
        fits_read_col(fptr, TSTRING, cspectrum, row+1, 1, 1,
                      "", spectrum, &anynul, status);
        if (EXIT_SUCCESS!=*status) {
          SIMPUT_ERROR("failed reading spectrum references from light curve");
          break;
        }

	// We allocate spectra accordingly whether the first phase in
	// the fits file is phase 0.0 or not:
	if (first_phase_no0){
	  lc->spectrum[row+1]=
            (char*)malloc((strlen(spectrum[0])+1)*sizeof(char));
          CHECK_NULL_BREAK(lc->spectrum[row+1], *status,
                           "memory allocation for spectrum string failed");
	          strcpy(lc->spectrum[row+1], spectrum[0]);
	} else {
	  lc->spectrum[row]=
	    (char*)malloc((strlen(spectrum[0])+1)*sizeof(char));
	  CHECK_NULL_BREAK(lc->spectrum[row], *status,
			   "memory allocation for spectrum string failed");
	  strcpy(lc->spectrum[row], spectrum[0]);
	}

	// Now we add the extra spectrums at the end and/or at the begining
	if(row==0 && last_phase_no1){
	  lc->spectrum[lc->nentries-1]=(char*)malloc((strlen(spectrum[0])+1)*sizeof(char));
	  CHECK_NULL_BREAK(lc->spectrum[lc->nentries-1], *status,
			   "memory allocation for spectrum string failed");
	  strcpy(lc->spectrum[lc->nentries-1], spectrum[0]);
	}
	else if (row==nfitsentries-1 && first_phase_no0){
	  lc->spectrum[0]=(char*)malloc((strlen(spectrum[0])+1)*sizeof(char));
          CHECK_NULL_BREAK(lc->spectrum[0], *status,
                           "memory allocation for spectrum string failed");
          strcpy(lc->spectrum[0], spectrum[0]);
	}
      }


      CHECK_STATUS_BREAK(*status);
    }


    // IMAGE
    if (cimage>0) {
      long row;
      // We take the images in the fits file:
      for (row=0; row<nfitsentries; row++) {
        fits_read_col(fptr, TSTRING, cimage, row+1, 1, 1,
                      "", image, &anynul, status);
        if (EXIT_SUCCESS!=*status) {
          SIMPUT_ERROR("failed reading image references from light curve");
          break;
        }

        lc->image[row]=
          (char*)malloc((strlen(image[0])+1)*sizeof(char));
        CHECK_NULL_BREAK(lc->image[row], *status,
                         "memory allocation for image string failed");
        strcpy(lc->image[row], image[0]);
	// Now we add the first image again at the end:
	if(row==0 && last_phase_no1){
	  lc->image[lc->nentries-1]=
	    (char*)malloc((strlen(image[0])+1)*sizeof(char));
	  CHECK_NULL_BREAK(lc->image[row], *status,
			   "memory allocation for image string failed");
	  strcpy(lc->image[lc->nentries-1], image[0]);
	}
      }
      CHECK_STATUS_BREAK(*status);
    }
    // END of reading the data from the FITS table.

    // Store the file reference to the timing extension for later comparisons.
    lc->fileref=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_BREAK(lc->fileref, *status,
		     "memory allocation for file reference failed");
    strcpy(lc->fileref, filename);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=spectrum[0]) {
    free(spectrum[0]);
    spectrum[0]=NULL;
  }
  if (NULL!=image[0]) {
    free(image[0]);
    image[0]=NULL;
  }

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, lc);

  return(lc);
}


void saveSimputLC(SimputLC* const lc, const char* const filename,
		  char* const extname, int extver, int* const status)
{
  fitsfile* fptr=NULL;

  // String buffers.
  char* spectrum[1]={NULL};
  char* image[1]={NULL};

  int ncolumns=0;
  char **ttype=NULL;
  char **tform=NULL;
  char **tunit=NULL;

  do { // Error handling loop.

    // Check if the given light curve either contains a time
    // or a phase column, but not both.
    if ((NULL==lc->time) && (NULL==lc->phase)) {
      SIMPUT_ERROR("light curve does not contain TIME or PHASE column");
      *status=EXIT_FAILURE;
      break;
    }
    if ((NULL!=lc->time) && (NULL!=lc->phase)) {
      SIMPUT_ERROR("light curve contains both TIME and PHASE column");
      *status=EXIT_FAILURE;
      break;
    }
    if (NULL==lc->flux) {
      SIMPUT_ERROR("light curve does not contain FLUX column");
      *status=EXIT_FAILURE;
      break;
    }

    // Check if the EXTNAME has been specified.
    if (NULL==extname) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }
    if (0==strlen(extname)) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }

    // Check if the specified file exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "check failed whether file '%s' exists", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not open file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // Try to move to the specified extension.
      int status2=EXIT_SUCCESS;
      fits_write_errmark();
      fits_movnam_hdu(fptr, BINARY_TBL, extname, extver, &status2);
      fits_clear_errmark();
      if (BAD_HDU_NUM!=status2) {
	// If that works, the extension already exists.
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "extension '%s' with EXTVER=%d already exists",
		extname, extver);
	SIMPUT_ERROR(msg);
	*status=EXIT_FAILURE;
	break;
      }

    } else {
      // If no, create a new file.
      fits_create_file(&fptr, filename, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    // END of check, whether the specified file exists.

    // Create a new binary table.
    // Determine the number of columns.
    ncolumns=2;
    if (NULL!=lc->spectrum) ncolumns++;
    if (NULL!=lc->image)    ncolumns++;
    // Allocate memory for the format strings.
    ttype=(char**)malloc(ncolumns*sizeof(char*));
    tform=(char**)malloc(ncolumns*sizeof(char*));
    tunit=(char**)malloc(ncolumns*sizeof(char*));
    CHECK_NULL_BREAK(ttype, *status,
		     "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tform, *status,
		     "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tunit, *status,
		     "memory allocation for string buffer failed");
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      ttype[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      tform[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      tunit[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      CHECK_NULL_BREAK(ttype[ii], *status,
		       "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tform[ii], *status,
		       "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tunit[ii], *status,
		       "memory allocation for string buffer failed");
    }
    CHECK_STATUS_BREAK(*status);

    // Set up the table format.
    int ctime=0, cphase=0, cflux=0, cspectrum=0, cimage=0;
    if (NULL!=lc->time) {
      ctime=1;
      strcpy(ttype[0], "TIME");
      strcpy(tform[0], "D");
      strcpy(tunit[0], "s");
    } else {
      cphase=1;
      strcpy(ttype[0], "PHASE");
      strcpy(tform[0], "E");
      strcpy(tunit[0], "");
    }
    cflux=2;
    strcpy(ttype[1], "FLUX");
    strcpy(tform[1], "E");
    strcpy(tunit[1], "");
    if (NULL!=lc->spectrum) {
      cspectrum=3;
      strcpy(ttype[2], "SPECTRUM");
      strcpy(tform[2], "");
      strcpy(tunit[2], "1PA");
    }
    if (NULL!=lc->image) {
      cimage=4;
      strcpy(ttype[3], "IMAGE");
      strcpy(tform[3], "");
      strcpy(tunit[3], "1PA");
    }

    // Create the table.
    fits_create_tbl(fptr, BINARY_TBL, 0, ncolumns,
		    ttype, tform, tunit, extname, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not create binary table for source catalog in file '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Write header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "LIGHTCURVE", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);
    fits_write_key(fptr, TDOUBLE, "MJDREF", &lc->mjdref, "", status);
    fits_write_key(fptr, TDOUBLE, "TIMEZERO", &lc->timezero, "", status);
    // The keyword FLUXSCAL is deprecated. Therefore it is only stored in the
    // output file if its value is different from 1.
    if (fabs(lc->fluxscal-1.0)>1.0e-6) {
      fits_write_key(fptr, TFLOAT,  "FLUXSCAL", &lc->fluxscal, "", status);
    }
    int periodic=0;
    if (cphase>0) {
      // Only for periodic light curves.
      periodic=1;
      fits_write_key(fptr, TDOUBLE, "PHASE0" , &lc->phase0 , "", status);
      fits_write_key(fptr, TDOUBLE, "PERIOD" , &lc->period , "", status);
      fits_write_key(fptr, TDOUBLE, "DPERIOD", &lc->dperiod, "", status);
    }
    fits_write_key(fptr, TINT,  "PERIODIC", &periodic, "", status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed writing FITS keywords in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Create new rows in the table and store the data of the spectrum in it.
    fits_insert_rows(fptr, 0, lc->nentries, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed appending new rows to light curve");
      break;
    }

    if (ctime>0) {
      fits_write_col(fptr, TDOUBLE, ctime, 1, 1, lc->nentries,
		     lc->time, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing time values to light curve");
	break;
      }
    } else {
      fits_write_col(fptr, TDOUBLE, cphase, 1, 1, lc->nentries,
		     lc->phase, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed writing phase values to light curve");
	break;
      }
    }

    fits_write_col(fptr, TFLOAT, cflux, 1, 1, lc->nentries,
		   lc->flux, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing flux values to light curve");
      break;
    }

    if (cspectrum>0) {
      long row;
      for (row=0; row<lc->nentries; row++) {
	fits_write_col(fptr, TSTRING, cspectrum, row+1, 1, 1,
		       &lc->spectrum[row], status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed writing spectrum reference to light curve");
	  break;
	}
      }
      CHECK_STATUS_BREAK(*status);
    }
    if (cimage>0) {
      long row;
      for (row=0; row<lc->nentries; row++) {
	fits_write_col(fptr, TSTRING, cimage, row+1, 1, 1,
		       &lc->image[row], status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed writing image reference to light curve");
	  break;
	}
      }
      CHECK_STATUS_BREAK(*status);
    }

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=spectrum[0]) {
    free(spectrum[0]);
    spectrum[0]=NULL;
  }
  if (NULL!=image[0]) {
    free(image[0]);
    image[0]=NULL;
  }

  if (NULL!=ttype) {
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      if (NULL!=ttype[ii]) free(ttype[ii]);
    }
    free(ttype);
  }
  if (NULL!=tform) {
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      if (NULL!=tform[ii]) free(tform[ii]);
    }
    free(tform);
  }
  if (NULL!=tunit) {
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      if (NULL!=tunit[ii]) free(tunit[ii]);
    }
    free(tunit);
  }

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


SimputPSD* loadSimputPSD(const char* const filename, int* const status)
{
  // Get an empty SimputPSD data structure.
  SimputPSD* psd=newSimputPSD(status);
  CHECK_STATUS_RET(*status, psd);

  // Open the specified FITS file. The filename must uniquely identify
  // the PSD contained in a binary table via the extended filename
  // syntax. Therefore we do not have to care about the HDU number.
  fitsfile* fptr=NULL;
  fits_open_table(&fptr, filename, READONLY, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not open FITS table in file '%s'", filename);
    SIMPUT_ERROR(msg);
    return(psd);
  }

  do { // Error handling loop.

    // Get the column names.
    int cfrequency=0, cpower=0;
    // Required columns:
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "FREQUENCY", &cfrequency, status);
    fits_clear_errmark();
    if (EXIT_SUCCESS!=*status) {
      *status=EXIT_SUCCESS;
      fits_get_colnum(fptr, CASEINSEN, "FREQUENC", &cfrequency, status);
    }
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'FREQUENCY' in PSD '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    fits_get_colnum(fptr, CASEINSEN, "POWER", &cpower, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not find column 'POWER' in PSD '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Determine the unit conversion factors.
    float ffrequency=0.;
    char ufrequency[SIMPUT_MAXSTR];
    read_unit(fptr, cfrequency, ufrequency, status);
    CHECK_STATUS_BREAK(*status);
    ffrequency = unit_conversion_Hz(ufrequency);
    if (0.==ffrequency) {
      SIMPUT_ERROR("unknown units of 'FREQUENCY' column");
      *status=EXIT_FAILURE;
      break;
    }

    float fpower=0.;
    char upower[SIMPUT_MAXSTR];
    read_unit(fptr, cpower, upower, status);
    CHECK_STATUS_BREAK(*status);
    fpower = unit_conversion_s(upower);
    if (0.==fpower) {
      SIMPUT_ERROR("unknown units of 'POWER' column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.


    // Determine the number of rows in the table.
    fits_get_num_rows(fptr, &psd->nentries, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("could not determine number of entries in PSD");
      break;
    }

    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "PSD '%s' contains %ld data points\n",
	    filename, psd->nentries);
    SIMPUT_INFO(msg);

    // Allocate memory for the arrays.
    psd->frequency=(float*)malloc(psd->nentries*sizeof(float));
    CHECK_NULL_BREAK(psd->frequency, *status,
		     "memory allocation for PSD failed");
    psd->power    =(float*)malloc(psd->nentries*sizeof(float));
    CHECK_NULL_BREAK(psd->power, *status,
		     "memory allocation for PSD failed");

    // Read the data from the table.
    int anynul=0;
    // FREQUENCY
    fits_read_col(fptr, TFLOAT, cfrequency, 1, 1, psd->nentries,
		  0, psd->frequency, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading 'FREQUENCY' from PSD");
      break;
    }

    // POWER
    fits_read_col(fptr, TFLOAT, cpower, 1, 1, psd->nentries,
		  0, psd->power, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading 'POWER' from PSD");
      break;
    }

    // END of reading the data from the FITS table.

    // Store the file reference to the PSD for later comparisons.
    psd->fileref=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_BREAK(psd->fileref, *status,
		     "memory allocation for file reference failed");
    strcpy(psd->fileref, filename);

  } while(0); // END of error handling loop.

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, psd);

  return(psd);
}


void saveSimputPSD(SimputPSD* const psd, const char* const filename,
		   char* const extname, int extver, int* const status)
{
  fitsfile* fptr=NULL;

  // String buffers.
  int ncolumns=0;
  char **ttype=NULL;
  char **tform=NULL;
  char **tunit=NULL;

  do { // Error handling loop.

    // Check if the EXTNAME has been specified.
    if (NULL==extname) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }
    if (0==strlen(extname)) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }

    // Check if the specified file exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "check failed whether file '%s' exists", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not open file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // Try to move to the specified extension.
      int status2=EXIT_SUCCESS;
      fits_write_errmark();
      fits_movnam_hdu(fptr, BINARY_TBL, extname, extver, &status2);
      fits_clear_errmark();
      if (BAD_HDU_NUM!=status2) {
	// If that works, the extension already exists.
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "extension '%s' with EXTVER=%d already exists",
		extname, extver);
	SIMPUT_ERROR(msg);
	*status=EXIT_FAILURE;
	break;
      }

    } else {
      // If not, create a new file.
      fits_create_file(&fptr, filename, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    // END of check, whether the specified file exists.

    // Create a new binary table.
    // Determine the number of columns.
    ncolumns=2;
    // Allocate memory for the format strings.
    ttype=(char**)malloc(ncolumns*sizeof(char*));
    tform=(char**)malloc(ncolumns*sizeof(char*));
    tunit=(char**)malloc(ncolumns*sizeof(char*));
    CHECK_NULL_BREAK(ttype, *status,
		     "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tform, *status,
		     "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tunit, *status,
		     "memory allocation for string buffer failed");
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      ttype[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      tform[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      tunit[ii]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      CHECK_NULL_BREAK(ttype[ii], *status,
           "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tform[ii], *status,
           "memory allocation for string buffer failed");
      CHECK_NULL_BREAK(tunit[ii], *status,
           "memory allocation for string buffer failed");
    }
    CHECK_STATUS_BREAK(*status);

    // Set up the table format.
    int cfreq = 1, cpower = 2;
      strcpy(ttype[0], "FREQUENCY");
      strcpy(tform[0], "1E");
      strcpy(tunit[0], "Hz");
      strcpy(ttype[1], "POWER");
      strcpy(tform[1], "1E");
      strcpy(tunit[1], "1/Hz");

    // Create the table.
    fits_create_tbl(fptr, BINARY_TBL, 0, ncolumns,
        ttype, tform, tunit, extname, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not create binary table for spectrum in file '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Write header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "POWSPEC", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed writing FITS keywords in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Create new rows in the table and store the data of the power
    // spectrum in it.
    fits_insert_rows(fptr, 0, psd->nentries, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed appending new rows to PSD");
      break;
    }

    fits_write_col(fptr, TFLOAT, cfreq, 1, 1, psd->nentries,
		   psd->frequency, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing frequency values to PSD");
      break;
    }

    fits_write_col(fptr, TFLOAT, cpower, 1, 1, psd->nentries,
		   psd->power, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed writing power values to PSD");
      break;
    }

  } while(0); // END of error handling loop.

  // Release allocated memory.

  if (NULL!=ttype) {
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      if (NULL!=ttype[ii]) free(ttype[ii]);
    }
    free(ttype);
  }
  if (NULL!=tform) {
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      if (NULL!=tform[ii]) free(tform[ii]);
    }
    free(tform);
  }
  if (NULL!=tunit) {
    int ii;
    for (ii=0; ii<ncolumns; ii++) {
      if (NULL!=tunit[ii]) free(tunit[ii]);
    }
    free(tunit);
  }

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


SimputImg* loadSimputImg(const char* const filename, int* const status)
{
  // Image input buffer.
  double* image1d=NULL;

  // String buffer for FITS header.
  char* headerstr=NULL;

  // File pointer.
  fitsfile* fptr=NULL;

  // Get an empty SimputImg data structure.
  SimputImg* img=newSimputImg(status);
  CHECK_STATUS_RET(*status, img);

  do { // Error handling loop.

    // Open the specified FITS file. The filename must uniquely identify
    // the light curve contained in a binary table via the extended filename
    // syntax. Therefore we do not have to care about the HDU number.
    fits_open_image(&fptr, filename, READONLY, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not open FITS image in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Read the WCS header keywords.
    // Read the entire header into the string buffer.
    int nkeys;
    fits_hdr2str(fptr, 1, NULL, 0, &headerstr, &nkeys, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed reading FITS header of file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Parse the header string and store the data in the wcsprm data
    // structure.
    int nreject, nwcs;
    if (0!=wcspih(headerstr, nkeys, 0, 0, &nreject, &nwcs, &img->wcs)) {
      SIMPUT_ERROR("parsing of WCS header failed");
      *status=EXIT_FAILURE;
      break;
    }
    if (nreject>0) {
      SIMPUT_ERROR("parsing of WCS header failed");
      *status=EXIT_FAILURE;
      break;
    }

    // Determine the image dimensions.
    int naxis;
    fits_get_img_dim(fptr, &naxis, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed reading image dimensions in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    if (2!=naxis) {
      SIMPUT_ERROR("specified FITS HDU does not contain a 2-dimensional image");
      *status=EXIT_FAILURE;
      break;
    }

    long naxes[2];
    fits_get_img_size(fptr, naxis, naxes, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed reading image size in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    img->naxis1=naxes[0];
    img->naxis2=naxes[1];

    // Allocate memory for the image.
    img->dist=(double**)malloc(img->naxis1*sizeof(double*));
    CHECK_NULL_BREAK(img->dist, *status,
		     "memory allocation for source image failed");
    long ii;
    for (ii=0; ii<img->naxis1; ii++) {
      img->dist[ii] = (double*)malloc(img->naxis2*sizeof(double));
      CHECK_NULL_BREAK(img->dist[ii], *status,
		       "memory allocation for source image failed");
    }
    CHECK_STATUS_BREAK(*status);

    // Allocate memory for the image input buffer.
    image1d=(double*)malloc(img->naxis1*img->naxis2*sizeof(double));
    CHECK_NULL_BREAK(image1d, *status,
		     "memory allocation for source image buffer failed");

    // Read the image from the file.
    int anynul;
    double null_value=1e-308;   // lowest possible value
    long fpixel[2]={1, 1};   // Lower left corner.
    //              |--|--> FITS coordinates start at (1,1).
    long lpixel[2]={img->naxis1, img->naxis2}; // Upper right corner.
    long inc[2]   ={1, 1};
    fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value,
		     image1d, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed reading image from file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure and generate a probability distribution function,
    // i.e., sum up the pixels.
    double sum=0.;
    for(ii=0; ii<img->naxis1; ii++) {
    	long jj;
    	for(jj=0; jj<img->naxis2; jj++) {
    		sum+=image1d[ii+ img->naxis1*jj];
    		img->dist[ii][jj]=sum;
    	}
    }

    // Store the file reference to the image for later comparisons.
    img->fileref=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_BREAK(img->fileref, *status,
		     "memory allocation for file reference failed");
    strcpy(img->fileref, filename);

  } while(0); // END of error handling loop.

  // Release memory for buffer.
  if (NULL!=image1d) free(image1d);
  if (NULL!=headerstr) free(headerstr);

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, img);

  return(img);
}


void saveSimputImg(SimputImg* const img,
		   const char* const filename,
		   char* const extname, int extver,
		   int* const status)
{
  fitsfile* fptr=NULL;

  // Buffers.
  double* image1d=NULL;
  char* headerstr=NULL;

  do { // Error handling loop.

    // Check if the EXTNAME has been specified.
    if (NULL==extname) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }
    if (0==strlen(extname)) {
      SIMPUT_ERROR("EXTNAME not specified");
      *status=EXIT_FAILURE;
      break;
    }

    // Check if the specified file exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "check failed whether file '%s' exists", filename);
	SIMPUT_ERROR(msg);
	break;
      }

      // Try to move to the specified extension.
      int status2=EXIT_SUCCESS;
      fits_write_errmark();
      fits_movnam_hdu(fptr, IMAGE_HDU, extname, extver, &status2);
      fits_clear_errmark();
      if (BAD_HDU_NUM!=status2) {
	// If that works, the extension already exists.
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "extension '%s' with EXTVER=%d already exists",
		extname, extver);
	SIMPUT_ERROR(msg);
	*status=EXIT_FAILURE;
	break;
      }

    } else {
      // If no, create a new file.
      fits_create_file(&fptr, filename, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "could not create file '%s'", filename);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    // END of check, whether the specified file exists.


    // Allocate memory for the 1-dimensional image buffer
    // (required for output to FITS file).
    image1d=(double*)malloc(img->naxis1*img->naxis2*sizeof(double));
    CHECK_NULL_BREAK(image1d, *status,
		     "memory allocation for source image buffer failed");

    // Store the source image in the 1-dimensional buffer to handle it
    // to the FITS routine.
    long ii;
    for (ii=0; ii<img->naxis1; ii++) {
      long jj;
      for (jj=0; jj<img->naxis2; jj++) {
	image1d[ii+ img->naxis1*jj] = img->dist[ii][jj];
      }
    }

    // The data in the 1-dimensional image buffer still represents the
    // probability distribution stored in the image data structure.
    // However, in the output file we want to store the actual image,
    // NOT the distribution function. Therefore we have to inverted the
    // summing process.
    double sum=0.;
    for (ii=0; ii<img->naxis1; ii++) {
      long jj;
      for (jj=0; jj<img->naxis2; jj++) {
	double buffer = image1d[ii+ img->naxis1*jj];
	image1d[ii+ img->naxis1*jj] -= sum;
	sum = buffer;
      }
    }


    // Create a new FITS image.
    long naxes[2] = {img->naxis1, img->naxis2};
    fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status);
    //                                |-> naxis
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not create FITS image in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }
    // The image has been appended at the end if the FITS file.

    // Write header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "IMAGE", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
    fits_write_key(fptr, TSTRING, "EXTNAME", extname, "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed writing FITS keywords in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

    // Write WCS header keywords.
    int nkeyrec;
    if (0!=wcshdo(0, img->wcs, &nkeyrec, &headerstr)) {
      SIMPUT_ERROR("construction of WCS header failed");
      *status=EXIT_FAILURE;
      break;
    }

    char* strptr=headerstr;
    for (ii=0; ii<nkeyrec; ii++) {
    	if (strlen(strptr)==0) {
    		SIMPUT_ERROR("WCS header unexpectedly ended");
    		*status=EXIT_FAILURE;
    		break;
    	}

    	char strbuffer[81];
    	strncpy(strbuffer, strptr, 80);
    	strbuffer[80]='\0';
    	fits_write_record(fptr, strbuffer, status);
    	CHECK_STATUS_BREAK(*status);
    	strptr+=80;
    }

    CHECK_STATUS_BREAK(*status);

    // Store the image in the new extension.
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1)
    // Upper right corner.
    long lpixel[2] = {img->naxis1, img->naxis2};
    fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, image1d, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "failed writing FITS image in file '%s'", filename);
      SIMPUT_ERROR(msg);
      break;
    }

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=image1d) free(image1d);
  if (NULL!=headerstr) free(headerstr);

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


SimputPhList* openSimputPhList(const char* const filename,
			       const int mode,
			       int* const status)
{
  SimputPhList* phl=newSimputPhList(status);
  CHECK_STATUS_RET(*status, phl);

  // Open the photon list.
  fits_open_table(&phl->fptr, filename, mode, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not open FITS table in file '%s'", filename);
    SIMPUT_ERROR(msg);
    return(phl);
  }

  // Determine the total number of photons.
  fits_get_num_rows(phl->fptr, &phl->nphs, status);
  if (EXIT_SUCCESS!=*status) {
    SIMPUT_ERROR("could not determine number of entries in photon list");
    return(phl);
  }

  // Determine the reference RA and DEC values from the
  // header keywords.
  double refra=0., refdec=0.;
  char comment[SIMPUT_MAXSTR];
  fits_read_key(phl->fptr, TDOUBLE, "REFRA", &refra, comment, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not read FITS keyword 'REFRA' from photon list '%s'",
	    filename);
    SIMPUT_ERROR(msg);
    return(phl);
  }
  fits_read_key(phl->fptr, TDOUBLE, "REFDEC", &refdec, comment, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not read FITS keyword 'REFRA' from photon list '%s'",
	    filename);
    SIMPUT_ERROR(msg);
    return(phl);
  }
  if ((0.!=refra)||(0.!=refdec)) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("in current implementation photon lists must have "
		 "REFRA=0.0 and REFDEC=0.0");
    return(phl);
  }

  // Determine the column numbers.
  fits_get_colnum(phl->fptr, CASEINSEN, "RA", &phl->cra, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not find column 'RA' in photon list '%s'", filename);
    SIMPUT_ERROR(msg);
    return(phl);
  }

  fits_get_colnum(phl->fptr, CASEINSEN, "DEC", &phl->cdec, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not find column 'DEC' in photon list '%s'", filename);
    SIMPUT_ERROR(msg);
    return(phl);
  }

  fits_get_colnum(phl->fptr, CASEINSEN, "ENERGY", &phl->cenergy, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not find column 'ENERGY' in photon list '%s'", filename);
    SIMPUT_ERROR(msg);
    return(phl);
  }

  // Determine the index of the optional TIME column.
  int opt_status=EXIT_SUCCESS;
  fits_write_errmark();
  fits_get_colnum(phl->fptr, CASEINSEN, "TIME", &phl->ctime, &opt_status);
  if (EXIT_SUCCESS!=opt_status) {
    phl->ctime=0;
    opt_status=EXIT_SUCCESS;
  }
  fits_clear_errmark();

  // Determine the unit conversion factors.
  char ura[SIMPUT_MAXSTR];
  read_unit(phl->fptr, phl->cra, ura, status);
  CHECK_STATUS_RET(*status, phl);
  phl->fra=unit_conversion_rad(ura);
  if (0.==phl->fra) {
    SIMPUT_ERROR("unknown units of 'RA' column");
    *status=EXIT_FAILURE;
    return(phl);
  }

  char udec[SIMPUT_MAXSTR];
  read_unit(phl->fptr, phl->cdec, udec, status);
  CHECK_STATUS_RET(*status, phl);
  phl->fdec=unit_conversion_rad(udec);
  if (0.==phl->fdec) {
    SIMPUT_ERROR("unknown units of 'DEC' column");
    *status=EXIT_FAILURE;
    return(phl);
  }

  char uenergy[SIMPUT_MAXSTR];
  read_unit(phl->fptr, phl->cenergy, uenergy, status);
  CHECK_STATUS_RET(*status, phl);
  phl->fenergy=unit_conversion_keV(uenergy);
  if (0.==phl->fenergy) {
    SIMPUT_ERROR("unknown units of 'ENERGY' column");
    *status=EXIT_FAILURE;
    return(phl);
  }

  if (phl->ctime>0) {
    char utime[SIMPUT_MAXSTR];
    read_unit(phl->fptr, phl->ctime, utime, status);
    CHECK_STATUS_RET(*status, phl);
    phl->ftime=unit_conversion_s(utime);
    if (0.==phl->ftime) {
      SIMPUT_ERROR("unknown units of 'TIME' column");
      *status=EXIT_FAILURE;
      return(phl);
    }
  }

  // Determine timing keywords if TIME column is present.
  if (phl->ctime>0) {
    char comment[SIMPUT_MAXSTR];
    fits_read_key(phl->fptr, TDOUBLE, "MJDREF", &phl->mjdref, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'MJDREF' from photon list '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      return(phl);
    }

    fits_read_key(phl->fptr, TDOUBLE, "TIMEZERO", &phl->timezero, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'TIMEZERO' from photon list '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      return(phl);
    }

    fits_read_key(phl->fptr, TDOUBLE, "TSTART", &phl->tstart, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'TSTART' from photon list '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      return(phl);
    }

    fits_read_key(phl->fptr, TDOUBLE, "TSTOP", &phl->tstop, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'TSTOP' from photon list '%s'",
	      filename);
      SIMPUT_ERROR(msg);
      return(phl);
    }

    // Check if the covered interval is > 0.
    if (phl->tstop<=phl->tstart) {
      *status=EXIT_FAILURE;
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "TSTART (%e) >= TSTOP (%e) in photon list '%s'",
	      phl->tstart, phl->tstop, filename);
      SIMPUT_ERROR(msg);
      return(phl);
    }
  }

  // Store the file reference to the photon lists for later use.
  phl->fileref=
    (char*)malloc((strlen(filename)+1)*sizeof(char));
  CHECK_NULL_RET(phl->fileref, *status,
		 "memory allocation for file reference failed",
		 phl);
  strcpy(phl->fileref, filename);

  return(phl);
}


int getSimputExtType(SimputCtlg* const cat,
		     const char* const filename,
		     int* const status)
{
  // Check if there is any reference at all.
  if (0==strlen(filename)) {
    return(EXTTYPE_NONE);
  }

  // Cut the filename after the first occurence of a ']'.
  // This is sufficient to specify the extension.
  // If we did not apply this cut, the extension type for
  // each individual entry in a spectrum HDU would require
  // an opening of the file. Therefore this routine would be
  // much slower.
  char fileref[SIMPUT_MAXSTR];
  strcpy(fileref, filename);
  char* firstbracket=strchr(fileref, ']');
  if (NULL!=firstbracket) {
    firstbracket[1]='\0';
  }

  // Keep an internal cache of extension types in order to avoid
  // continuous re-opening.

  // Search if the required extension is available in the storage.
  int type=searchSimputExttypeBuffer(cat->extbuff, fileref);
  if (EXTTYPE_NONE!=type) {
    return(type);
  }

  // The extension is not contained in the cache. Therefore
  // we have to open it and check the header keywords.
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, fileref, READONLY, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not open file '%s'", fileref);
    SIMPUT_ERROR(msg);
    return(EXTTYPE_NONE);
  }

  // In case the filename refers to the file only WITHOUT any
  // explicit specification of the extension, check if the
  // primary extension is empty.
  char comment[SIMPUT_MAXSTR];
  if (NULL==firstbracket) {
    int naxis;
    fits_read_key(fptr, TINT, "NAXIS", &naxis, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "could not read FITS keyword 'NAXIS' from file '%s'", fileref);
      SIMPUT_ERROR(msg);
    }
    // If the primary extension is empty, move to the next HDU.
    if (0==naxis) {
      fits_movabs_hdu(fptr, 2, NULL, status);
      if (EXIT_SUCCESS!=*status) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "failed moving to 2nd HDU in file '%s'", fileref);
	SIMPUT_ERROR(msg);
      }
    }
  }

  // Read the HDUCLAS1 and HDUCLAS2 header keywords.
  char hduclas1[SIMPUT_MAXSTR];
  char hduclas2[SIMPUT_MAXSTR];
  fits_read_key(fptr, TSTRING, "HDUCLAS1", &hduclas1, comment, status);
  if (EXIT_SUCCESS!=*status) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "could not read FITS keyword 'HDUCLAS1' from file '%s'", fileref);
    SIMPUT_ERROR(msg);
  }
  // (Don't do an error checking here! Otherwise the file
  //  might not be closed after an error occurred with reading
  //  the header keyword.)

  // Read optional header keyword (is not used in SIMPUT
  // version >= 1.1.0).
  int opt_status=EXIT_SUCCESS;
  fits_write_errmark();
  fits_read_key(fptr, TSTRING, "HDUCLAS2", &hduclas2, comment, &opt_status);
  fits_clear_errmark();
  if (EXIT_SUCCESS!=opt_status) {
    hduclas2[0]='\0';
    opt_status=EXIT_SUCCESS;
  }

  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, EXTTYPE_NONE);


  // Check for the different extension types.

  if
    // SIMPUT version 1.0.0.
    (((0==strcmp(hduclas1, "SIMPUT")) &&
      (0==strcmp(hduclas2, "SPECTRUM"))) ||
     // SIMPUT version 1.1.0.
     (0==strcmp(hduclas1, "SPECTRUM"))) {
    type=EXTTYPE_MIDPSPEC;
  }

  else if
    // SIMPUT version 1.0.0.
    (((0==strcmp(hduclas1, "SIMPUT")) &&
      (0==strcmp(hduclas2, "LIGHTCUR"))) ||
     // SIMPUT version 1.1.0.
     (0==strcmp(hduclas1, "LIGHTCURVE"))) {
    type=EXTTYPE_LC;
  }

  else if
    // SIMPUT version 1.0.0.
    (((0==strcmp(hduclas1, "SIMPUT")) &&
      (0==strcmp(hduclas2, "POWSPEC"))) ||
     // SIMPUT version 1.1.0.
     (0==strcmp(hduclas1, "POWSPEC"))) {
    type=EXTTYPE_PSD;
  }

  else if
    // SIMPUT version 1.0.0.
    (((0==strcmp(hduclas1, "SIMPUT")) &&
      (0==strcmp(hduclas2, "IMAGE"))) ||
     // SIMPUT version 1.1.0.
     (0==strcmp(hduclas1, "IMAGE"))) {
    type=EXTTYPE_IMAGE;
  }

  else if
    // SIMPUT version 1.0.0.
    (((0==strcmp(hduclas1, "SIMPUT")) &&
      (0==strcmp(hduclas2, "PHOTONS"))) ||
     // SIMPUT version 1.1.0.
     (0==strcmp(hduclas1, "PHOTONS"))) {
    type=EXTTYPE_PHLIST;
  }

  else {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "extension type '%s' (HDUCLAS1) not supported", hduclas1);
    SIMPUT_ERROR(msg);
    *status=EXIT_FAILURE;
    return(EXTTYPE_NONE);
  }

  // Store the extension type in the internal cache.
  insertSimputExttypeBuffer(&(cat->extbuff), fileref, type, status);
  CHECK_STATUS_RET(*status, EXTTYPE_NONE);

  return(type);
}



void read_isisSpec_fits_file(char *fname, SimputMIdpSpec* simputspec,
		char *ISISFile, float Emin, float Emax,
		float plFlux, float bbFlux, float flFlux, float rflFlux,
		int *status){
	// Allocate memory for the buffer.
	SimputMIdpSpec* simputspecbuffer=newSimputMIdpSpec(status); // free it
	CHECK_STATUS_VOID(*status);

	fitsfile* fptr=NULL;

	// Loop over the different components of the spectral model.
	long nrows=0;
	int ii;
	for (ii=0; ii<4; ii++) {

		// Determine the file name.
		char filename[SIMPUT_MAXSTR];
		sprintf(filename, "%s.spec%d", fname, ii);
		fits_open_table(&fptr, filename, READONLY, status);
		CHECK_STATUS_VOID(*status);

		// Load the data from the file.
		int anynull;
		if (0==ii) {
			// Determine the number of rows.
			fits_get_num_rows(fptr, &nrows, status);
			CHECK_STATUS_VOID(*status);

			// Allocate memory and initialize.
			simputspec->nentries=nrows;
			simputspec->energy=(float*)malloc(nrows*sizeof(float));
			CHECK_NULL_VOID(simputspec->energy, *status,
					"memory allocation failed");
			simputspec->fluxdensity=(float*)malloc(nrows*sizeof(float));
			CHECK_NULL_VOID(simputspec->fluxdensity, *status,
					"memory allocation failed");
			long jj;
			for (jj=0; jj<nrows; jj++) {
				simputspec->energy[jj]=0.;
				simputspec->fluxdensity[jj]=0.;
			}

			simputspecbuffer->nentries=nrows;
			simputspecbuffer->energy=(float*)malloc(nrows*sizeof(float));
			CHECK_NULL_BREAK(simputspecbuffer->energy, *status,
					"memory allocation failed");
			simputspecbuffer->fluxdensity=(float*)malloc(nrows*sizeof(float));
			CHECK_NULL_BREAK(simputspecbuffer->fluxdensity, *status,
					"memory allocation failed");
			for (jj=0; jj<nrows; jj++) {
				simputspecbuffer->energy[jj]=0.;
				simputspecbuffer->fluxdensity[jj]=0.;
			}

			// Read the energy column.
			fits_read_col(fptr, TFLOAT, 1, 1, 1, nrows, 0,
					simputspec->energy, &anynull, status);
			CHECK_STATUS_VOID(*status);
			for (jj=0; jj<nrows; jj++) {
				simputspecbuffer->energy[jj]=simputspec->energy[jj];
			}

		} else {
			// Check whether the number of entries is
			// consistent with the previous files.
			fits_get_num_rows(fptr, &nrows, status);
			CHECK_STATUS_VOID(*status);

			if (nrows!=simputspec->nentries) {
				SIMPUT_ERROR("inconsistent sizes of spectra");
				*status=EXIT_FAILURE;
				return;
			}
		}

		// Read the flux column.
		fits_read_col(fptr, TFLOAT, 2, 1, 1, nrows, 0,
				simputspecbuffer->fluxdensity, &anynull, status);
		CHECK_STATUS_VOID(*status);

		fits_close_file(fptr, status);
		fptr=NULL;
		CHECK_STATUS_VOID(*status);

		// If the spectrum is given via individual components, they
		// have to be normalized according to their respective fluxes.
		if (strlen(ISISFile)==0) {
			// Determine the required flux in the reference band.
			float shouldflux=0.;
			switch(ii) {
			case 0:
				shouldflux=plFlux;
				break;
			case 1:
				shouldflux=bbFlux;
				break;
			case 2:
				shouldflux=flFlux;
				break;
			case 3:
				shouldflux=rflFlux;
				break;
			default:
				*status=EXIT_FAILURE;
				break;
			}
			CHECK_STATUS_VOID(*status);

			float factor;
			if (shouldflux==0.) {
				factor=0.;
			} else {
				// Determine the factor between the actual flux in the reference
				// band and the required flux.
				factor=shouldflux/
				getSimputMIdpSpecBandFlux(simputspecbuffer,Emin, Emax);
			}

			// Add the normalized component to the total spectrum.
			if (factor>0.) {
				long jj;
				for (jj=0; jj<nrows; jj++) {
					simputspec->fluxdensity[jj]+=simputspecbuffer->fluxdensity[jj]*factor;
				}
			}

		} else {
			// The spectral model is given in an ISIS parameter file.
			// Therefore we do not have to normalize it, but can directly
			// add it to the SIMPUT spectrum.
			long jj;
			for (jj=0; jj<nrows; jj++) {
				simputspec->fluxdensity[jj]=simputspecbuffer->fluxdensity[jj];
			}

			// Since there are no further components, we can skip
			// the further processing of the loop.
			break;
		}
	}
	CHECK_STATUS_VOID(*status);
	// END of loop over the different spectral components.
}

static int check_xspec_linecont(char* buf){

	if (buf && *buf && buf[strlen(buf)-1]=='\n' && strlen(buf)>1 && buf[strlen(buf)-2]=='-' ){
		return 1;
	} else {
		return 0;
	}

}

void read_xspecSpec_file(char *fname, SimputMIdpSpec* simputspec, int *status){
	// The spectrum is contained in a .qdp file produced by XSPEC/PLT,
	// and has to be loaded from there.

	FILE* xspecfile=NULL;

	// Open the file.
	char filename[SIMPUT_MAXSTR];
	sprintf(filename, "%s.qdp", fname);
	xspecfile=fopen(filename, "r");
	CHECK_NULL_VOID(xspecfile, *status, "could not open XSPEC .qdp file");

	// Determine the number of rows.
	long nlines=0;
	char c=0;
	char prev_c=0;
	while(!feof(xspecfile)) {
		c=fgetc(xspecfile);
		// if we have a line continuatiion, length is one less
		if ('\n'==c && prev_c!='-'){
			nlines++;
		}
		prev_c = c;
	}
	// Check if the last line has been empty.
	if('\n'==c) {
		nlines--;
	}

	// The first 3 lines do not contain data.
	nlines-=3;

	// Allocate memory.
	simputspec->nentries=nlines;
	simputspec->energy=(float*)malloc(nlines*sizeof(float));
	CHECK_NULL_VOID(simputspec->energy, *status, "memory allocation failed");
	simputspec->fluxdensity=(float*)malloc(nlines*sizeof(float));
	CHECK_NULL_VOID(simputspec->energy, *status, "memory allocation failed");

	// Reset the file pointer, read the data, and store them in
	// the SimputMIdpSpec data structure.
	rewind(xspecfile);
	// Read the first three lines.
	char sbuffer1[SIMPUT_MAXSTR], sbuffer2[SIMPUT_MAXSTR];
	int ibuffer;
	if (fscanf(xspecfile, "%s %s %d\n", sbuffer1, sbuffer2, &ibuffer)<3) {
		SIMPUT_ERROR("failed reading data from ASCII file");
		*status=EXIT_FAILURE;
		return;
	}
	if (fscanf(xspecfile, "%s\n", sbuffer1)<1) {
		SIMPUT_ERROR("failed reading data from ASCII file");
		*status=EXIT_FAILURE;
		return;
	}
	if (fscanf(xspecfile, "%s\n", sbuffer1)<1) {
		SIMPUT_ERROR("failed reading data from ASCII file");
		*status=EXIT_FAILURE;
		return;
	}
	// Read the actual data.
	long ii=0;
	char linebuffer[SIMPUT_MAXSTR];
	int linecont=0;
	float fbuffer;

	while( (fgets(linebuffer, SIMPUT_MAXSTR, xspecfile)!=NULL) && (ii<nlines) ){

		// was there a line continuation before?
		if(linecont==0){
			if (sscanf(linebuffer, "%f %f %f",
					&(simputspec->energy[ii]),
					&fbuffer,
					&(simputspec->fluxdensity[ii]))!=3) {
				SIMPUT_ERROR("failed reading data from ASCII file");
				*status=EXIT_FAILURE;
				return;
			}
			ii++;
		}
		linecont=check_xspec_linecont(linebuffer);
	}

	CHECK_STATUS_VOID(*status);

	// Close the file.
	fclose(xspecfile);
	xspecfile=NULL;
}

void write_isisSpec_fits_file(char *fname, char *ISISFile, char *ISISPrep,
		char *ISISPostCmd, float Elow, float Eup, int nbins, int logegrid,
		float plPhoIndex, float bbkT, float flSigma, float rflSpin, float NH,
		int *status){

	FILE* cmdfile=NULL;
	char cmdfilename[L_tmpnam]="";

	// Open the ISIS command file.
	if (NULL==tmpnam(cmdfilename)) {
		SIMPUT_ERROR("failed getting temporary filename for ISIS command file");
		*status=EXIT_FAILURE;
		return;
	}
	cmdfile=fopen(cmdfilename,"w");
	CHECK_NULL_VOID(cmdfile, *status, "opening temporary file failed");

	// Write the header.
	fprintf(cmdfile, "require(\"isisscripts\");\n");
	fprintf(cmdfile, "()=xspec_abund(\"wilm\");\n");
	fprintf(cmdfile, "use_localmodel(\"relxill\");\n");

	// Define the energy grid.
	fprintf(cmdfile, "variable lo, hi; \n");

	if (logegrid) {
		fprintf(cmdfile, "(lo,hi) = log_grid(%e,%e,%i); ", Elow, Eup, nbins);
	} else {
		fprintf(cmdfile, "(lo,hi) = linear_grid(%e,%e,%i); ", Elow, Eup, nbins);
	}


	fprintf(cmdfile, "variable fluxdensity;\n");
	fprintf(cmdfile, "variable spec;\n");

	// Distinguish whether the individual spectral components or
	// an ISIS spectral parameter file should be used.
	if (strlen(ISISFile)==0) {

		// Loop over the different components of the spectral model.
		int ii;
		for (ii=0; ii<4; ii++) {

			// Define the spectral model and set the parameters.
			switch(ii) {
			case 0:
				fprintf(cmdfile, "fit_fun(\"tbabs(1)*powerlaw(1)\");\n");
				fprintf(cmdfile, "set_par(\"powerlaw(1).PhoIndex\", %e);\n",
						plPhoIndex);
				break;
			case 1:
				fprintf(cmdfile, "fit_fun(\"tbabs(1)*bbody(1)\");\n");
				fprintf(cmdfile, "set_par(\"bbody(1).kT\", %e);\n", bbkT);
				break;
			case 2:
				fprintf(cmdfile, "fit_fun(\"tbabs(1)*egauss(1)\");\n");
				fprintf(cmdfile, "set_par(\"egauss(1).center\", 6.4);\n");
				fprintf(cmdfile, "set_par(\"egauss(1).sigma\", %e);\n", flSigma);
				break;
			case 3:
				fprintf(cmdfile, "fit_fun(\"tbabs(1)*relline(1)\");\n");
				fprintf(cmdfile, "set_par(\"relline(1).lineE\", 6.4);\n");
				fprintf(cmdfile, "set_par(\"relline(1).a\", %f);\n", rflSpin);
				break;
			default:
				*status=EXIT_FAILURE;
				break;
			}
			CHECK_STATUS_VOID(*status);

			// Absorption is the same for all spectral components.
			fprintf(cmdfile, "set_par(\"tbabs(1).nH\", %e);\n", NH);

			// Evaluate the spectral model and store the data in a temporary
			// FITS file.
			fprintf(cmdfile, "fluxdensity=eval_fun_keV(lo, hi)/(hi-lo);\n");
			fprintf(cmdfile, "spec=struct{ENERGY=0.5*(lo+hi), FLUXDENSITY=fluxdensity};\n");
			fprintf(cmdfile,
					"fits_write_binary_table(\"%s.spec%d\",\"SPECTRUM\", spec);\n",
					fname, ii);
		}
		CHECK_STATUS_VOID(*status);
		// END of loop over the different spectral components.

	} else {
		// An ISIS parameter file with an explizit spectral
		// model is given.
		if(strlen(ISISPrep)>0){
			fprintf(cmdfile, "require(\"%s\");\n", ISISPrep);
		}
		fprintf(cmdfile, "load_par(\"%s\");\n", ISISFile);
		if(strlen(ISISPostCmd) > 0){
			fprintf(cmdfile, "%s\n",ISISPostCmd);
		}
		fprintf(cmdfile, "fluxdensity=eval_fun_keV(lo, hi)/(hi-lo);\n");
		fprintf(cmdfile, "print(sum(fluxdensity)); list_par;\n");
		fprintf(cmdfile, "spec=struct{ENERGY=0.5*(lo+hi), FLUXDENSITY=fluxdensity};\n");
		fprintf(cmdfile,
				"fits_write_binary_table(\"%s.spec0\",\"SPECTRUM\", spec);\n",
				fname);
	}
	// END of using an explicit spectral model given in an ISIS
	// parameter file.

	fprintf(cmdfile, "exit;\n");

	// End of writing the ISIS command file.
	fclose(cmdfile);
	cmdfile=NULL;

	// Construct the shell command to run ISIS.
	char command[SIMPUT_MAXSTR];
	strcpy(command, "isis ");
	strcat(command, cmdfilename);

	// Run ISIS.
	*status=system(command);
	CHECK_STATUS_VOID(*status);

}

void write_xspecSpec_file(char *fname, char *XSPECFile, char *XSPECPrep, char *XSPECPostCmd, float Elow,
		float Eup,	int nbins, int logegrid,  int *status){

	FILE* cmdfile=NULL;
	char cmdfilename[L_tmpnam]="";

   	// Open the command file.
    	if (NULL==tmpnam(cmdfilename)) {
    		SIMPUT_ERROR("failed getting temporary filename for Xspec command file");
    		*status=EXIT_FAILURE;
    		return;
    	}
    	cmdfile=fopen(cmdfilename,"w");
    	CHECK_NULL_VOID(cmdfile, *status, "opening temporary file failed");

    	// Write the header.
    	fprintf(cmdfile, "query yes\n");

        if ( (strlen(XSPECPrep)!=0) && (strcmp(XSPECPrep,"none")!=0) && (strcmp(XSPECPrep,"NONE")!=0) ) {
                fprintf(cmdfile, "@%s\n", XSPECPrep);
        }

    	fprintf(cmdfile, "@%s\n", XSPECFile);

    	if (strlen(XSPECPostCmd)!=0) {
        	fprintf(cmdfile, "%s\n", XSPECPostCmd);
    	}

    	if (logegrid){
    		fprintf(cmdfile, "dummyrsp %f %f %d log\n", Elow, Eup, nbins);
    	} else {
    		fprintf(cmdfile, "dummyrsp %f %f %d lin\n", Elow, Eup, nbins);
    	}
    	fprintf(cmdfile, "setplot device /null\n");

    	fprintf(cmdfile, "setplot command wdata %s.qdp\n", fname);
    	fprintf(cmdfile, "plot model\n");
    	fprintf(cmdfile, "quit\n");

    	// End of writing the command file.
    	fclose(cmdfile);
    	cmdfile=NULL;

    	// Construct the shell command to run Xspec.
    	char command[SIMPUT_MAXSTR];
    	strcpy(command, "xspec ");
    	strcat(command, cmdfilename);

    	// Run Xspec.
    	*status=system(command);
    	CHECK_STATUS_VOID(*status);
}



/*
 * Initialize the spectrum cache
 */
void initSpecCache()
{
  if ( SPEC_MAX_CACHE < 1 )
  {
    headas_chat(5, "Spectrum caching disabled at compile time\n");
    SpecCache = NULL;
    return ;
  }
  int status = 0;
  SpecCache = newSimputSpecExtCache(&status);
  if ( status != EXIT_SUCCESS )
  {
    SIMPUT_WARNING("Could not allocate spectrum cache");
    SpecCache = NULL;
  }
}


/*
 * Destroy the spectrum cache
 */
void destroySpecCache()
{
  destroySpecCacheBuff(SpecCache);
  SpecCache = NULL;
}

/*
 * Test whether a fits file is already openend in the cache
 */
long specIsCached(char *fname, char *extname, int extver)
{
  headas_chat(5, "Testing whether extension %s version %d of fitsfile %s is already in openend\n", extname, extver, fname);
  for (long ii=0; ii<SpecCache->n; ii++)
  {
    if (
	extver == SpecCache->extver[ii] &&
	strcmp(fname, SpecCache->filename[ii]) == 0 &&
	strcmp(extname, SpecCache->extname[ii]) == 0
	)
       {
	  return ii;
       }
  }
  return -1;
}


/*
 * function to sort the names column in the order strcmp, necessary for the bsearch
 * Note: The exact alignment of the strings afterwards is not known,
 * however, that doesn't matter because bsearch then also uses strcmp do find the
 * correct string again
 */
static int cmpNameCol(const void *a, const void *b)
{
  long i1 = *((long *) a);
  long i2 = *((long *) b);

  return strcmp( NamePtr[i1],
      NamePtr[i2] );
}

/*
 * Just a wrapper for a strcmp function in order to be compatible with bsearch
 */
static int myStrCmp(const void *s1, const void *s2) {
    return strcmp(*(char **) s1, *(char **) s2);
}


/*
 * Open a fits file with the specified extension in the nth spectrum 'slot'
 */
void openNthSpecCache(char *fname, char *extname, int extver, long n, int *status)
{
  headas_chat(5, "Opening extension %s of file %s into spectrum cache\n", extname, fname);

  long numrows;
  int namelen = 0;
  int anynull = 0;
  int typedecode;
  long nenergy, nfluxdensity, width;
  char uenergy[SIMPUT_MAXSTR], uflux[SIMPUT_MAXSTR];

  SpecCache->extver[n] = extver;

  SpecCache->filename[n] = (char *) malloc( (strlen(fname)+1) * sizeof(char));
  CHECK_NULL_VOID(SpecCache->filename[n], *status, "Could not allocate string buffer");
  strcpy(SpecCache->filename[n], fname);

  SpecCache->extname[n] = (char *) malloc( (strlen(extname)+1) * sizeof(char));
  CHECK_NULL_VOID(SpecCache->extname[n], *status, "Could not allocate string buffer");
  strcpy(SpecCache->extname[n], extname);

  headas_chat(5, "Opening %s\n", fname);
  fits_open_file(&(SpecCache->ext[n]), fname, READONLY, status);
  FITSERROR;

  headas_chat(5, "Moving ffptr to extension %s\n", extname);
  fits_movnam_hdu(SpecCache->ext[n], BINARY_TBL, extname, SpecCache->extver[n], status);
  FITSERROR;

  headas_chat(5, "Getting column number of the NAME column\n");
  fits_get_colnum(SpecCache->ext[n], CASEINSEN, "NAME", &(SpecCache->cname[n]), status);
  if ( *status )
  {
    headas_chat(5, "Column NAME in the spectrum extension %s of file %s doesn't exist, ignoring\n", extname, fname);
    headas_chat(5, "Resetting cfitsio error stack\n");
    fits_clear_errmsg();
    *status = EXIT_SUCCESS;
    SpecCache->cname[n] = 0;
  } else {
    int naxis;
    long naxes;
    fits_read_tdim(SpecCache->ext[n], SpecCache->cname[n], 2, &naxis, &naxes, status);
    FITSERROR;
    namelen = naxes+1;
    headas_chat(5, "Length of the names in this file: %ld\n", namelen);
  }

  headas_chat(5, "Getting column number of the ENERGY column\n");
  fits_get_colnum(SpecCache->ext[n], CASEINSEN, "ENERGY", &(SpecCache->cenergy[n]), status);
  FITSERROR;

  headas_chat(5, "Getting column number of the FLUXDENSITY column\n");
  fits_get_colnum(SpecCache->ext[n], CASEINSEN, "FLUXDENSITY", &(SpecCache->cflux[n]), status);
  if ( *status != EXIT_SUCCESS )
  {
    *status = EXIT_SUCCESS;
    fits_clear_errmsg();
    char *msg = (char *) malloc(SIMPUT_MAXSTR * sizeof(char));
    if ( msg == NULL )
    {
      fprintf(stderr, "Couldn't allocate memory for error message. That's dumb. I will die now\n");
      exit (EXIT_FAILURE);
    }
    sprintf(msg, "The FLUXDENSITY column in file %s extension %s is not called FLUXDENSITY, trying FLUX now", fname, extname);
    SIMPUT_WARNING(msg);
    free(msg);

    fits_get_colnum(SpecCache->ext[n], CASEINSEN, "FLUX", &(SpecCache->cflux[n]), status);
    FITSERROR;
    headas_chat(5, "Success\n");
  }

  headas_chat(5, "Getting number of entries in this spectrum extension:\n");
  fits_get_num_rows(SpecCache->ext[n], &numrows, status);
  FITSERROR;
  headas_chat(5, "%ld\n", numrows);
  SpecCache->nspec[n] = numrows;


  headas_chat(5, "Determining table format\n");

  fits_get_coltype(SpecCache->ext[n], SpecCache->cenergy[n], &typedecode,
      &nenergy, &width, status);
  FITSERROR;

  fits_get_coltype(SpecCache->ext[n], SpecCache->cflux[n], &typedecode,
      &nfluxdensity, &width, status);
  FITSERROR;

  // Variable length data type, different method, see the loadSimputMIdpSpec function
  if ( nenergy == 1 )
  {
    long offset;
    fits_read_descript(SpecCache->ext[n], SpecCache->cenergy[n], 1, &nenergy,
	&offset, status);
    FITSERROR;

    fits_read_descript(SpecCache->ext[n], SpecCache->cflux[n], 1, &nfluxdensity,
	&offset, status);
    FITSERROR;
  }

  if (nenergy!=nfluxdensity) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "number of energy and flux entries in spectrum '%s%s' is "
            "not equivalent", fname, extname);
    SIMPUT_ERROR(msg);
    *status=EXIT_FAILURE;
    return;
  }

  SpecCache->nbins[n] = nenergy;
  char msg[SIMPUT_MAXSTR];
  sprintf(msg, "spectrum '%s' contains %ld data points",
          fname, SpecCache->nbins[n]);
  SIMPUT_INFO(msg);

  headas_chat(5, "Determining conversion factors\n");
  read_unit(SpecCache->ext[n], SpecCache->cenergy[n], uenergy, status);
  FITSERROR;
  SpecCache->fenergy[n] = unit_conversion_keV(uenergy);
  if ( SpecCache->fenergy[n] == 0. )
  {
    SIMPUT_ERROR("unknown units of 'ENERGY' column");
    *status = EXIT_FAILURE;
    return ;
  }

  read_unit(SpecCache->ext[n], SpecCache->cflux[n], uflux, status);
  FITSERROR;
  SpecCache->fflux[n] = unit_conversion_keV(uenergy);
  if ( SpecCache->fflux[n] == 0. )
  {
    SIMPUT_ERROR("unknown units of 'FLUXDENSITY' column");
    *status = EXIT_FAILURE;
    return ;
  }


  if ( namelen )
  {
    headas_chat(5, "Allocating the namecol struct for %ld entries with %d chars\n", numrows, namelen);
    SpecCache->namecol[n] = newSpecNameCol(numrows, namelen, status);
    CHECK_NULL_VOID(SpecCache->namecol[n], *status, "Error creating the namecol struct");
  } else {
    SpecCache->namecol[n] = NULL;
  }

  // Actually reading the columns
  if ( namelen )
  {
    char **tmpstr;
    SpecCache->namecol[n]->n = numrows;
    headas_chat(5, "Reading the names column\n");
    fits_read_col(SpecCache->ext[n], TSTRING, SpecCache->cname[n], 1, 1, numrows, NULL, SpecCache->namecol[n]->name, &anynull, status);
    FITSERROR;

    // Fill up the rows for sorting and comparing
    for (long ii=0; ii<SpecCache->namecol[n]->n; ii++)
    {
      SpecCache->namecol[n]->row[ii] = ii;
    }

    headas_chat(5, "Sorting the names\n");
    NamePtr = SpecCache->namecol[n]->name;
    qsort(SpecCache->namecol[n]->row,
	SpecCache->namecol[n]->n,
	sizeof(long),
	&cmpNameCol);
    NamePtr = NULL;

    tmpstr = (char **) malloc(SpecCache->namecol[n]->n * sizeof(char *));
    CHECK_NULL_VOID(tmpstr, *status, "Could not allocate temporary string array");
    for (long ii=0; ii<SpecCache->namecol[n]->n; ii++)
    {
      tmpstr[ii] = (char *) malloc(SpecCache->namecol[n]->namelen * sizeof(char));
      CHECK_NULL_VOID(tmpstr[ii], *status, "Could not allocate temporary string array");
      strcpy(tmpstr[ii], SpecCache->namecol[n]->name[SpecCache->namecol[n]->row[ii]]);
      free(SpecCache->namecol[n]->name[SpecCache->namecol[n]->row[ii]]);
    }
    free(SpecCache->namecol[n]->name);
    SpecCache->namecol[n]->name = tmpstr;
  }
  if ( SpecCache->n < SPEC_MAX_CACHE  )
  {
    SpecCache->n++;
  }
  SpecCache->last = n;

/*
  printf("filenamen[%ld]: %s\n", n, SpecCache->filename[n]);
  printf("extname[%ld]: %s\n", n, SpecCache->extname[n]);
  printf("extver[%ld]: %d\n", n, SpecCache->extver[n]);
  printf("cname[%ld]: %d\n", n, SpecCache->cname[n]);
  printf("cflux[%ld]: %d\n", n, SpecCache->cflux[n]);
  printf("nspec[%ld]: %ld\n", n, SpecCache->nspec[n]);
  printf("fenergy[%ld]: %f\n", n, SpecCache->fenergy[n]);
  printf("fflux[%ld]: %f\n", n, SpecCache->fflux[n]);
  printf("nbins[%ld]: %ld\n", n, SpecCache->nbins[n]);
  printf("namecol:\n");
  for ( long ii=0; ii<SpecCache->namecol[n]->n; ii++)
  {
    printf("\tnamecol[%ld]->name[%ld]: %s\n", n, ii, SpecCache->namecol[n]->name[ii]);
    printf("\tnamecol[%ld]->row[%ld]: %ld\n", n, ii, SpecCache->namecol[n]->row[ii]);
  }
  */
}

/*
 * Function to scan the entire filename of a fits file, including the extended filename syntax,
 * and write the real filename, the extension name and the extension version referred to by the
 * extended filename syntax to the addresses pointed to respectively, and in the end return
 * a pointer to an allocated string containing the expression of the extended filename syntax
 * used to filter the row in the spectrum extension, NULL upon error
 */
char *scanSpecFileName(char *filename, char **basename, char **extname, int *extver, int *status)
{
  // Another string buffer
  char strbuff[SIMPUT_MAXSTR];
  // Temporary pointer
  char *tmpstr;
  // string holding an expression to work on
  char *expr;

  headas_chat(5, "\nObtaining the real filename: ");
  tmpstr = strndup(filename, strchr(filename, '[') - filename);
  CHECK_NULL_RET(tmpstr,*status,"scanSpecFileName: Error allocating tmpstr",NULL);

  *basename = (char *) malloc( (strlen(tmpstr) + 1) * sizeof(char));
  CHECK_NULL_RET(*basename, *status, "Error allocating string buffer", NULL);
  strcpy(*basename, tmpstr);
  headas_chat(5, "%s\n", *basename);

  expr = strndup( strchr(filename, '[')+1, strchr( strchr(filename, '[') , ']') - strchr(filename, '[') - 1 );

  CHECK_NULL_RET(expr,*status,"scanSpecFileName: Error allocating expr",NULL);

  headas_chat(5, "Testing whether the extension version is provided ... ");
  if ( strchr(expr, ',') != NULL ) {
     headas_chat(5, "yes\n");
     if ( sscanf( strchr(expr, ',') + 1, "%d", extver) == EOF ) {
	sprintf(strbuff, "Error scanning extended filesyntax for extension version in %s\n", filename);
	SIMPUT_ERROR(strbuff);
	return NULL;
     } else {
	headas_chat(5, "Extension version: %d\n", *extver);
	*(strchr(expr, ',')) = '\0';
	headas_chat(5, "Scanning for extension name\n");
	if ( sscanf(expr, "%s", strbuff) == EOF ) {
	   sprintf(strbuff, "Error scanning extended filesyntax for extension name in %s\n", filename);
	   SIMPUT_ERROR(strbuff);
	   return NULL;
	} else {
	   *extname = (char *) malloc( (strlen(strbuff) + 1) * sizeof(char));
	   CHECK_NULL_RET(*extname, *status, "Error allocating string buffer", NULL);
	   strcpy(*extname, strbuff);
	   headas_chat(5, "Extension name: %s\n", *extname);
	}
     }
  } else {
     headas_chat(5, "no, falling back to extension version 0\n");
     *extver = 0;
     headas_chat(5, "Scanning for extension name\n");
     if ( sscanf(expr, "%s", strbuff) == EOF ) {
	sprintf(strbuff, "Error scanning extended filesyntax for extension version in %s\n", filename);
	SIMPUT_ERROR(strbuff);
	return NULL;
     } else {
	*extname = (char *) malloc( (strlen(strbuff) + 1) * sizeof(char));
	CHECK_NULL_RET(*extname, *status, "Error allocating string buffer", NULL);
	strcpy(*extname, strbuff);
	headas_chat(5, "Extension name: %s\n", *extname);
     }
  }

  free(expr);
  free(tmpstr);

  char *expr2 = strchr(filename, '[');
  // shouldn't happen, but better safe than sorry
  CHECK_NULL_RET(expr2,*status,"File does not conform to extended syntax",NULL);
  expr2++; // move to start of selection string

  char *expr3=strchr(expr2,'[');

  char *retval;
  if (expr3==NULL) {
     // file has only one [] in the FITS selection string

     // bail out the string is [] or ] is missing,
     // (note C's shortcircuiting, check strlen first, so the
     // 2nd check will only happen if strlen>=1)
     if (strlen(expr2)<1 || expr2[strlen(expr2)-1]!=']') {
	SIMPUT_ERROR("Malformed FITS selection string");
	*status=EXIT_FAILURE;
	return(NULL);
     }
     retval=strndup(expr2,strlen(expr2)-1);
  } else {
	  // expr3 points at [. So character before that must be ]
	  if (*(expr3-1)!=']') {
		  SIMPUT_ERROR("Malformed FITS selection string");
		  *status=EXIT_FAILURE;
		  return(NULL);
	  }
	  retval=strndup(expr3,(size_t) (expr3-1));
  }

  CHECK_NULL_RET(retval,*status,"scanSpecFileName: Error allocating retval",NULL);

  headas_chat(5, "Extracted string: %s\n", retval);

  return retval;
}

/*
 * Read a MId spectrum from the indth cached spectrum from the row in row row
 */
SimputMIdpSpec *readCacheSpec(long ind, long row, char *fname, int *status)
{
  SimputMIdpSpec *spec;
  char buff[SIMPUT_MAXSTR];
  int anynull = 0;

  if ( SpecCache == NULL )
  {
    sprintf(buff, "Error reading cached spectrum: Cache not initialized\nFalling back to old method\n");
    SIMPUT_WARNING(buff);
    return NULL;
  }

  spec = newSimputMIdpSpec(status);
  CHECK_STATUS_RET(*status, NULL);

  spec->nentries = SpecCache->nbins[ind];

  spec->energy = (float *) malloc(spec->nentries * sizeof(float));
  CHECK_NULL_RET(spec->energy, *status, "Error allocating space for spectrum", NULL);
  spec->fluxdensity = (float *) malloc(spec->nentries * sizeof(float));
  CHECK_NULL_RET(spec->fluxdensity, *status, "Error allocating space for spectrum", NULL);

  // Read the data
  fits_read_col_flt(SpecCache->ext[ind], SpecCache->cenergy[ind], row, 1, SpecCache->nbins[ind], 0., spec->energy, &anynull, status);
  FITSERROR;

  fits_read_col_flt(SpecCache->ext[ind], SpecCache->cflux[ind], row, 1, SpecCache->nbins[ind], 0., spec->fluxdensity, &anynull, status);
  FITSERROR;

  if ( SpecCache->cname[ind] )
  {
    spec->name = (char *) calloc(SpecCache->namecol[ind]->namelen, sizeof(char));
    CHECK_NULL_RET(spec->name, *status, "Error allocating string buffer", NULL);
    fits_read_col_str(SpecCache->ext[ind], SpecCache->cname[ind], row,
	1, 1, "\0", &(spec->name), &anynull, status);
    FITSERROR;
  }

  for (long ii=0; ii<spec->nentries; ii++)
  {
    spec->energy[ii] *= SpecCache->fenergy[ind];
  }
  for (long ii=0; ii<spec->nentries; ii++)
  {
    spec->fluxdensity[ii] *= SpecCache->fflux[ind];
  }

  spec->fileref = (char *) malloc( (strlen(fname)+1) * sizeof(char));
  CHECK_NULL_RET(spec->fileref, *status, "Error allocating string buffer", NULL);
  strcpy(spec->fileref, fname);

  return spec;
}


long getSpecRow(char *expr, long ind)
{
  long row;
  char *pos;
  char msg[SIMPUT_MAXSTR];
  char *name;

  headas_chat(5, "Checking whether row or name was provided ... ");
  if ( (pos = strstr(expr, "#row==")) != NULL )
  {
    headas_chat(5, "row\n");
    pos += strlen("#row==");
    if ( sscanf(pos, "%ld", &row) == EOF )
    {
      sprintf(msg, "Error getting row number");
      SIMPUT_ERROR(msg);
      return -1;
    }
    headas_chat(5,"Row is %ld\n", row);
    return row;
  } else {
	  // let's see if a string like "NAME=='spec1'" is given
	  pos = strstr(expr, "NAME==");
	  // try also, if lowercase "name" is given if we did not match yet
	  if (pos == NULL){
		  pos = strstr(expr, "name==");
	  }

	  if ( pos != NULL ) {
		  char **ii;
		  headas_chat(5, "name\n");
		  name = pos + strlen("NAME==") + 1;
		  name[strlen(name)-2] = '\0';  // -2 for deleting the characters " '] "
		  headas_chat(5, "Name to search for: \"%s\"\n", name);

		  NamePtr = SpecCache->namecol[ind]->name;
		  ii = (char **) bsearch((const void *) &name,
				  (const void *) SpecCache->namecol[ind]->name,
				  (size_t) SpecCache->namecol[ind]->n,
				  sizeof(char *),
				  myStrCmp);
		  NamePtr = NULL;

		  row = SpecCache->namecol[ind]->row[ii - SpecCache->namecol[ind]->name];
		  headas_chat(5, "Position: %ld\n", row);
		  return row+1;
	  } else {
		  headas_chat(5, "nothing, returning -1\n");
		  return -1;
	  }
  }

  return -1;
}

long getNextSpecCache()
{
  if ( SpecCache->last < SPEC_MAX_CACHE - 1)
  {
    return SpecCache->last+1;
  }
  return 0;
}
