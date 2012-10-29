#include "common.h"


static void read_unit(fitsfile* const fptr, const int column, 
		      char* unit, int* const status)
{
  // Read the header keyword.
  char keyword[SIMPUT_MAXSTR], comment[SIMPUT_MAXSTR];
  sprintf(keyword, "TUNIT%d", column);
  fits_read_key(fptr, TSTRING, keyword, unit, comment, status);
  CHECK_STATUS_VOID(*status);
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
    CHECK_STATUS_BREAK(*status);

    // Split rootname into the file path and the file name.
    char* lastslash = strrchr(rootname, '/');
    if (NULL==lastslash) {
      cat->filepath=(char*)malloc(sizeof(char));
      CHECK_NULL_BREAK(cat->filepath, *status, 
		       "memory allocation for filepath failed");
      cat->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
      CHECK_NULL_BREAK(cat->filename, *status, 
		       "memory allocation for filename failed");
      strcpy(cat->filepath, "");
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
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // The file already exists => open it.
      fits_open_file(&cat->fptr, filename, mode, status);
      CHECK_STATUS_BREAK(*status);

    } else if (READWRITE==mode) {
      // The file does not exist, but it shall be created.
      fits_create_file(&cat->fptr, filename, status);
      CHECK_STATUS_BREAK(*status);

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
	if (EXIT_SUCCESS!=*status) break;
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
	CHECK_STATUS_BREAK(*status);
	
	// Write the necessary header keywords.
	fits_write_key(cat->fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
	fits_write_key(cat->fptr, TSTRING, "HDUCLAS1", "SRC_CAT", "", status);
	fits_write_key(cat->fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
	fits_write_key(cat->fptr, TSTRING, "RADESYS", "FK5", "", status);
	float equinox=2000.0;
	fits_update_key(cat->fptr, TFLOAT, "EQUINOX", &equinox, "", status);
	CHECK_STATUS_BREAK(*status);

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
    fits_get_colnum(cat->fptr, CASEINSEN, "RA", &cat->cra, status);
    fits_get_colnum(cat->fptr, CASEINSEN, "DEC", &cat->cdec, status);
    fits_get_colnum(cat->fptr, CASEINSEN, "E_MIN", &cat->ce_min, status);
    fits_get_colnum(cat->fptr, CASEINSEN, "E_MAX", &cat->ce_max, status);
    fits_get_colnum(cat->fptr, CASEINSEN, "FLUX", &cat->cflux, status);
    fits_get_colnum(cat->fptr, CASEINSEN, "SPECTRUM", &cat->cspectrum, status);
    CHECK_STATUS_BREAK(*status);
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
      SIMPUT_ERROR("unknown units in RA column");
      *status=EXIT_FAILURE;
      break;
    }

    char udec[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->cdec, udec, status);
    CHECK_STATUS_BREAK(*status);
    cat->fdec = unit_conversion_rad(udec);
    if (0.==cat->fdec) {
      SIMPUT_ERROR("unknown units in DEC column");
      *status=EXIT_FAILURE;
      break;
    }

    if (cat->cimgrota>0) {
      char uimgrota[SIMPUT_MAXSTR];
      read_unit(cat->fptr, cat->cimgrota, uimgrota, status);
      CHECK_STATUS_BREAK(*status);
      cat->fimgrota = unit_conversion_rad(uimgrota);
      if (0.==cat->fimgrota) {
	SIMPUT_ERROR("unknown units in IMGROTA column");
	*status=EXIT_FAILURE;
	break;
      }
    }

    char ue_min[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->ce_min, ue_min, status);
    CHECK_STATUS_BREAK(*status);
    cat->fe_min = unit_conversion_keV(ue_min);
    if (0.==cat->fe_min) {
      SIMPUT_ERROR("unknown units in E_MIN column");
      *status=EXIT_FAILURE;
      break;
    }

    char ue_max[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->ce_max, ue_max, status);
    CHECK_STATUS_BREAK(*status);
    cat->fe_max = unit_conversion_keV(ue_max);
    if (0.==cat->fe_max) {
      SIMPUT_ERROR("unknown units in E_MAX column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(cat->fptr, cat->cflux, uflux, status);
    CHECK_STATUS_BREAK(*status);
    cat->fflux = unit_conversion_ergpspcm2(uflux);
    if (0.==cat->fflux) {
      SIMPUT_ERROR("unknown units in FLUX column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries.
    fits_get_num_rows(cat->fptr, &cat->nentries, status);
    CHECK_STATUS_BREAK(*status);

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
    
    if (cat->csrc_name>0) {
      fits_read_col_str(cat->fptr, cat->csrc_name, row, 1, 1, 
			"", src_name, &anynul, status);
    } else {
      strcpy(src_name[0], "");
    }

    fits_read_col(cat->fptr, TDOUBLE, cat->cra, row, 1, 1, 
		  &ra, &ra, &anynul, status);
    ra *=cat->fra;  // Convert to [rad].
    fits_read_col(cat->fptr, TDOUBLE, cat->cdec, row, 1, 1, 
		  &dec, &dec, &anynul, status);
    dec*=cat->fdec; // Convert to [rad].

    if (cat->cimgrota>0) {
      fits_read_col(cat->fptr, TFLOAT, cat->cimgrota, row, 1, 1, 
		    &imgrota, &imgrota, &anynul, status);
      imgrota*=cat->fimgrota; // Convert to [rad].
    }
    if (cat->cimgscal>0) {
      fits_read_col(cat->fptr, TFLOAT, cat->cimgscal, row, 1, 1, 
		    &imgscal, &imgscal, &anynul, status);
    }
    
    fits_read_col(cat->fptr, TFLOAT, cat->ce_min, row, 1, 1, 
		  &e_min, &e_min, &anynul, status);
    e_min*=cat->fe_min; // Convert to [keV].
    fits_read_col(cat->fptr, TFLOAT, cat->ce_max, row, 1, 1, 
		  &e_max, &e_max, &anynul, status);
    e_max*=cat->fe_max; // Convert to [keV].
    fits_read_col(cat->fptr, TFLOAT, cat->cflux, row, 1, 1, 
		  &flux, &flux, &anynul, status);
    flux *=cat->fflux; // Convert to [erg/s/cm**2].
    if (flux > 1.e20) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "flux (%e erg/s/cm**2) exceeds maximum value, "
	      "therefore reset to 0", flux);
      SIMPUT_WARNING(msg);
      flux=0.;
    }

    fits_read_col(cat->fptr, TSTRING, cat->cspectrum, row, 1, 1, 
		  "", spectrum, &anynul, status);
    if (cat->cimage>0) {
      fits_read_col(cat->fptr, TSTRING, cat->cimage, row, 1, 1, 
		    "", image, &anynul, status);
    } else {
      strcpy(image[0], "");
    }

    if (cat->ctiming>0) {
      fits_read_col(cat->fptr, TSTRING, cat->ctiming, row, 1, 1, 
		    "", timing, &anynul, status);
    } else {
      strcpy(timing[0], "");
    }
    CHECK_STATUS_BREAK(*status);

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
  CHECK_STATUS_VOID(*status);
  cat->nentries++;

  fits_write_col(cat->fptr, TLONG, cat->csrc_id, cat->nentries, 1, 1, 
		 &src->src_id, status);
  fits_write_col(cat->fptr, TSTRING, cat->csrc_name, cat->nentries, 1, 1, 
		 &src->src_name, status);
  double ra = src->ra*180./M_PI;
  fits_write_col(cat->fptr, TDOUBLE, cat->cra, cat->nentries, 1, 1, 
		 &ra, status);
  double dec = src->dec*180./M_PI;
  fits_write_col(cat->fptr, TDOUBLE, cat->cdec, cat->nentries, 1, 1, 
		 &dec, status);
  float imgrota = src->imgrota*180./M_PI;
  fits_write_col(cat->fptr, TFLOAT, cat->cimgrota, cat->nentries, 1, 1, 
		 &imgrota, status);
  fits_write_col(cat->fptr, TFLOAT, cat->cimgscal, cat->nentries, 1, 1, 
		 &src->imgscal, status);
  fits_write_col(cat->fptr, TFLOAT, cat->ce_min, cat->nentries, 1, 1, 
		 &src->e_min, status);
  fits_write_col(cat->fptr, TFLOAT, cat->ce_max, cat->nentries, 1, 1, 
		 &src->e_max, status);
  fits_write_col(cat->fptr, TFLOAT, cat->cflux, cat->nentries, 1, 1, 
		 &src->eflux, status);
  fits_write_col(cat->fptr, TSTRING, cat->cspectrum, cat->nentries, 1, 1, 
		 &src->spectrum, status);
  fits_write_col(cat->fptr, TSTRING, cat->cimage, cat->nentries, 1, 1, 
		 &src->image, status);
  fits_write_col(cat->fptr, TSTRING, cat->ctiming, cat->nentries, 1, 1, 
		 &src->timing, status);
  CHECK_STATUS_VOID(*status);
}


void appendSimputSrcBlock(SimputCtlg* const cat,
			  SimputSrc** const src,
			  const long nsources,
			  int* const status)
{
  // Insert new rows.
  fits_insert_rows(cat->fptr, cat->nentries, nsources, status);
  CHECK_STATUS_VOID(*status);

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
      fits_write_col(cat->fptr, TSTRING, cat->cspectrum, first+ii, 1, 1, 
		     &(src[ii]->spectrum), status);
      fits_write_col(cat->fptr, TSTRING, cat->cimage, first+ii, 1, 1, 
		     &(src[ii]->image), status);
      fits_write_col(cat->fptr, TSTRING, cat->ctiming, first+ii, 1, 1, 
		     &(src[ii]->timing), status);
      CHECK_STATUS_BREAK(*status);
    }
    CHECK_STATUS_BREAK(*status);
    // END of loop over all sources.

    // Write the buffers to the file.
    fits_write_col(cat->fptr, TLONG, cat->csrc_id, first, 1, nsources, 
		   src_id, status);
    fits_write_col(cat->fptr, TDOUBLE, cat->cra, first, 1, nsources, 
		   ra, status);
    fits_write_col(cat->fptr, TDOUBLE, cat->cdec, first, 1, nsources, 
		   dec, status);
    fits_write_col(cat->fptr, TFLOAT, cat->cimgrota, first, 1, nsources, 
		   imgrota, status);
    fits_write_col(cat->fptr, TFLOAT, cat->cimgscal, first, 1, nsources, 
		   imgscal, status);
    fits_write_col(cat->fptr, TFLOAT, cat->ce_min, first, 1, nsources, 
		   e_min, status);
    fits_write_col(cat->fptr, TFLOAT, cat->ce_max, first, 1, nsources, 
		   e_max, status);
    fits_write_col(cat->fptr, TFLOAT, cat->cflux, first, 1, nsources, 
		   eflux, status);
    CHECK_STATUS_BREAK(*status);

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


SimputMIdpSpec* loadSimputMIdpSpec(const char* const filename,
				   int* const status)
{
  // String buffer.
  char* name[1]={NULL};

  SimputMIdpSpec* spec=newSimputMIdpSpec(status);
  CHECK_STATUS_RET(*status, spec);

  // Open the specified FITS file. The filename must uniquely identify
  // the spectrum contained in a binary table via the extended filename 
  // syntax. It must even specify the row, in which the spectrum is
  // contained. Therefore we do not have to care about the HDU or row
  // number.
  fitsfile* fptr=NULL;
  fits_open_table(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, spec);

  do { // Error handling loop.

    // Get the column names.
    int cenergy=0, cflux=0, cname=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
    CHECK_STATUS_BREAK(*status);
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
    float fenergy = unit_conversion_keV(uenergy);
    if (0.==fenergy) {
      SIMPUT_ERROR("unknown units in ENERGY column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(fptr, cflux, uflux, status);
    CHECK_STATUS_BREAK(*status);
    float fflux = unit_conversion_phpspcm2pkeV(uflux);
    if (0.==fflux) {
      SIMPUT_ERROR("unknown units in FLUX column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries in the 2 vector columns.
    int typecode;
    long nenergy, nflux, width;
    fits_get_coltype(fptr, cenergy, &typecode, &nenergy, &width, status);
    fits_get_coltype(fptr, cflux,   &typecode, &nflux,   &width, status);
    CHECK_STATUS_BREAK(*status);

    // If the columns are of variable-length data type, the returned repeat
    // value is 1. In that case we have to use another routine to get the
    // number of elements in a particular row.
    if (1==nenergy) {
      long offset;
      fits_read_descript(fptr, cenergy, 1, &nenergy, &offset, status);
      fits_read_descript(fptr, cflux  , 1, &nflux  , &offset, status);
      CHECK_STATUS_BREAK(*status);
    }

    // The number of energy bins and of flux entries must be identical.
    if (nenergy!=nflux) {
      SIMPUT_ERROR("number of energy and flux entries in spectrum is "
		   "not equivalent");
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
    spec->pflux =(float*)malloc(spec->nentries*sizeof(float));
    CHECK_NULL_BREAK(spec->pflux, *status, 
		     "memory allocation for spectrum failed");

    // Allocate memory for string buffer.
    name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(name[0], *status, 
		     "memory allocation for string buffer failed");

    // Read the data from the table.
    int anynul=0;
    fits_read_col(fptr, TFLOAT, cenergy, 1, 1, spec->nentries, 
		  NULL, spec->energy, &anynul, status);
    fits_read_col(fptr, TFLOAT, cflux, 1, 1, spec->nentries, 
		  NULL, spec->pflux, &anynul, status);

    if (cname>0) {
      fits_read_col(fptr, TSTRING, cname, 1, 1, 1, "", name, &anynul, status);
    } else { 
      strcpy(name[0], "");
    }

    CHECK_STATUS_BREAK(*status);

    // Multiply with unit scaling factor.
    long ii;
    for (ii=0; ii<spec->nentries; ii++) {
      spec->energy[ii] *= fenergy;
      spec->pflux[ii]  *= fflux;
    }

    // Copy the name (ID) of the spectrum from the string buffer
    // to the data structure.
    spec->name = (char*)malloc((strlen(name[0])+1)*sizeof(char));
    CHECK_NULL_BREAK(spec->name, *status, 
		     "memory allocation for name string failed");
    strcpy(spec->name, name[0]);

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

  do { // Error handling loop.

    // Check if the filename refers to a binary table extension
    // containing mission-independent spectra.
    int exttype=getExtType(cat, filename, status);
    CHECK_STATUS_BREAK(*status);

    if (EXTTYPE_MIDPSPEC!=exttype) {
      // Only mission-independent spectra can be pre-loaded into 
      // the cache.
      break;
    }

    // Check if the source catalog contains a spectrum buffer.
    if (NULL==cat->midpspecbuff) {
      cat->midpspecbuff=newSimputMIdpSpecBuffer(status);
      CHECK_STATUS_BREAK(*status);
    }

    // Convert the void* pointer to the spectrum buffer into the right
    // format.
    struct SimputMIdpSpecBuffer* sb=
      (struct SimputMIdpSpecBuffer*)cat->midpspecbuff;

    // In case there are no spectra available at all, allocate 
    // memory for the array (storage for spectra).
    if (NULL==sb->spectra) {
      sb->spectra=
	(SimputMIdpSpec**)malloc(MAXMIDPSPEC*sizeof(SimputMIdpSpec*));
      CHECK_NULL_BREAK(sb->spectra, *status, 
		       "memory allocation for spectra failed");
    }

    // Open the specified FITS file. The filename must uniquely identify
    // the extension containing the spectra via the extended filename 
    // syntax.
    fits_open_table(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);

    // Determine the number of columns in the table.
    long nrows;
    fits_get_num_rows(fptr, &nrows, status);
    CHECK_STATUS_BREAK(*status);

    // Check if as many spectra can be stored in the cache.
    if (sb->nspectra+nrows>MAXMIDPSPEC) {
      *status=EXIT_FAILURE;
      SIMPUT_ERROR("cache too small to store all spectra");
      break;
    }

    // Get the column names.
    int cenergy=0, cflux=0, cname=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
    CHECK_STATUS_BREAK(*status);
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
      SIMPUT_ERROR("unknown units in ENERGY column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(fptr, cflux, uflux, status);
    CHECK_STATUS_BREAK(*status);
    float fflux=unit_conversion_phpspcm2pkeV(uflux);
    if (0.==fflux) {
      SIMPUT_ERROR("unknown units in FLUX column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries in the 2 vector columns.
    int typecode;
    long nenergy, nflux, width;
    fits_get_coltype(fptr, cenergy, &typecode, &nenergy, &width, status);
    fits_get_coltype(fptr, cflux,   &typecode, &nflux,   &width, status);
    CHECK_STATUS_BREAK(*status);

    // If the columns are of variable-length data type, the returned repeat
    // value is 1. In that case we have to use another routine to get the
    // number of elements in a particular row.
    if (1==nenergy) {
      long offset;
      fits_read_descript(fptr, cenergy, 1, &nenergy, &offset, status);
      fits_read_descript(fptr, cflux  , 1, &nflux  , &offset, status);
      CHECK_STATUS_BREAK(*status);
    }

    // The number of energy bins and of flux entries must be identical.
    if (nenergy!=nflux) {
      SIMPUT_ERROR("number of energy and flux entries in spectrum is "
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
      spec->pflux =(float*)malloc(spec->nentries*sizeof(float));
      CHECK_NULL_BREAK(spec->pflux, *status, 
		       "memory allocation for spectrum failed");

      // Read the data from the table.
      int anynul=0;
      fits_read_col(fptr, TFLOAT, cenergy, jj+1, 1, spec->nentries, 
		    NULL, spec->energy, &anynul, status);
      fits_read_col(fptr, TFLOAT, cflux, jj+1, 1, spec->nentries, 
		    NULL, spec->pflux, &anynul, status);
      if (cname>0) {
	fits_read_col(fptr, TSTRING, cname, jj+1, 1, 1, "", name, 
		      &anynul, status);
      } else { 
	strcpy(name[0], "");
      }
      CHECK_STATUS_BREAK(*status);

      // Multiply with unit scaling factor.
      long ii;
      for (ii=0; ii<spec->nentries; ii++) {
	spec->energy[ii]*=fenergy;
	spec->pflux[ii] *=fflux;
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

      // Add the spectrum to the cache.
      sb->spectra[sb->nspectra++]=spec;
    }
    CHECK_STATUS_BREAK(*status);
    // END of reading all spectra.

  } while(0); // END of error handling loop.
  
  // Release allocated memory.
  if (NULL!=name[0]) free(name[0]);

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
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      CHECK_STATUS_BREAK(*status);
    } else {
      // If no, create a new file. 
      fits_create_file(&fptr, filename, status);
      CHECK_STATUS_BREAK(*status);
    }
    // END of check, whether the specified file exists.


    // Try to move to the specified extension.
    int cenergy=0, cflux=0, cname=0;
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

      strcpy(ttype[1], "FLUX");
      sprintf(tform[1], "1PE");
      strcpy(tunit[1], "photon/s/cm**2/keV");

      strcpy(ttype[2], "NAME");
      strcpy(tform[2], "48A");
      strcpy(tunit[2], "");

      // Create the table.
      fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, extname, status);
      CHECK_STATUS_BREAK(*status);

      // Write header keywords.
      fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS1", "SPECTRUM", "", status);
      fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
      fits_write_key(fptr, TINT, "EXTVER", &extver, "", status);
      CHECK_STATUS_BREAK(*status);

      // The new table contains now data up to now.
      nrows=0;

    } else {
      // The extension already exists.
      // Determine the number of contained rows.
      fits_get_num_rows(fptr, &nrows, status);
      CHECK_STATUS_BREAK(*status);
    }
    // END of check, whether the specified extension exists.


    // Determine the column numbers.
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, status);
    fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
    CHECK_STATUS_BREAK(*status);
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
	  SIMPUT_ERROR("NAME value of spectrum contains more than 48 characters");
	  *status=EXIT_FAILURE;
	  break;
	}

	// Check if the NAME column is present.
	if (0==cname) {
	  SIMPUT_ERROR("spectrum extension does not contain a NAME column");
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
	  if (0==strcmp(name[0], spec->name)) {
	    SIMPUT_ERROR("name in spectrum data structure is not unique");
	    *status=EXIT_FAILURE;
	    break;
	  }
	}
	CHECK_STATUS_BREAK(*status);
      }
    }

    // Create a new row in the table and store the data of the spectrum in it.
    fits_insert_rows(fptr, nrows++, 1, status);
    CHECK_STATUS_BREAK(*status);
    fits_write_col(fptr, TFLOAT, cenergy, nrows, 1, spec->nentries, 
		   spec->energy, status);
    fits_write_col(fptr, TFLOAT, cflux, nrows, 1, spec->nentries, 
		   spec->pflux, status);
    if ((cname>0) && (NULL!=spec->name)) {
      fits_write_col(fptr, TSTRING, cname, nrows, 1, 1, &spec->name, status);
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
    CHECK_STATUS_BREAK(*status);

    // Get an empty SimputLC data structure.
    lc=newSimputLC(status);
    CHECK_STATUS_BREAK(*status);

    // Get the column names.
    int ctime=0, cphase=0, cflux=0, cspectrum=0, cimage=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
    CHECK_STATUS_BREAK(*status);
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

    // Print warnings.
    if (cspectrum>0) {
      SIMPUT_WARNING("the use of spectral variation via light curves "
		     "is not fully implemented");
    }
    if (cspectrum>0) {
      SIMPUT_WARNING("the use of spatial variation via light curves "
		     "is not fully implemented");
    }

    // Check, whether there is either a TIME or a PHASE column (but not both).
    if ((0==ctime)&&(0==cphase)) {
      SIMPUT_ERROR("table extension contains neither TIME nor PHASE column");
      *status=EXIT_FAILURE;
      return(lc);
    } else if ((ctime>0)&&(cphase>0)) {
      SIMPUT_ERROR("table extension contains both TIME and PHASE column");
      *status=EXIT_FAILURE;
      return(lc);
    }

    // Determine the unit conversion factors.
    float ftime=0.;
    if (ctime>0) {
      char utime[SIMPUT_MAXSTR];
      read_unit(fptr, ctime, utime, status);
      CHECK_STATUS_BREAK(*status);
      ftime = unit_conversion_s(utime);
      if (0.==ftime) {
	SIMPUT_ERROR("unknown units in TIME column");
	*status=EXIT_FAILURE;
	break;
      }
    }
    // END of determine unit conversion factors.


    // Read the header keywords.
    char comment[SIMPUT_MAXSTR];
    // Required keywords.
    fits_read_key(fptr, TDOUBLE, "MJDREF",   &lc->mjdref,   comment, status);
    fits_read_key(fptr, TDOUBLE, "TIMEZERO", &lc->timezero, comment, status);
    if (cphase>0) {
      // Only for periodic light curves.
      fits_read_key(fptr, TDOUBLE, "PHASE0", &lc->phase0, comment, status);
      fits_read_key(fptr, TDOUBLE, "PERIOD", &lc->period, comment, status);
    } else {
      lc->phase0 = 0.;
      lc->period = 0.;
    }
    CHECK_STATUS_BREAK(*status);
    // Optional keywords.
    opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_read_key(fptr, TFLOAT, "FLUXSCAL", &lc->fluxscal, comment, &opt_status);
    if (EXIT_SUCCESS!=opt_status) {
      // FLUXSCAL is not given in the FITS header. We therefore assume
      // that it has a value of 1.
      lc->fluxscal = 1.;
      opt_status=EXIT_SUCCESS;
    }
    fits_clear_errmark();


    // Determine the number of rows in the table.
    fits_get_num_rows(fptr, &lc->nentries, status);
    CHECK_STATUS_BREAK(*status);

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
    int anynul=0;

    // TIME
    if (ctime>0) {
      fits_read_col(fptr, TDOUBLE, ctime, 1, 1, lc->nentries, 
		    0, lc->time, &anynul, status);
      CHECK_STATUS_BREAK(*status);
      // Multiply with unit scaling factor.
      long row;
      for (row=0; row<lc->nentries; row++) {
	lc->time[row]*=ftime;
      }
    }

    // PHASE
    if (cphase>0) {
      fits_read_col(fptr, TDOUBLE, cphase, 1, 1, lc->nentries, 
		    0, lc->phase, &anynul, status);
      CHECK_STATUS_BREAK(*status);
    }

    // FLUX
    fits_read_col(fptr, TFLOAT, cflux, 1, 1, lc->nentries, 
		  0, lc->flux, &anynul, status);
    CHECK_STATUS_BREAK(*status);

    // SPECTRUM
    if (cspectrum>0) {
      long row;
      for (row=0; row<lc->nentries; row++) {
	fits_read_col(fptr, TSTRING, cspectrum, row+1, 1, 1, 
		      "", spectrum, &anynul, status);
	CHECK_STATUS_BREAK(*status);
	lc->spectrum[row]=
	  (char*)malloc((strlen(spectrum[0])+1)*sizeof(char));
	CHECK_NULL_BREAK(lc->spectrum[row], *status,
			 "memory allocation for spectrum string failed");
	strcpy(lc->spectrum[row], spectrum[0]);
      }      
      CHECK_STATUS_BREAK(*status);
    }

    // IMAGE
    if (cimage>0) {
      long row;
      for (row=0; row<lc->nentries; row++) {
	fits_read_col(fptr, TSTRING, cimage, row+1, 1, 1, 
		      "", image, &anynul, status);
	CHECK_STATUS_BREAK(*status);
	lc->image[row]=
	  (char*)malloc((strlen(image[0])+1)*sizeof(char));
	CHECK_NULL_BREAK(lc->image[row], *status,
			 "memory allocation for image string failed");
	strcpy(lc->image[row], image[0]);
      }      
      CHECK_STATUS_BREAK(*status);
    }
    // END of reading the data from the FITS table.

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
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      CHECK_STATUS_BREAK(*status);
      
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
      CHECK_STATUS_BREAK(*status);
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
    CHECK_STATUS_BREAK(*status);

    // Write header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "LIGHTCURVE", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);
    fits_write_key(fptr, TDOUBLE, "MJDREF", &lc->mjdref, "", status);
    fits_write_key(fptr, TDOUBLE, "TIMEZERO", &lc->timezero, "", status);
    fits_write_key(fptr, TFLOAT,  "FLUXSCAL", &lc->fluxscal, "", status);
    int periodic=0;
    if (cphase>0) {
      // Only for periodic light curves.
      periodic=1;
      fits_write_key(fptr, TDOUBLE, "PHASE0", &lc->phase0, "", status);
      fits_write_key(fptr, TDOUBLE, "PERIOD", &lc->period, "", status);
    }
    fits_write_key(fptr, TINT,  "PERIODIC", &periodic, "", status);
    CHECK_STATUS_BREAK(*status);

    // Create new rows in the table and store the data of the spectrum in it.
    fits_insert_rows(fptr, 0, lc->nentries, status);
    CHECK_STATUS_BREAK(*status);
    
    if (ctime>0) {
      fits_write_col(fptr, TDOUBLE, ctime, 1, 1, lc->nentries, 
		     lc->time, status);
      CHECK_STATUS_BREAK(*status);
    } else {
      fits_write_col(fptr, TDOUBLE, cphase, 1, 1, lc->nentries, 
		     lc->phase, status);
      CHECK_STATUS_BREAK(*status);
    }
    fits_write_col(fptr, TFLOAT, cflux, 1, 1, lc->nentries, 
		   lc->flux, status);
    CHECK_STATUS_BREAK(*status);
    if (cspectrum>0) {
      long row;
      for (row=0; row<lc->nentries; row++) {
	fits_write_col(fptr, TSTRING, cspectrum, row+1, 1, 1, 
		       &lc->spectrum[row], status);
	CHECK_STATUS_BREAK(*status);
      }
      CHECK_STATUS_BREAK(*status);
    }
    if (cimage>0) {
      long row;
      for (row=0; row<lc->nentries; row++) {
	fits_write_col(fptr, TSTRING, cimage, row+1, 1, 1, 
		       &lc->image[row], status);
	CHECK_STATUS_BREAK(*status);
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
  CHECK_STATUS_RET(*status, psd);

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
    CHECK_STATUS_BREAK(*status);

    fits_get_colnum(fptr, CASEINSEN, "POWER", &cpower, status);
    CHECK_STATUS_BREAK(*status);

    // Determine the unit conversion factors.
    float ffrequency=0.;
    char ufrequency[SIMPUT_MAXSTR];
    read_unit(fptr, cfrequency, ufrequency, status);
    CHECK_STATUS_BREAK(*status);
    ffrequency = unit_conversion_Hz(ufrequency);
    if (0.==ffrequency) {
      SIMPUT_ERROR("unknown units in FREQUENCY column");
      *status=EXIT_FAILURE;
      break;
    }

    float fpower=0.;
    char upower[SIMPUT_MAXSTR];
    read_unit(fptr, cpower, upower, status);
    CHECK_STATUS_BREAK(*status);
    fpower = unit_conversion_s(upower);
    if (0.==fpower) {
      SIMPUT_ERROR("unknown units in POWER column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.


    // Determine the number of rows in the table.
    fits_get_num_rows(fptr, &psd->nentries, status);
    CHECK_STATUS_BREAK(*status);

    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "PSD '%s' contains %ld data points\n", 
	    filename, psd->nentries);
    SIMPUT_INFO(msg);

    // Allocate memory for the arrays.
    psd->frequency= (float*)malloc(psd->nentries*sizeof(float));
    CHECK_NULL_BREAK(psd->frequency, *status, 
		     "memory allocation for PSD failed");
    psd->power    = (float*)malloc(psd->nentries*sizeof(float));
    CHECK_NULL_BREAK(psd->power, *status, 
		     "memory allocation for PSD failed");

    // Read the data from the table.
    int anynul=0;
    // FREQUENC
    fits_read_col(fptr, TFLOAT, cfrequency, 1, 1, psd->nentries, 
		  0, psd->frequency, &anynul, status);
    CHECK_STATUS_BREAK(*status);
    // POWER
    fits_read_col(fptr, TFLOAT, cpower, 1, 1, psd->nentries, 
		  0, psd->power, &anynul, status);
    CHECK_STATUS_BREAK(*status);
    // END of reading the data from the FITS table.

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
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // If yes, open it.
      fits_open_file(&fptr, filename, READWRITE, status);
      CHECK_STATUS_BREAK(*status);

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
      CHECK_STATUS_BREAK(*status);
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
    CHECK_STATUS_BREAK(*status);

    // Write header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "POWSPEC", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);

    // Create new rows in the table and store the data of the power 
    // spectrum in it.
    fits_insert_rows(fptr, 0, psd->nentries, status);
    CHECK_STATUS_BREAK(*status);

    fits_write_col(fptr, TFLOAT, cfreq, 1, 1, psd->nentries,
       psd->frequency, status);
    CHECK_STATUS_BREAK(*status);
    fits_write_col(fptr, TFLOAT, cpower, 1, 1, psd->nentries,
       psd->power, status);
    CHECK_STATUS_BREAK(*status);

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

  // Get an empty SimputImg data structure.
  SimputImg* img=newSimputImg(status);
  CHECK_STATUS_RET(*status, img);

  // Open the specified FITS file. The filename must uniquely identify
  // the light curve contained in a binary table via the extended filename 
  // syntax. Therefore we do not have to care about the HDU number.
  fitsfile* fptr=NULL;
  fits_open_image(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, img);

  do { // Error handling loop.

    // Read the WCS header keywords.
    // Read the entire header into the string buffer.
    int nkeys;
    fits_hdr2str(fptr, 1, NULL, 0, &headerstr, &nkeys, status);
    CHECK_STATUS_BREAK(*status);
    // Parse the header string and store the data in the wcsprm data
    // structure.
    int nreject, nwcs;
    if (0!=wcspih(headerstr, nkeys, 0, 3, &nreject, &nwcs, &img->wcs)) {
      SIMPUT_ERROR("parsing of WCS header failed");
      *status = EXIT_FAILURE;
      break;
    }
    if (nreject>0) {
      SIMPUT_ERROR("parsing of WCS header failed");
      *status = EXIT_FAILURE;
      break;
    }

    // Determine the image dimensions.
    int naxis;
    fits_get_img_dim(fptr, &naxis, status);
    CHECK_STATUS_BREAK(*status);
    if (2!=naxis) {
      SIMPUT_ERROR("specified FITS HDU does not contain a 2-dimensional image");
      *status=EXIT_FAILURE;
      break;
    }
    long naxes[2];
    fits_get_img_size(fptr, naxis, naxes, status);
    CHECK_STATUS_BREAK(*status);
    img->naxis1 = naxes[0];
    img->naxis2 = naxes[1];

    // Allocate memory for the image.
    img->dist = (double**)malloc(img->naxis1*sizeof(double*));
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
    double null_value=0.;
    long fpixel[2] = {1, 1};   // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1).
    long lpixel[2] = {img->naxis1, img->naxis2}; // Upper right corner.
    long inc[2]    = {1, 1};
    fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value, 
		     image1d, &anynul, status);
    CHECK_STATUS_BREAK(*status);
    
    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure and generate a probability distribution function,
    // i.e., sum up the pixels.
    double sum=0.;
    for(ii=0; ii<img->naxis1; ii++) {
      long jj;
      for(jj=0; jj<img->naxis2; jj++) {
	sum += image1d[ii+ img->naxis1*jj];
	img->dist[ii][jj] = sum;
      }
    }

    // Read the optional FLUXSCAL header keyword.
    char comment[SIMPUT_MAXSTR];
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_read_key(fptr, TFLOAT, "FLUXSCAL", &img->fluxscal, 
		  comment, &opt_status);
    if (EXIT_SUCCESS!=opt_status) {
      // FLUXSCAL is not given in the FITS header. Therefore it is 
      // set to the default value of 1.
      img->fluxscal = 1.;
    }
    fits_clear_errmark();

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
      CHECK_STATUS_BREAK(*status);
      
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
      CHECK_STATUS_BREAK(*status);
    }
    // END of check, whether the specified file exists.


    // Allocate memory for the 1-dimensional image buffer 
    // (required for output to FITS file).
    image1d = (double*)malloc(img->naxis1*img->naxis2*sizeof(double));
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
    CHECK_STATUS_BREAK(*status);
    // The image has been appended at the end if the FITS file.

    // Write header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC/SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "IMAGE", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.1.0", "", status);
    fits_write_key(fptr, TSTRING, "EXTNAME", extname, "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);
    fits_write_key(fptr, TFLOAT,  "FLUXSCAL", &img->fluxscal, "", status);
    CHECK_STATUS_BREAK(*status);

    // Write WCS header keywords.
    int nkeyrec;
    if (0!=wcshdo(0, img->wcs, &nkeyrec, &headerstr)) {
      SIMPUT_ERROR("construction of WCS header failed");
      *status=EXIT_FAILURE;
      break;
    }
    char* strptr=headerstr;
    while (strlen(strptr)>0) {
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
    CHECK_STATUS_BREAK(*status);

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
  CHECK_STATUS_RET(*status, phl);

  // Determine the total number of photons.
  fits_get_num_rows(phl->fptr, &phl->nphs, status);
  CHECK_STATUS_RET(*status, phl);

  // Determine the reference RA and DEC values from the
  // header keywords.
  double refra=0., refdec=0.;
  char comment[SIMPUT_MAXSTR];
  fits_read_key(phl->fptr, TDOUBLE, "REFRA", &refra, comment, status);
  fits_read_key(phl->fptr, TDOUBLE, "REFDEC", &refdec, comment, status);
  CHECK_STATUS_RET(*status, phl);
  if ((0.!=refra)||(0.!=refdec)) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("in current implementation photon lists must have "
		 "REFRA=0.0 and REFDEC=0.0");
    return(phl);
  }

  // Determine the column numbers.
  fits_get_colnum(phl->fptr, CASEINSEN, "RA", &phl->cra, status);
  fits_get_colnum(phl->fptr, CASEINSEN, "DEC", &phl->cdec, status);
  fits_get_colnum(phl->fptr, CASEINSEN, "ENERGY", &phl->cenergy, status);
  CHECK_STATUS_RET(*status, phl);

  // Determine the unit conversion factors.
  char ura[SIMPUT_MAXSTR];
  read_unit(phl->fptr, phl->cra, ura, status);
  CHECK_STATUS_RET(*status, phl);
  phl->fra=unit_conversion_rad(ura);
  if (0.==phl->fra) {
    SIMPUT_ERROR("unknown units in RA column");
    *status=EXIT_FAILURE;
    return(phl);
  }

  char udec[SIMPUT_MAXSTR];
  read_unit(phl->fptr, phl->cdec, udec, status);
  CHECK_STATUS_RET(*status, phl);
  phl->fdec=unit_conversion_rad(udec);
  if (0.==phl->fdec) {
    SIMPUT_ERROR("unknown units in DEC column");
    *status=EXIT_FAILURE;
    return(phl);
  }

  char uenergy[SIMPUT_MAXSTR];
  read_unit(phl->fptr, phl->cenergy, uenergy, status);
  CHECK_STATUS_RET(*status, phl);
  phl->fenergy=unit_conversion_keV(uenergy);
  if (0.==phl->fenergy) {
    SIMPUT_ERROR("unknown units in ENERGY column");
    *status=EXIT_FAILURE;
    return(phl);
  }

  return(phl);
}


int getExtType(SimputCtlg* const cat, 
	       const char* const filename, 
	       int* const status)
{
  int type=EXTTYPE_NONE;

  // Check if there is any reference at all.
  if (0==strlen(filename)) {
    return(EXTTYPE_NONE);
  }

  // Keep an internal cache of extension types in order to avoid 
  // continuous re-opening.

  // Maximum number of extensions in the cache.
  const int maxhdus=100000; 

  // Check if the source catalog contains an extension type buffer.
  if (NULL==cat->extbuff) {
    cat->extbuff=newSimputExttypeBuffer(status);
    CHECK_STATUS_RET(*status, type);
  }

  // Convert the void* pointer to the extension type buffer 
  // into the right format.
  struct SimputExttypeBuffer* eb=
    (struct SimputExttypeBuffer*)cat->extbuff;

  // In case there are no extensions available at all, allocate 
  // memory for the array (storage for extension types).
  if (NULL==eb->hdus) {
    eb->hdus=(int*)malloc(maxhdus*sizeof(int));
    CHECK_NULL_RET(eb->hdus, *status, 
		   "memory allocation for extension types failed", type);
    eb->filenames=(char**)malloc(maxhdus*sizeof(char*));
    CHECK_NULL_RET(eb->filenames, *status, 
		   "memory allocation for extension types failed", type);
    long ii;
    for (ii=0; ii<maxhdus; ii++) {
      eb->hdus[ii]=EXTTYPE_NONE;
      eb->filenames[ii]=NULL;
    }
  }

  // Search if the required extension is available in the storage.
  long ii;
  for (ii=0; ii<eb->nhdus; ii++) {
    // Check if the extension is equivalent to the required one.
    if (0==strcmp(eb->filenames[ii], filename)) {
      // If yes, return the extension type.
      return(eb->hdus[ii]);
    }
  }


  // The extension is not contained in the cache. Therefore 
  // we have to open it an check the header keywords.
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, type);

  // Read the HDUCLAS1 and HDUCLAS2 header keywords.
  char comment[SIMPUT_MAXSTR];
  char hduclas1[SIMPUT_MAXSTR];
  char hduclas2[SIMPUT_MAXSTR];
  fits_read_key(fptr, TSTRING, "HDUCLAS1", &hduclas1, comment, status);
  // (Don't do an error checking here! Otherwise the file
  //  might not be closed after an error occurred with reading
  //  the header keyword.)

  // Read optional header keyword (is not used in SIMPUT
  // version >= 1.1.0).
  int opt_status=EXIT_SUCCESS;
  fits_write_errmark();
  fits_read_key(fptr, TSTRING, "HDUCLAS2", &hduclas2, comment, &opt_status);
  fits_clear_errmark();
  if (opt_status!=EXIT_SUCCESS) {
    strcpy(hduclas2, "");
    opt_status=EXIT_SUCCESS;
  }
  
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, type);


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
    sprintf(msg, "extension type '%s' not supported", hduclas1);
    SIMPUT_ERROR(msg);
    *status=EXIT_FAILURE;
    return(type);
  }

  // Store the extension type in the internal cache.
  // Determine the storage position.
  if (eb->nhdus<maxhdus) {
    eb->chdu=eb->nhdus;
    eb->nhdus++;
  } else {
    eb->chdu++;
    if (eb->chdu>=maxhdus) {
      eb->chdu=0;
    }
    free(eb->filenames[eb->chdu]);
  }
  eb->hdus[eb->chdu]=type;
  eb->filenames[eb->chdu]=
    (char*)malloc((strlen(filename)+1)*sizeof(char));
  CHECK_NULL_RET(eb->filenames[eb->chdu], *status, 
		 "memory allocation for file reference failed", 
		 eb->hdus[eb->chdu]);
  strcpy(eb->filenames[eb->chdu], filename);


  return(type);
}


