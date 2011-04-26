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


/////////////////////////////////////////////////////////////////
// Macros.
/////////////////////////////////////////////////////////////////


/** Common string length. */
#define SIMPUT_MAXSTR (1025)


#define SIMPUT_ERROR(msg) \
  fprintf(stderr, "Error in %s: %s!\n", __func__, msg)

#define CHECK_STATUS_RET(a,b) \
  if (EXIT_SUCCESS!=a) return(b)

#define CHECK_STATUS_VOID(a) \
  if (EXIT_SUCCESS!=a) return

#define CHECK_STATUS_BREAK(a) \
  if (EXIT_SUCCESS!=a) break;

#define CHECK_NULL_RET(a,status,msg,ret) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return(ret);\
  }

#define CHECK_NULL_VOID(a,status,msg) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return;\
  }

#define CHECK_NULL_BREAK(a,status,msg) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    break;\
  }

/** Macro returning the maximum of 2 values. */
#define MAX(a, b) ( (a)>(b) ? (a) : (b) )

/** Macro returning the minimum of 2 values. */
#define MIN(a, b) ( (a)<(b) ? (a) : (b) )


/////////////////////////////////////////////////////////////////
// Structures.
/////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////
// Static Variables.
/////////////////////////////////////////////////////////////////


/** Instrument ARF. */
static struct ARF* static_arf=NULL;

/** Random number generator. */
static double(*static_rndgen)(void)=NULL;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Check if the HDU referred to by filename is a binary table
    extension. */
static int check_if_btbl(const char* const filename, int* const status);

/** Read the TUNITn keyword from the header of the current FITS HDU,
    where n denotes the specfied column. */
static void read_unit(fitsfile* const fptr, const int column, 
		      char* unit, int* const status);

/** Convert a string into lower case letters. The string has to be
    terminated by a '\0' mark. */
static void strtolower(char* const string);

/** Determine the factor required to convert the specified unit into
    [rad]. If the conversion is not possible or implemented, the
    function return value is 0. */
static float unit_conversion_rad(const char* const unit);

/** Determine the factor required to convert the specified unit into
    [keV]. If the conversion is not possible or implemented, the
    function return value is 0. */
static float unit_conversion_keV(const char* const unit);

/** Determine the factor required to convert the specified unit into
    [erg/s/cm**2]. If the conversion is not possible or implemented,
    the function return value is 0. */
static float unit_conversion_ergpspcm2(const char* const unit);

/** Determine the factor required to convert the specified unit into
    [ph/s/cm**2/keV]. If the conversion is not possible or
    implemented, the function return value is 0. */
static float unit_conversion_phpspcm2pkeV(const char* const unit);

/** Return the requested spectrum. Keeps a certain number of spectra
    in an internal storage. If the requested spectrum is not located
    in the internal storage, it is loaded from the reference given in
    the source catalog. */
static SimputMissionIndepSpec* 
returnSimputMissionIndepSpec(const SimputSourceEntry* const src,
			     int* const status);

/** Determine a random photon energy according to the specified
    spectral distribution. */
static float getRndPhotonEnergy(const SimputMissionIndepSpec* const spec,
				      int* const status);

/** Determine the source flux in [erg/s/cm**2] within a certain energy
    band for the particular spectrum. */
static float getEbandFlux(const SimputSourceEntry* const src,
			  const float emin, const float emax,
			  int* const status);

/** Determine the photon rate in [photon/s] within a certain energy
    band for the particular spectrum. */
static float getEbandRate(const SimputSourceEntry* const src,
			  const float emin, const float emax,
			  int* const status);


/////////////////////////////////////////////////////////////////
// Function Definitions.
/////////////////////////////////////////////////////////////////


SimputSourceEntry* getSimputSourceEntry(int* const status)
{
  SimputSourceEntry* entry=(SimputSourceEntry*)malloc(sizeof(SimputSourceEntry));
  CHECK_NULL_RET(entry, *status, 
		 "memory allocation for SimputSourceEntry failed", entry);

  // Initialize elements.
  entry->src_id  =0;
  entry->src_name=NULL;
  entry->ra      =0.;
  entry->dec     =0.;
  entry->imgrota =0.;
  entry->imgscal =1.;
  entry->e_min   =0.;
  entry->e_max   =0.;
  entry->flux    =0.;
  entry->spectrum=NULL;
  entry->image   =NULL;
  entry->lightcur=NULL;
  entry->filename=NULL;
  entry->filepath=NULL;

  return(entry);
}


SimputSourceEntry* getSimputSourceEntryV(const int src_id, 
					 const char* const src_name,
					 const double ra,
					 const double dec,
					 const float imgrota,
					 const float imgscal,
					 const float e_min,
					 const float e_max,
					 const float flux,
					 const char* const spectrum,
					 const char* const image,
					 const char* const lightcur,
					 int* const status)
{
  SimputSourceEntry* entry=getSimputSourceEntry(status);
  CHECK_STATUS_RET(*status, entry);

  // Initialize with the given values.
  entry->src_id = src_id;

  entry->src_name=(char*)malloc((strlen(src_name)+1)*sizeof(char));
  CHECK_NULL_RET(entry->src_name, *status,
		 "memory allocation for source name failed", entry);
  strcpy(entry->src_name, src_name);

  entry->ra      = ra;
  entry->dec     = dec;
  entry->imgrota = imgrota;
  entry->imgscal = imgscal;
  entry->e_min   = e_min;
  entry->e_max   = e_max;
  entry->flux    = flux;

  entry->spectrum=(char*)malloc((strlen(spectrum)+1)*sizeof(char));
  CHECK_NULL_RET(entry->spectrum, *status,
		 "memory allocation for source name failed", entry);
  strcpy(entry->spectrum, spectrum);

  entry->image  =(char*)malloc((strlen(image)+1)*sizeof(char));
  CHECK_NULL_RET(entry->image, *status,
		 "memory allocation for source name failed", entry);
  strcpy(entry->image, image);

  entry->lightcur=(char*)malloc((strlen(lightcur)+1)*sizeof(char));
  CHECK_NULL_RET(entry->lightcur, *status,
		 "memory allocation for source name failed", entry);
  strcpy(entry->lightcur, lightcur);
  

  return(entry);
}


void freeSimputSourceEntry(SimputSourceEntry** const entry)
{
  if (NULL!=*entry) {
    if (NULL!=(*entry)->src_name) {
      free((*entry)->src_name);
    }
    if (NULL!=(*entry)->spectrum) {
      free((*entry)->spectrum);
    }
    if (NULL!=(*entry)->image) {
      free((*entry)->image);
    }
    if (NULL!=(*entry)->lightcur) {
      free((*entry)->lightcur);
    }
    free(*entry);
    *entry=NULL;
  }
}


SimputSourceCatalog* getSimputSourceCatalog(int* const status)
{
  SimputSourceCatalog* catalog=(SimputSourceCatalog*)malloc(sizeof(SimputSourceCatalog));
  CHECK_NULL_RET(catalog, *status, 
		 "memory allocation for SimputSourceCatalog failed", catalog);

  // Initialize elements.
  catalog->nentries=0;
  catalog->entries =NULL;
  catalog->filepath=NULL;
  catalog->filename=NULL;

  return(catalog);
}


void freeSimputSourceCatalog(SimputSourceCatalog** const catalog)
{
  if (NULL!=*catalog) {
    if ((*catalog)->nentries>0) {
      int ii;
      for (ii=0; ii<(*catalog)->nentries; ii++) {
	if (NULL!=(*catalog)->entries[ii]) {
	  freeSimputSourceEntry(&((*catalog)->entries[ii]));
	}
      }
    }
    if (NULL!=(*catalog)->entries) {
      free((*catalog)->entries);
    }
    if (NULL!=(*catalog)->filepath) {
      free((*catalog)->filepath);
    }
    if (NULL!=(*catalog)->filename) {
      free((*catalog)->filename);
    }
    free(*catalog);
    *catalog=NULL;
  }
}


SimputSourceCatalog* loadSimputSourceCatalog(const char* const filename,
					     int* const status)
{
  SimputSourceCatalog* catalog = getSimputSourceCatalog(status);
  CHECK_STATUS_RET(*status, catalog);

  // Store the filename and filepath of the FITS file containing
  // the source catalog.
  char cfilename[SIMPUT_MAXSTR];
  char rootname[SIMPUT_MAXSTR];
  // Make a local copy of the filename variable in order to avoid
  // compiler warnings due to discarded const qualifier at the 
  // subsequent function call.
  strcpy(cfilename, filename);
  fits_parse_rootname(cfilename, rootname, status);
  CHECK_STATUS_RET(*status, catalog);

  // Split rootname into the file path and the file name.
  char* lastslash = strrchr(rootname, '/');
  if (NULL==lastslash) {
    catalog->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(catalog->filepath, *status, 
		   "memory allocation for filepath failed", catalog);
    catalog->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(catalog->filename, *status, 
		   "memory allocation for filename failed", catalog);
    strcpy(catalog->filepath, "");
    strcpy(catalog->filename, rootname);
  } else {
    lastslash++;
    catalog->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(catalog->filename, *status, 
		   "memory allocation for filename failed", catalog);
    strcpy(catalog->filename, lastslash);
      
    *lastslash='\0';
    catalog->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(catalog->filepath, *status, 
		   "memory allocation for filepath failed", catalog);
    strcpy(catalog->filepath, rootname);
  }
  // END of storing the filename and filepath.

  // Open the specified FITS file.
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, catalog);

  // Move to the right extension.
  fits_movnam_hdu(fptr, BINARY_TBL, "SRC_CAT", 0, status);
  CHECK_STATUS_RET(*status, catalog);

  do { // Error handling loop.
    // Get the column names.
    int csrc_id=0, csrc_name=0, cra=0, cdec=0, cimgrota=0, cimgscal=0,
      ce_min=0, ce_max=0, cflux=0, cspectrum=0, cimage=0, clightcur=0;
    // Required columns:
    fits_get_colnum(fptr, CASEINSEN, "SRC_ID", &csrc_id, status);
    fits_get_colnum(fptr, CASEINSEN, "RA", &cra, status);
    fits_get_colnum(fptr, CASEINSEN, "DEC", &cdec, status);
    fits_get_colnum(fptr, CASEINSEN, "E_MIN", &ce_min, status);
    fits_get_colnum(fptr, CASEINSEN, "E_MAX", &ce_max, status);
    fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
    fits_get_colnum(fptr, CASEINSEN, "SPECTRUM", &cspectrum, status);
    fits_get_colnum(fptr, CASEINSEN, "IMAGE", &cimage, status);
    fits_get_colnum(fptr, CASEINSEN, "LIGHTCUR", &clightcur, status);
    CHECK_STATUS_BREAK(*status);
    // Optional columns:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(fptr, CASEINSEN, "SRC_NAME", &csrc_name, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(fptr, CASEINSEN, "IMGROTA", &cimgrota, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(fptr, CASEINSEN, "IMGSCAL", &cimgscal, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Take care of the units. Determine conversion factors.
    char ura[SIMPUT_MAXSTR];
    read_unit(fptr, cra, ura, status);
    CHECK_STATUS_BREAK(*status);
    float fra = unit_conversion_rad(ura);
    if (0.==fra) {
      SIMPUT_ERROR("unknown units in RA column");
      *status=EXIT_FAILURE;
      break;
    }

    char udec[SIMPUT_MAXSTR];
    read_unit(fptr, cdec, udec, status);
    CHECK_STATUS_BREAK(*status);
    float fdec = unit_conversion_rad(udec);
    if (0.==fdec) {
      SIMPUT_ERROR("unknown units in DEC column");
      *status=EXIT_FAILURE;
      break;
    }

    float fimgrota=0;
    if (cimgrota>0) {
      char uimgrota[SIMPUT_MAXSTR];
      read_unit(fptr, cimgrota, uimgrota, status);
      CHECK_STATUS_BREAK(*status);
      fimgrota = unit_conversion_rad(uimgrota);
      if (0.==fimgrota) {
	SIMPUT_ERROR("unknown units in IMGROTA column");
	*status=EXIT_FAILURE;
	break;
      }
    }

    char ue_min[SIMPUT_MAXSTR];
    read_unit(fptr, ce_min, ue_min, status);
    CHECK_STATUS_BREAK(*status);
    float fe_min = unit_conversion_keV(ue_min);
    if (0.==fe_min) {
      SIMPUT_ERROR("unknown units in E_MIN column");
      *status=EXIT_FAILURE;
      break;
    }

    char ue_max[SIMPUT_MAXSTR];
    read_unit(fptr, ce_max, ue_max, status);
    CHECK_STATUS_BREAK(*status);
    float fe_max = unit_conversion_keV(ue_max);
    if (0.==fe_max) {
      SIMPUT_ERROR("unknown units in E_MAX column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(fptr, cflux, uflux, status);
    CHECK_STATUS_BREAK(*status);
    float fflux = unit_conversion_ergpspcm2(uflux);
    if (0.==fflux) {
      SIMPUT_ERROR("unknown units in FLUX column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of required entries.
    long nrows;
    fits_get_num_rows(fptr, &nrows, status);
    CHECK_STATUS_BREAK(*status);
    catalog->entries  = (SimputSourceEntry**)malloc(nrows*sizeof(SimputSourceEntry*));
    CHECK_NULL_BREAK(catalog->entries, *status, 
		     "memory allocation for catalog entries failed");
    catalog->nentries = (int)nrows;

    // Allocate memory for string buffers.
    char* src_name[1]={NULL};
    char* spectrum[1]={NULL};
    char* image[1]   ={NULL};
    char* lightcur[1]={NULL};
    src_name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(src_name[0], *status, 
		     "memory allocation for string buffer failed");
    spectrum[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(spectrum[0], *status, 
		     "memory allocation for string buffer failed");
    image[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(image[0], *status, 
		     "memory allocation for string buffer failed");
    lightcur[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(lightcur[0], *status, 
		     "memory allocation for string buffer failed");

    // Loop over all rows in the table.
    long ii;
    for (ii=0; ii<nrows; ii++) {
      int src_id=0;
      double ra=0., dec=0.;
      float imgrota=0., imgscal=1.;
      float e_min=0., e_max=0., flux=0.;
      
      // Read the data from the table.
      int anynul=0;
      fits_read_col(fptr, TUINT, csrc_id, ii+1, 1, 1, 
		    &src_id, &src_id, &anynul, status);

      if (csrc_name>0) {
	fits_read_col_str(fptr, csrc_name, ii+1, 1, 1, 
			  "", src_name, &anynul, status);
      } else {
	strcpy(src_name[0], "");
      }

      fits_read_col(fptr, TDOUBLE, cra, ii+1, 1, 1, &ra, &ra, &anynul, status);
      ra *= fra; // Convert to [rad].
      fits_read_col(fptr, TDOUBLE, cdec, ii+1, 1, 1, &dec, &dec, &anynul, status);
      dec *= fdec; // Convert to [rad].

      if (cimgrota>0) {
	fits_read_col(fptr, TFLOAT, cimgrota, ii+1, 1, 1, 
		      &imgrota, &imgrota, &anynul, status);
	imgrota *= fimgrota; // Convert to [rad].
      }
      if (cimgscal>0) {
	fits_read_col(fptr, TFLOAT, cimgscal, ii+1, 1, 1, 
		      &imgscal, &imgscal, &anynul, status);
      }
      
      fits_read_col(fptr, TFLOAT, ce_min, ii+1, 1, 1, &e_min, &e_min, &anynul, status);
      e_min *= fe_min; // Convert to [keV].
      fits_read_col(fptr, TFLOAT, ce_max, ii+1, 1, 1, &e_max, &e_max, &anynul, status);
      e_max *= fe_max; // Convert to [keV].
      fits_read_col(fptr, TFLOAT, cflux, ii+1, 1, 1, &flux, &flux, &anynul, status);
      flux  *= fflux; // Convert to [erg/s/cm**2].

      fits_read_col(fptr, TSTRING, cspectrum, ii+1, 1, 1, 
		    "", spectrum, &anynul, status);
      fits_read_col(fptr, TSTRING, cimage, ii+1, 1, 1, 
		    "", image, &anynul, status);
      fits_read_col(fptr, TSTRING, clightcur, ii+1, 1, 1, 
		    "", lightcur, &anynul, status);

      CHECK_STATUS_BREAK(*status);

      // Add a new entry to the catalog.
      catalog->entries[ii] = 
	getSimputSourceEntryV(src_id, src_name[0], ra, dec, imgrota, imgscal, 
			      e_min, e_max, flux, spectrum[0], image[0],
			      lightcur[0], status);
      CHECK_STATUS_BREAK(*status);

      // Set the pointers to the filename and filepath in the catalog
      // data structure.
      catalog->entries[ii]->filepath = &catalog->filepath;
      catalog->entries[ii]->filename = &catalog->filename;

    }
    // END of loop over all rows in the table.

    // Release allocated memory.
    if (NULL!=src_name[0]) free(src_name[0]);
    if (NULL!=spectrum[0]) free(spectrum[0]);
    if (NULL!=image[0])    free(image[0]);
    if (NULL!=lightcur[0]) free(lightcur[0]);

    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.
  
  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, catalog);

  return(catalog);
}


void saveSimputSourceCatalog(const SimputSourceCatalog* const catalog,
			     const char* const filename,
			     int* const status)
{
  fitsfile* fptr=NULL;
  
  do { // Error handling loop.

    // Check if the file already exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // The file already exists.
      // Open the file and check, whether it contains a source catalog.
      fits_open_file(&fptr, filename, READWRITE, status);
      CHECK_STATUS_BREAK(*status);
      
      int status2=EXIT_SUCCESS;
      fits_write_errmark();
      fits_movnam_hdu(fptr, BINARY_TBL, "SRC_CAT", 0, &status2);
      if (BAD_HDU_NUM!=status2) {
	// The file alreay contains a source catalog.
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "the file '%s' already contains a source catalog", filename);
	SIMPUT_ERROR(msg);
	*status=EXIT_FAILURE;
	break;
      }
      fits_clear_errmark();

    } else {
      // The does not exist yet.
      // Create and open a new empty FITS file.
      fits_create_file(&fptr, filename, status);
      CHECK_STATUS_BREAK(*status);
    }


    // Create a binary table.
    const int csrc_id   = 1;
    const int csrc_name = 2;
    const int cra       = 3;
    const int cdec      = 4;
    const int cimgrota  = 5;
    const int cimgscal  = 6;
    const int ce_min    = 7;
    const int ce_max    = 8;
    const int cflux     = 9;
    const int cspectrum = 10;
    const int cimage    = 11;
    const int clightcur = 12;
    char *ttype[] = { "SRC_ID", "SRC_NAME", "RA", "DEC", "IMGROTA", "IMGSCAL", 
		      "E_MIN", "E_MAX", "FLUX", "SPECTRUM", "IMAGE", "LIGHTCUR" };
    char *tform[] = { "I", "1PA", "D", "D", "E", "E", 
		      "E", "E", "E", "1PA", "1PA", "1PA" };
    char *tunit[] = { "", "", "deg", "deg", "deg", "",  
		      "keV", "keV", "erg/s/cm**2", "", "", "" };
    // Provide option to use different units?
    fits_create_tbl(fptr, BINARY_TBL, 0, 12, ttype, tform, tunit, "SRC_CAT", status);
    CHECK_STATUS_BREAK(*status);

    // Write the necessary header keywords.
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS2", "SRC_CAT", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0", "", status);
    fits_write_key(fptr, TSTRING, "RADESYS", "FK5", "", status);
    float equinox=2000.0;
    fits_update_key(fptr, TFLOAT, "EQUINOX", &equinox, "", status);
    CHECK_STATUS_BREAK(*status);

    // Write the data.
    fits_insert_rows(fptr, 0, catalog->nentries, status);
    CHECK_STATUS_BREAK(*status);
    int ii;
    for (ii=0; ii<catalog->nentries; ii++) {
      fits_write_col(fptr, TUINT, csrc_id, ii+1, 1, 1, 
		     &catalog->entries[ii]->src_id, status);
      fits_write_col(fptr, TSTRING, csrc_name, ii+1, 1, 1, 
		     &catalog->entries[ii]->src_name, status);
      double ra = catalog->entries[ii]->ra*180./M_PI;
      fits_write_col(fptr, TDOUBLE, cra, ii+1, 1, 1, &ra, status);
      double dec = catalog->entries[ii]->dec*180./M_PI;
      fits_write_col(fptr, TDOUBLE, cdec, ii+1, 1, 1, &dec, status);
      float imgrota = catalog->entries[ii]->imgrota*180./M_PI;
      fits_write_col(fptr, TFLOAT, cimgrota, ii+1, 1, 1, &imgrota, status);
      fits_write_col(fptr, TFLOAT, cimgscal, ii+1, 1, 1, 
		     &catalog->entries[ii]->imgscal, status);
      fits_write_col(fptr, TFLOAT, ce_min, ii+1, 1, 1, 
		     &catalog->entries[ii]->e_min, status);
      fits_write_col(fptr, TFLOAT, ce_max, ii+1, 1, 1, 
		     &catalog->entries[ii]->e_max, status);
      fits_write_col(fptr, TFLOAT, cflux, ii+1, 1, 1, 
		     &catalog->entries[ii]->flux, status);
      fits_write_col(fptr, TSTRING, cspectrum, ii+1, 1, 1, 
		     &catalog->entries[ii]->spectrum, status);
      fits_write_col(fptr, TSTRING, cimage, ii+1, 1, 1, 
		     &catalog->entries[ii]->image, status);
      fits_write_col(fptr, TSTRING, clightcur, ii+1, 1, 1, 
		     &catalog->entries[ii]->lightcur, status);
      CHECK_STATUS_BREAK(*status);
    }
    CHECK_STATUS_BREAK(*status);
    // END of loop over all entries in the catalog.

  } while(0); // END of error handling loop.

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


static int check_if_btbl(const char* const filename,
			 int* const status)
{
  // Check if this is a FITS image.
  fitsfile* fptr=NULL;
  int hdutype;

  do { // Beginning of error handling loop.

    fits_open_file(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);
    
    fits_get_hdu_type(fptr, &hdutype, status);
    CHECK_STATUS_BREAK(*status);

  } while (0); // END of error handling loop.

  if (NULL!=fptr) fits_close_file(fptr, status);

  if (BINARY_TBL==hdutype) {
    return(1);
  } else {
    return(0);
  }
}


static void read_unit(fitsfile* const fptr, const int column, 
		      char* unit, int* const status)
{
  // Read the header keyword.
  char keyword[SIMPUT_MAXSTR], comment[SIMPUT_MAXSTR];
  sprintf(keyword, "TUNIT%d", column);
  fits_read_key(fptr, TSTRING, keyword, unit, comment, status);
  CHECK_STATUS_VOID(*status);
}


static void strtolower(char* const string)
{
  int ii=0;
  while (string[ii]!='\0') {
    string[ii] = tolower(string[ii]);
    ii++;
  };
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
  } else {
    // Unknown units.
    return(0.);
  }
}


SimputMissionIndepSpec* getSimputMissionIndepSpec(int* const status)
{
  SimputMissionIndepSpec* spec=
    (SimputMissionIndepSpec*)malloc(sizeof(SimputMissionIndepSpec));
  CHECK_NULL_RET(spec, *status, 
		 "memory allocation for SimputMissionIndepSpec failed", spec);

  // Initialize elements.
  spec->nentries=0;
  spec->energy  =NULL;
  spec->flux    =NULL;
  spec->distribution=NULL;
  spec->name    =NULL;
  spec->fileref =NULL;

  return(spec);  
}


void freeSimputMissionIndepSpec(SimputMissionIndepSpec** const spec)
{
  if (NULL!=*spec) {
    if (NULL!=(*spec)->energy) {
      free((*spec)->energy);
    }
    if (NULL!=(*spec)->flux) {
      free((*spec)->flux);
    }
    if (NULL!=(*spec)->distribution) {
      free((*spec)->distribution);
    }
    if (NULL!=(*spec)->name) {
      free((*spec)->name);
    }
    if (NULL!=(*spec)->fileref) {
      free((*spec)->fileref);
    }
    free(*spec);
    *spec=NULL;
  }
}


SimputMissionIndepSpec* loadSimputMissionIndepSpec(const char* const filename,
						   int* const status)
{
  // Allocate memory for string buffers.
  char* name[1]={NULL};

  SimputMissionIndepSpec* spec = getSimputMissionIndepSpec(status);
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
    if (nenergy!=nflux) {
      SIMPUT_ERROR("number of energy and flux entries in spectrum is not equivalent");
      *status=EXIT_FAILURE;
      break;
    }
    spec->nentries = (int)nenergy;
    printf("spectrum '%s' contains %d data points\n", 
	   filename, spec->nentries);

    // Allocate memory for the arrays.
    spec->energy  = (float*)malloc(spec->nentries*sizeof(float));
    CHECK_NULL_BREAK(spec->energy, *status, 
		     "memory allocation for spectrum failed");
    spec->flux    = (float*)malloc(spec->nentries*sizeof(float));
    CHECK_NULL_BREAK(spec->flux, *status, 
		     "memory allocation for spectrum failed");

    // Allocate memory for string buffer.
    name[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(name[0], *status, 
		     "memory allocation for string buffer failed");

    // Read the data from the table.
    int anynul=0;
    fits_read_col(fptr, TFLOAT, cenergy, 1, 1, spec->nentries, 
		  0, spec->energy, &anynul, status);
    fits_read_col(fptr, TFLOAT, cflux, 1, 1, spec->nentries, 
		  0, spec->flux, &anynul, status);

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
      spec->flux[ii]   *= fflux;
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



void saveSimputMissionIndepSpec(SimputMissionIndepSpec* const spec,
				const char* const filename,
				char* const extname,
				int extver,
				int* const status)
{
  fitsfile* fptr=NULL;
  
  // String buffer.
  char* name[1]={NULL}; 

  do { // Error handling loop.

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
      char *ttype[] = { "ENERGY", "FLUX", "NAME" };
      char *tform[] = { "E", "E", "1PA" };
      char *tunit[] = { "keV", "photon/s/cm**2/keV", "" };
      fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, extname, status);
      CHECK_STATUS_BREAK(*status);

      // Write header keywords.
      fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC", "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS1", "SIMPUT", "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS2", "SPECTRUM", "", status);
      fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0", "", status);
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
    fits_write_col(fptr, TFLOAT, cenergy, nrows, 1, (long)spec->nentries, 
		   &spec->energy, status);
    fits_write_col(fptr, TFLOAT, cflux, nrows, 1, (long)spec->nentries, 
		   &spec->flux, status);
    if ((cname>0) && (NULL!=spec->name)) {
      fits_write_col(fptr, TSTRING, cname, nrows, 1, 1, &spec->name, status);
    }
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=name[0]) free(name[0]);

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


void simputSetARF(struct ARF* const arf)
{
  static_arf = arf;
}


void simputSetRndGen(double(*rndgen)(void))
{
  static_rndgen=rndgen;
}


static SimputMissionIndepSpec* 
returnSimputMissionIndepSpec(const SimputSourceEntry* const src,
			     int* const status)
{
  const int maxspectra=10;
  static int nspectra=0;
  static SimputMissionIndepSpec** spectra=NULL;

  // Check, whether the source refers to a spectrum.
  if (NULL==src->spectrum) {
    SIMPUT_ERROR("source does not refer to a spectrum");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // In case there are no spectra available at all, allocate 
  // memory for the array (storage for spectra).
  if (NULL==spectra) {
    spectra = 
      (SimputMissionIndepSpec**)malloc(maxspectra*sizeof(SimputMissionIndepSpec*));
    CHECK_NULL_RET(spectra, *status, 
		   "memory allocation for spectra failed", NULL);
  }

  // Search if the required spectrum is available in the storage.
  int ii;
  for (ii=0; ii<nspectra; ii++) {
    // Check if the spectrum is equivalent to the required one.
    if (0==strcmp(spectra[ii]->fileref, src->spectrum)) {
      // If yes, determine a random photon energy from the spectral distribution.
      return(spectra[ii]);
    }
  }

  // The required spectrum is not contained in the storage.
  // Therefore we must load it from the specified location.
  if (nspectra>=maxspectra) {
    SIMPUT_ERROR("too many spectra in the internal storage");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // Load the mission-independent spectrum.
  char filename[SIMPUT_MAXSTR];
  if ('['==src->spectrum[0]) {
    strcpy(filename, *src->filepath);
    strcat(filename, *src->filename);
    strcat(filename, src->spectrum);
  } else {
    if ('/'!=src->spectrum[0]) {
      strcpy(filename, *src->filepath);
    } else {
      strcpy(filename, "");
    }
    strcat(filename, src->spectrum);
  }
  spectra[nspectra]=loadSimputMissionIndepSpec(filename, status);
  CHECK_STATUS_RET(*status, spectra[nspectra]);
  nspectra++;

  // Store the file reference to the spectrum for later comparisons.
  spectra[nspectra-1]->fileref = 
    (char*)malloc((strlen(src->spectrum)+1)*sizeof(char));
  CHECK_NULL_RET(spectra[nspectra-1]->fileref, *status, 
		 "memory allocation for file reference failed", 
		 spectra[nspectra-1]);
  strcpy(spectra[nspectra-1]->fileref, src->spectrum);

  // Multiply it by the ARF in order to obtain the spectral distribution.
  convSimputMissionIndepSpecWithARF(spectra[nspectra-1], status);
  CHECK_STATUS_RET(*status, spectra[nspectra-1]);
   
  return(spectra[nspectra-1]);
}


float getSimputPhotonEnergy(const SimputSourceEntry* const src,
			    int* const status)
{
  SimputMissionIndepSpec* spec=returnSimputMissionIndepSpec(src, status);
  CHECK_STATUS_RET(*status, 0.);

  // Determine a random photon energy from the spectral distribution.
  return(getRndPhotonEnergy(spec, status));
}


static float getRndPhotonEnergy(const SimputMissionIndepSpec* const spec,
				      int* const status) 
{
  int upper=spec->nentries, lower=0, mid;
  
  // Get a random number in the interval [0,1].
  float rnd = (float)static_rndgen();
  assert(rnd>=0.);
  assert(rnd<=1.);

  if (NULL==spec->distribution) {
    SIMPUT_ERROR("spectral distribution undefined");
    *status=EXIT_FAILURE;
    return(0.);
  }

  // Multiply with the total photon rate.
  rnd *= spec->distribution[spec->nentries-1];

  // Determine the energy of the photon (using binary search).
  while (upper>lower) {
    mid = (lower+upper)/2;
    if (spec->distribution[mid]<rnd) {
      lower = mid+1;
    } else {
      upper = mid;
    }
  }

  // Return the corresponding energy.
  if (0==lower) {
    return(spec->energy[0]*(float)static_rndgen());
  } else {
    return(spec->energy[lower-1] + 
	   (spec->energy[lower]-spec->energy[lower-1])*
	   (float)static_rndgen());
  }
}


void convSimputMissionIndepSpecWithARF(SimputMissionIndepSpec* const spec, 
				       int* const status)
{  
  // Check if the ARF is defined.
  CHECK_NULL_VOID(static_arf, *status, "instrument ARF undefined");

  // Allocate memory.
  spec->distribution = (float*)malloc(spec->nentries*sizeof(float));
  CHECK_NULL_VOID(spec->distribution, *status,
		 "memory allocation for spectral distribution failed");

  // Multiply each bin by the ARF and the width of the bin.
  // [photon/s/cm^2/keV] -> [photon/s]
  // The ARF contribution corresponding to a particular spectral bin 
  // is obtained by interpolation.
  float last_energy=0.; // [keV]
  int ii;
  for (ii=0; ii<spec->nentries; ii++) {
    // Determine the ARF contribution by interpolation.
    float arf_contribution=0.;
    long jj;
    for (jj=0; jj<static_arf->NumberEnergyBins; jj++) {
      if ((static_arf->LowEnergy[jj]<spec->energy[ii]) && 
	  (static_arf->HighEnergy[jj]>last_energy)) {
	float emin = MAX(static_arf->LowEnergy[jj], last_energy);
	float emax = MIN(static_arf->HighEnergy[jj], spec->energy[ii]);
	assert(emax>emin);
	arf_contribution += static_arf->EffArea[jj] * (emax-emin);
      }
    }

    spec->distribution[ii]=spec->flux[ii]*arf_contribution;
    last_energy           =spec->energy[ii];

    // Create the spectral distribution noramlized to the total 
    // photon rate [photon/s]. 
    if (ii>0) {
      spec->distribution[ii] += spec->distribution[ii-1];
    }
  }
}


static float getEbandFlux(const SimputSourceEntry* const src,
			  const float emin, const float emax,
			  int* const status)
{
  // Conversion factor from [keV] -> [erg].
  const float keV2erg = 1.602e-9;

  SimputMissionIndepSpec* spec=returnSimputMissionIndepSpec(src, status);
  CHECK_STATUS_RET(*status, 0.);

  int ii;
  float flux = 0.;
  float last_energy=0.;
  for (ii=0; ii<spec->nentries; ii++) {
    if ((spec->energy[ii]<emax) && (spec->energy[ii]>emin)) {
      float min = MAX(last_energy, emin);
      float max = MIN(spec->energy[ii], emax);
      assert(max>min);
      flux += (max-min) * spec->flux[ii] * spec->energy[ii];
    }
    last_energy=spec->energy[ii];
  }

  // Convert units of 'flux' from [keV/s/cm^2] -> [erg/s/cm^2].
  flux *= keV2erg;

  return(flux);
}


static float getEbandRate(const SimputSourceEntry* const src,
			  const float emin, const float emax,
			  int* const status)
{
  SimputMissionIndepSpec* spec=returnSimputMissionIndepSpec(src, status);
  CHECK_STATUS_RET(*status, 0.);

  int ii, upper=0;
  
  for (ii=spec->nentries-1; ii>=0; ii--) {
    if ((0==upper) && (spec->energy[ii]<=emax) && (spec->energy[ii]>emin)) {
      upper=ii;
    }
    if ((upper>0) && (spec->energy[ii]<emin)) {
      return(spec->distribution[upper]-spec->distribution[ii+1]);
    }
  }
  if (upper>0) {
    return(spec->distribution[upper]);
  }

  return(0.);
}


float getSimputPhotonRate(const SimputSourceEntry* const src,
			  int* const status)
{
  SimputMissionIndepSpec* spec=returnSimputMissionIndepSpec(src, status);
  CHECK_STATUS_RET(*status, 0.);

  return(src->flux / 
	 getEbandFlux(src, src->e_min, src->e_max, status) *
	 spec->distribution[spec->nentries-1]);
}


SimputLC* getSimputLC(int* const status)
{
  SimputLC* lc=(SimputLC*)malloc(sizeof(SimputLC));
  CHECK_NULL_RET(lc, *status, 
		 "memory allocation for SimputLC failed", lc);

  // Initialize elements.
  lc->nentries=0;
  lc->time    =NULL;
  lc->phase   =NULL;
  lc->flux    =NULL;
  lc->spectrum=NULL;
  lc->image   =NULL;
  lc->mjdref  =0.;
  lc->timezero=0.;
  lc->phase0  =0.;
  lc->period  =0.;
  lc->fluxscal=0.;
  lc->fileref =NULL;

  return(lc);
}


void freeSimputLC(SimputLC** const lc)
{
  if (NULL!=*lc) {
    if ((*lc)->nentries>0) {
      if (NULL!=(*lc)->spectrum) {
	int ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->spectrum[ii]) {
	    free((*lc)->spectrum[ii]);
	  }
	}
      }
      if (NULL!=(*lc)->image) {
	int ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->image[ii]) {
	    free((*lc)->image[ii]);
	  }
	}
      }
    }
    if (NULL!=(*lc)->time) {
      free((*lc)->time);
    }
    if (NULL!=(*lc)->phase) {
      free((*lc)->phase);
    }
    if (NULL!=(*lc)->flux) {
      free((*lc)->flux);
    }
    if (NULL!=(*lc)->fileref) {
      free((*lc)->fileref);
    }
    free(*lc);
    *lc=NULL;
  }
}


