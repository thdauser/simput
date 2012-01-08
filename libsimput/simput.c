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

/** The following macros are used to the store light curve and the PSD
    in the right format for the GSL routines. */
#define REAL(z,i) ((z)[(i)])
#define IMAG(z,i,n) ((z)[(n)-(i)])


/////////////////////////////////////////////////////////////////
// Structures.
/////////////////////////////////////////////////////////////////


struct SimputSourceBuffer {
  long nsources; // Current number of sources in the cash.
  long csource;  // Index of next position in cash that will be used.

  // Cash for the sources.
  SimputSource** sources;

  // This array contains the row numbers of the sources in the storage
  // given in the same order as the corresponding sources in the storage.
  // The array is used to replace the oldest source in the storage.
  long* rownumbers;

  // Array with a size corresponding to the number of entries 
  // in the catalog. Each entry in the array refers to the index
  // of the corresponding source in the storage. If the respective
  // source is not contained in the storage, the value in the array
  // is -1.
  long* rowmap;
};


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


/** Read the TUNITn keyword from the header of the current FITS HDU,
    where n denotes the specfied column. */
static void read_unit(fitsfile* const fptr, const int column, 
		      char* unit, int* const status);

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

/** Determine the factor required to convert the specified unit into
    [s]. If the conversion is not possible or implemented, the
    function return value is 0. */
static float unit_conversion_s(const char* const unit);

/** Determine the factor required to convert the specified unit into
    [Hz]. If the conversion is not possible or implemented, the
    function return value is 0. */
static float unit_conversion_Hz(const char* const unit);

/** Return the requested light curve. Keeps a certain number of light
    curves in an internal storage. If the requested light curve is not
    located in the internal storage, it is loaded from the reference
    given in the source catalog. If the the source does not refer to a
    light curve (i.e. it is a source with constant brightness) the
    function return value is NULL. */
static SimputLC* returnSimputLC(const SimputSource* const src,
				const double time, const double mjdref,
				int* const status);


/** Return the requested image. Keeps a certain number of images in an
    internal storage. If the requested image is not located in the
    internal storage, it is loaded from the reference given in the
    source catalog. If the the source does not refer to an image
    (i.e. it is a point-like source) the function return value is
    NULL. */
static SimputImg* returnSimputImg(const SimputSource* const src,
				  int* const status);


/////////////////////////////////////////////////////////////////
// Function Definitions.
/////////////////////////////////////////////////////////////////


SimputSource* getSimputSource(int* const status)
{
  SimputSource* entry=(SimputSource*)malloc(sizeof(SimputSource));
  CHECK_NULL_RET(entry, *status, 
		 "memory allocation for SimputSource failed", entry);

  // Initialize elements.
  entry->src_id  =0;
  entry->src_name=NULL;
  entry->ra      =0.;
  entry->dec     =0.;
  entry->imgrota =0.;
  entry->imgscal =1.;
  entry->e_min   =0.;
  entry->e_max   =0.;
  entry->eflux   =0.;
  entry->spectrum=NULL;
  entry->image   =NULL;
  entry->lightcur=NULL;
  entry->filename=NULL;
  entry->filepath=NULL;

  return(entry);
}


SimputSource* getSimputSourceV(const long src_id, 
			       const char* const src_name,
			       const double ra,
			       const double dec,
			       const float imgrota,
			       const float imgscal,
			       const float e_min,
			       const float e_max,
			       const float eflux,
			       const char* const spectrum,
			       const char* const image,
			       const char* const lightcur,
			       int* const status)
{
  SimputSource* entry=getSimputSource(status);
  CHECK_STATUS_RET(*status, entry);

  // Initialize with the given values.
  if (src_id<0) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "SRC_ID (%ld) must have a positive value", src_id);
    SIMPUT_ERROR(msg);
    *status=EXIT_FAILURE;
    return(entry);
  }
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
  entry->eflux   = eflux;

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


void freeSimputSource(SimputSource** const entry)
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


static struct SimputSourceBuffer* getSimputSourceBuffer(int* const status)
{
  struct SimputSourceBuffer *srcbuff = 
    (struct SimputSourceBuffer*)malloc(sizeof(struct SimputSourceBuffer));

  CHECK_NULL_RET(srcbuff, *status, 
		 "memory allocation for SimputSourceBuffer failed", srcbuff);
    
  srcbuff->nsources  =0;
  srcbuff->csource   =0;
  srcbuff->sources   =NULL;
  srcbuff->rownumbers=NULL;
  srcbuff->rowmap    =NULL;

  return(srcbuff);
}


static void freeSimputSourceBuffer(struct SimputSourceBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->sources) {
      long ii;
      for (ii=0; ii<(*sb)->nsources; ii++) {
	freeSimputSource(&((*sb)->sources[ii]));
      }
      free((*sb)->sources);
    }
    if (NULL!=(*sb)->rownumbers) {
      free((*sb)->rownumbers);
    }
    if (NULL!=(*sb)->rowmap) {
      free((*sb)->rowmap);
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputSource* loadSimputSource(SimputCatalog* const cf,
			       const long row,
			       int* const status)
{
  SimputSource* se=NULL;

  // String buffers.
  char* src_name[1]={NULL};
  char* spectrum[1]={NULL};
  char* image[1]   ={NULL};
  char* lightcur[1]={NULL};

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
    lightcur[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
    CHECK_NULL_BREAK(lightcur[0], *status, 
		     "memory allocation for string buffer failed");

    long src_id=0;
    double ra=0., dec=0.;
    float imgrota=0., imgscal=1.;
    float e_min=0., e_max=0., flux=0.;
      
    // Read the data from the table.
    int anynul=0;
    fits_read_col(cf->fptr, TLONG, cf->csrc_id, row, 1, 1, 
		  &src_id, &src_id, &anynul, status);
    
    if (cf->csrc_name>0) {
      fits_read_col_str(cf->fptr, cf->csrc_name, row, 1, 1, 
			"", src_name, &anynul, status);
    } else {
      strcpy(src_name[0], "");
    }

    fits_read_col(cf->fptr, TDOUBLE, cf->cra, row, 1, 1, &ra, &ra, &anynul, status);
    ra  *= cf->fra;  // Convert to [rad].
    fits_read_col(cf->fptr, TDOUBLE, cf->cdec, row, 1, 1, &dec, &dec, &anynul, status);
    dec *= cf->fdec; // Convert to [rad].

    if (cf->cimgrota>0) {
      fits_read_col(cf->fptr, TFLOAT, cf->cimgrota, row, 1, 1, 
		    &imgrota, &imgrota, &anynul, status);
      imgrota *= cf->fimgrota; // Convert to [rad].
    }
    if (cf->cimgscal>0) {
      fits_read_col(cf->fptr, TFLOAT, cf->cimgscal, row, 1, 1, 
		    &imgscal, &imgscal, &anynul, status);
    }
    
    fits_read_col(cf->fptr, TFLOAT, cf->ce_min, row, 1, 1, 
		  &e_min, &e_min, &anynul, status);
    e_min *= cf->fe_min; // Convert to [keV].
    fits_read_col(cf->fptr, TFLOAT, cf->ce_max, row, 1, 1, 
		  &e_max, &e_max, &anynul, status);
    e_max *= cf->fe_max; // Convert to [keV].
    fits_read_col(cf->fptr, TFLOAT, cf->cflux, row, 1, 1, 
		  &flux, &flux, &anynul, status);
    flux  *= cf->fflux; // Convert to [erg/s/cm**2].

    fits_read_col(cf->fptr, TSTRING, cf->cspectrum, row, 1, 1, 
		  "", spectrum, &anynul, status);
    fits_read_col(cf->fptr, TSTRING, cf->cimage, row, 1, 1, 
		  "", image, &anynul, status);
    fits_read_col(cf->fptr, TSTRING, cf->clightcur, row, 1, 1, 
		  "", lightcur, &anynul, status);

    CHECK_STATUS_BREAK(*status);

    // Create a new SimputSource data structure.
    se = getSimputSourceV(src_id, src_name[0], ra, dec, imgrota, imgscal, 
			  e_min, e_max, flux, spectrum[0], image[0],
			  lightcur[0], status);
    CHECK_STATUS_BREAK(*status);

    // Set the pointers to the filename and filepath in the
    // SimputCatalogFile data structure.
    se->filepath = &cf->filepath;
    se->filename = &cf->filename;

  } while(0); // END of error handling loop.

  // Release memory.
  if (NULL!=src_name[0]) free(src_name[0]);
  if (NULL!=spectrum[0]) free(spectrum[0]);
  if (NULL!=image[0])    free(image[0]);
  if (NULL!=lightcur[0]) free(lightcur[0]);

  return(se);
}


SimputSource* loadCacheSimputSource(SimputCatalog* const cf,
				    const long row,
				    int* const status)
{
  // Maximum number of sources in the cache.
  const long maxsources=1000000; 

  // Check if the source catalog contains a source buffer.
  if (NULL==cf->srcbuff) {
    cf->srcbuff = getSimputSourceBuffer(status);
    CHECK_STATUS_RET(*status, NULL);
  }

  // Convert the void* pointer to the source buffer into the right
  // format.
  struct SimputSourceBuffer* sb = (struct SimputSourceBuffer*)cf->srcbuff;


  // Allocate memory for the cache.
  if (NULL==sb->rowmap) {
    sb->rowmap=(long*)malloc(cf->nentries*sizeof(long));
    CHECK_NULL_RET(sb->rowmap, *status, 
		   "memory allocation for row map failed", NULL);
    long jj;
    for (jj=0; jj<cf->nentries; jj++) {
      sb->rowmap[jj] = -1;
    }
  }

  // Check if the cache already exists or whether this routine is 
  // called for the first time.
  if (NULL==sb->sources) {
    sb->sources=(SimputSource**)malloc(maxsources*sizeof(SimputSource*));
    CHECK_NULL_RET(sb->sources, *status,
		   "memory allocation for source cache failed", NULL);
  }
  if (NULL==sb->rownumbers) {
    sb->rownumbers=(long*)malloc(maxsources*sizeof(long));
    CHECK_NULL_RET(sb->rownumbers, *status,
		   "memory allocation for array of row numbers failed", NULL);
  }


  // Check if the requested row is already available in the cache.
  if (sb->rowmap[row-1]>0) {
    return(sb->sources[sb->rowmap[row-1]]);
  }

  // The requested source is not contained in the cache.
  // Therefore we must load it from the FITS file.
  
  // Check if the cache is already full.
  if (sb->nsources<maxsources) {
    sb->csource = sb->nsources;
    sb->nsources++;
  } else {
    sb->csource++;
    if (sb->csource>=maxsources) {
      sb->csource=0;
    }
    // Destroy the source that is currently stored at this place 
    // in the cache.
    freeSimputSource(&(sb->sources[sb->csource]));
    sb->rowmap[sb->rownumbers[sb->csource]-1] = -1;
    sb->rownumbers[sb->csource] = 0;
  }

  // Load the source from the FITS file.
  sb->sources[sb->csource]=loadSimputSource(cf, row, status);
  CHECK_STATUS_RET(*status, sb->sources[sb->csource]);
  sb->rownumbers[sb->csource]=row;
  sb->rowmap[row-1]=sb->csource;

  return(sb->sources[sb->csource]);
}


SimputCatalog* getSimputCatalog(int* const status)
{
  SimputCatalog* cf=
    (SimputCatalog*)malloc(sizeof(SimputCatalog));
  CHECK_NULL_RET(cf, *status, 
		 "memory allocation for SimputCatalog failed", cf);

  // Initialize elements.
  cf->fptr     =NULL;
  cf->nentries =0;
  cf->csrc_id  =0;
  cf->csrc_name=0;
  cf->cra      =0;
  cf->cdec     =0;
  cf->cimgrota =0;
  cf->cimgscal =0;
  cf->ce_min   =0;
  cf->ce_max   =0;
  cf->cflux    =0;
  cf->cspectrum=0;
  cf->cimage   =0;
  cf->clightcur=0;
  cf->fra      =0.;
  cf->fdec     =0.;
  cf->fimgrota =0.;
  cf->fflux    =0.;
  cf->fe_min   =0.;
  cf->fe_max   =0.;
  cf->filepath =NULL;
  cf->filename =NULL;
  cf->srcbuff  =NULL;

  return(cf);
}


void freeSimputCatalog(SimputCatalog** const cf,
		       int* const status)
{
  if (NULL!=*cf) {
    if (NULL!=(*cf)->fptr) {
      fits_close_file((*cf)->fptr, status);
    }
    if (NULL!=(*cf)->filepath) {
      free((*cf)->filepath);
    }
    if (NULL!=(*cf)->filename) {
      free((*cf)->filename);
    }
    if (NULL!=(*cf)->srcbuff) {
      freeSimputSourceBuffer((struct SimputSourceBuffer**)&((*cf)->srcbuff));
    }
    free(*cf);
    *cf=NULL;
  }
}


SimputCatalog* openSimputCatalog(const char* const filename,
				 const int mode,
				 int* const status)
{
  SimputCatalog* cf = getSimputCatalog(status);
  CHECK_STATUS_RET(*status, cf);

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
      cf->filepath=(char*)malloc(sizeof(char));
      CHECK_NULL_BREAK(cf->filepath, *status, 
		       "memory allocation for filepath failed");
      cf->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
      CHECK_NULL_BREAK(cf->filename, *status, 
		       "memory allocation for filename failed");
      strcpy(cf->filepath, "");
      strcpy(cf->filename, rootname);
    } else {
      lastslash++;
      cf->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
      CHECK_NULL_BREAK(cf->filename, *status, 
		       "memory allocation for filename failed");
      strcpy(cf->filename, lastslash);
      
      *lastslash='\0';
      cf->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
      CHECK_NULL_BREAK(cf->filepath, *status, 
		       "memory allocation for filepath failed");
      strcpy(cf->filepath, rootname);
    }
    // END of storing the filename and filepath.

    // Check if the file already exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    CHECK_STATUS_BREAK(*status);
    if (1==exists) {
      // The file already exists => open it.
      fits_open_file(&cf->fptr, filename, mode, status);
      CHECK_STATUS_BREAK(*status);

    } else if (READWRITE==mode) {
      // The file does not exist, but it shall be created.
      fits_create_file(&cf->fptr, filename, status);
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
    fits_movnam_hdu(cf->fptr, BINARY_TBL, "SRC_CAT", 0, &status2);
    fits_clear_errmark();
    if (BAD_HDU_NUM==status2) {
      if (READWRITE==mode) {
	// The file does not contain a source catalog => create one.
	char *ttype[] = { "SRC_ID", "SRC_NAME", "RA", "DEC", "IMGROTA", "IMGSCAL", 
			  "E_MIN", "E_MAX", "FLUX", "SPECTRUM", "IMAGE", "LIGHTCUR" };
	char *tform[] = { "J", "1PA", "D", "D", "E", "E", 
			  "E", "E", "E", "1PA", "1PA", "1PA" };
	char *tunit[] = { "", "", "deg", "deg", "deg", "",  
			  "keV", "keV", "erg/s/cm**2", "", "", "" };
	fits_create_tbl(cf->fptr, BINARY_TBL, 0, 12, ttype, tform, tunit, 
			"SRC_CAT", status);
	CHECK_STATUS_BREAK(*status);
	
	// Write the necessary header keywords.
	fits_write_key(cf->fptr, TSTRING, "HDUCLASS", "HEASARC", "", status);
	fits_write_key(cf->fptr, TSTRING, "HDUCLAS1", "SIMPUT", "", status);
	fits_write_key(cf->fptr, TSTRING, "HDUCLAS2", "SRC_CAT", "", status);
	fits_write_key(cf->fptr, TSTRING, "HDUVERS", "1.0.0", "", status);
	fits_write_key(cf->fptr, TSTRING, "RADESYS", "FK5", "", status);
	float equinox=2000.0;
	fits_update_key(cf->fptr, TFLOAT, "EQUINOX", &equinox, "", status);
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
    fits_get_colnum(cf->fptr, CASEINSEN, "SRC_ID", &cf->csrc_id, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "RA", &cf->cra, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "DEC", &cf->cdec, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "E_MIN", &cf->ce_min, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "E_MAX", &cf->ce_max, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "FLUX", &cf->cflux, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "SPECTRUM", &cf->cspectrum, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "IMAGE", &cf->cimage, status);
    fits_get_colnum(cf->fptr, CASEINSEN, "LIGHTCUR", &cf->clightcur, status);
    CHECK_STATUS_BREAK(*status);
    // Optional columns:
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_get_colnum(cf->fptr, CASEINSEN, "SRC_NAME", &cf->csrc_name, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(cf->fptr, CASEINSEN, "IMGROTA", &cf->cimgrota, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_get_colnum(cf->fptr, CASEINSEN, "IMGSCAL", &cf->cimgscal, &opt_status);
    opt_status=EXIT_SUCCESS;
    fits_clear_errmark();

    // Take care of the units. Determine conversion factors.
    char ura[SIMPUT_MAXSTR];
    read_unit(cf->fptr, cf->cra, ura, status);
    CHECK_STATUS_BREAK(*status);
    cf->fra = unit_conversion_rad(ura);
    if (0.==cf->fra) {
      SIMPUT_ERROR("unknown units in RA column");
      *status=EXIT_FAILURE;
      break;
    }

    char udec[SIMPUT_MAXSTR];
    read_unit(cf->fptr, cf->cdec, udec, status);
    CHECK_STATUS_BREAK(*status);
    cf->fdec = unit_conversion_rad(udec);
    if (0.==cf->fdec) {
      SIMPUT_ERROR("unknown units in DEC column");
      *status=EXIT_FAILURE;
      break;
    }

    if (cf->cimgrota>0) {
      char uimgrota[SIMPUT_MAXSTR];
      read_unit(cf->fptr, cf->cimgrota, uimgrota, status);
      CHECK_STATUS_BREAK(*status);
      cf->fimgrota = unit_conversion_rad(uimgrota);
      if (0.==cf->fimgrota) {
	SIMPUT_ERROR("unknown units in IMGROTA column");
	*status=EXIT_FAILURE;
	break;
      }
    }

    char ue_min[SIMPUT_MAXSTR];
    read_unit(cf->fptr, cf->ce_min, ue_min, status);
    CHECK_STATUS_BREAK(*status);
    cf->fe_min = unit_conversion_keV(ue_min);
    if (0.==cf->fe_min) {
      SIMPUT_ERROR("unknown units in E_MIN column");
      *status=EXIT_FAILURE;
      break;
    }

    char ue_max[SIMPUT_MAXSTR];
    read_unit(cf->fptr, cf->ce_max, ue_max, status);
    CHECK_STATUS_BREAK(*status);
    cf->fe_max = unit_conversion_keV(ue_max);
    if (0.==cf->fe_max) {
      SIMPUT_ERROR("unknown units in E_MAX column");
      *status=EXIT_FAILURE;
      break;
    }

    char uflux[SIMPUT_MAXSTR];
    read_unit(cf->fptr, cf->cflux, uflux, status);
    CHECK_STATUS_BREAK(*status);
    cf->fflux = unit_conversion_ergpspcm2(uflux);
    if (0.==cf->fflux) {
      SIMPUT_ERROR("unknown units in FLUX column");
      *status=EXIT_FAILURE;
      break;
    }
    // END of determine unit conversion factors.

    // Determine the number of entries.
    fits_get_num_rows(cf->fptr, &cf->nentries, status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  return(cf);
}


void appendSimputSource(SimputCatalog* const cf,
			SimputSource* const src,
			int* const status)
{
  // Write the data.
  fits_insert_rows(cf->fptr, cf->nentries, 1, status);
  CHECK_STATUS_VOID(*status);
  cf->nentries++;

  fits_write_col(cf->fptr, TLONG, cf->csrc_id, cf->nentries, 1, 1, 
		 &src->src_id, status);
  fits_write_col(cf->fptr, TSTRING, cf->csrc_name, cf->nentries, 1, 1, 
		 &src->src_name, status);
  double ra = src->ra*180./M_PI;
  fits_write_col(cf->fptr, TDOUBLE, cf->cra, cf->nentries, 1, 1, &ra, status);
  double dec = src->dec*180./M_PI;
  fits_write_col(cf->fptr, TDOUBLE, cf->cdec, cf->nentries, 1, 1, &dec, status);
  float imgrota = src->imgrota*180./M_PI;
  fits_write_col(cf->fptr, TFLOAT, cf->cimgrota, cf->nentries, 1, 1, &imgrota, status);
  fits_write_col(cf->fptr, TFLOAT, cf->cimgscal, cf->nentries, 1, 1, 
		 &src->imgscal, status);
  fits_write_col(cf->fptr, TFLOAT, cf->ce_min, cf->nentries, 1, 1, 
		 &src->e_min, status);
  fits_write_col(cf->fptr, TFLOAT, cf->ce_max, cf->nentries, 1, 1, 
		 &src->e_max, status);
  fits_write_col(cf->fptr, TFLOAT, cf->cflux, cf->nentries, 1, 1, 
		 &src->eflux, status);
  fits_write_col(cf->fptr, TSTRING, cf->cspectrum, cf->nentries, 1, 1, 
		 &src->spectrum, status);
  fits_write_col(cf->fptr, TSTRING, cf->cimage, cf->nentries, 1, 1, 
		 &src->image, status);
  fits_write_col(cf->fptr, TSTRING, cf->clightcur, cf->nentries, 1, 1, 
		 &src->lightcur, status);
  CHECK_STATUS_VOID(*status);
}


void appendSimputSourceBlock(SimputCatalog* const cf,
			     SimputSource** const src,
			     const long nsources,
			     int* const status)
{
  // Insert new rows.
  fits_insert_rows(cf->fptr, cf->nentries, nsources, status);
  CHECK_STATUS_VOID(*status);

  long first=cf->nentries+1;
  cf->nentries+=nsources;

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
    src_id = (long*)malloc(nsources*sizeof(long));
    CHECK_NULL_BREAK(src_id, *status, "memory allocation failed");
    ra  = (double*)malloc(nsources*sizeof(double));
    CHECK_NULL_BREAK(ra, *status, "memory allocation failed");
    dec = (double*)malloc(nsources*sizeof(double));
    CHECK_NULL_BREAK(dec, *status, "memory allocation failed");
    imgrota = (float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(imgrota, *status, "memory allocation failed");
    imgscal = (float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(imgscal, *status, "memory allocation failed");
    e_min   = (float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(e_min, *status, "memory allocation failed");
    e_max   = (float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(e_max, *status, "memory allocation failed");
    eflux   = (float*)malloc(nsources*sizeof(float));
    CHECK_NULL_BREAK(eflux, *status, "memory allocation failed");

    
    // Write the data.
    long ii;
    for (ii=0; ii<nsources; ii++) {

      // Copy values to buffers.
      src_id[ii] = src[ii]->src_id;
      ra[ii]     = src[ii]->ra *180./M_PI;
      dec[ii]    = src[ii]->dec*180./M_PI;
      imgrota[ii]= src[ii]->imgrota;
      imgscal[ii]= src[ii]->imgscal;
      e_min[ii]  = src[ii]->e_min;
      e_max[ii]  = src[ii]->e_max;
      eflux[ii]  = src[ii]->eflux;

      // Write strings.
      fits_write_col(cf->fptr, TSTRING, cf->csrc_name, first+ii, 1, 1, 
		     &(src[ii]->src_name), status);
      fits_write_col(cf->fptr, TSTRING, cf->cspectrum, first+ii, 1, 1, 
		     &(src[ii]->spectrum), status);
      fits_write_col(cf->fptr, TSTRING, cf->cimage, first+ii, 1, 1, 
		     &(src[ii]->image), status);
      fits_write_col(cf->fptr, TSTRING, cf->clightcur, first+ii, 1, 1, 
		     &(src[ii]->lightcur), status);
      CHECK_STATUS_BREAK(*status);
    }
    CHECK_STATUS_BREAK(*status);
    // END of loop over all sources.

    // Write the buffers to the file.
    fits_write_col(cf->fptr, TLONG, cf->csrc_id, first, 1, nsources, 
		   src_id, status);
    fits_write_col(cf->fptr, TDOUBLE, cf->cra, first, 1, nsources, 
		   ra, status);
    fits_write_col(cf->fptr, TDOUBLE, cf->cdec, first, 1, nsources, 
		   dec, status);
    fits_write_col(cf->fptr, TFLOAT, cf->cimgrota, first, 1, nsources, 
		   imgrota, status);
    fits_write_col(cf->fptr, TFLOAT, cf->cimgscal, first, 1, nsources, 
		   imgscal, status);
    fits_write_col(cf->fptr, TFLOAT, cf->ce_min, first, 1, nsources, 
		   e_min, status);
    fits_write_col(cf->fptr, TFLOAT, cf->ce_max, first, 1, nsources, 
		   e_max, status);
    fits_write_col(cf->fptr, TFLOAT, cf->cflux, first, 1, nsources, 
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
      (0==strcmp(unit, "Hz^-1")) || (0==strcmp(unit, "Hz**-1")) ){
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


SimputMIdpSpec* getSimputMIdpSpec(int* const status)
{
  SimputMIdpSpec* spec=
    (SimputMIdpSpec*)malloc(sizeof(SimputMIdpSpec));
  CHECK_NULL_RET(spec, *status, 
		 "memory allocation for SimputMIdpSpec failed", spec);

  // Initialize elements.
  spec->nentries=0;
  spec->energy  =NULL;
  spec->pflux   =NULL;
  spec->refflux =NULL;
  spec->distribution=NULL;
  spec->name    =NULL;
  spec->fileref =NULL;

  return(spec);  
}


void freeSimputMIdpSpec(SimputMIdpSpec** const spec)
{
  if (NULL!=*spec) {
    if (NULL!=(*spec)->energy) {
      free((*spec)->energy);
    }
    if (NULL!=(*spec)->pflux) {
      free((*spec)->pflux);
    }
    if (NULL!=(*spec)->refflux) {
      free((*spec)->refflux);
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


SimputMIdpSpec* loadSimputMIdpSpec(const char* const filename,
				   int* const status)
{
  // String buffer.
  char* name[1]={NULL};

  SimputMIdpSpec* spec=getSimputMIdpSpec(status);
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
      SIMPUT_ERROR("number of energy and flux entries in spectrum is not equivalent");
      *status=EXIT_FAILURE;
      break;
    }
    spec->nentries=nenergy;
    printf("spectrum '%s' contains %ld data points\n", 
	   filename, spec->nentries);

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
		  0, spec->energy, &anynul, status);
    fits_read_col(fptr, TFLOAT, cflux, 1, 1, spec->nentries, 
		  0, spec->pflux, &anynul, status);

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
      strcpy(tform[2], "32A");
      strcpy(tunit[2], "");

      // Create the table.
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
	// Check if the NAME string is too long.
	if (strlen(spec->name)>32) {
	  SIMPUT_ERROR("NAME value of spectrum contains more than 32 characters");
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


void getSimputSpectrumValue(const SimputMIdpSpec* const spec,
			    const long row,
			    float* const energy, 
			    float* const pflux,
			    int* const status)
{
  if (row>=spec->nentries) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("row number exceeds number of bins in spectrum");
    return;
  }

  *energy=spec->energy[row];
  *pflux =spec->pflux[row];
}


void simputSetARF(struct ARF* const arf)
{
  static_arf = arf;
}


void simputSetRndGen(double(*rndgen)(void))
{
  static_rndgen=rndgen;
}


/** Determine the time corresponding to a particular light curve bin
    [s]. The function takes into account, whether the light curve is
    periodic or not. For periodic light curves the specified number of
    periods is added to the time value. For non-periodic light curves
    the nperiod parameter is neglected. The returned value includes
    the MJDREF and TIMEZERO contributions. */
static double getLCTime(const SimputLC* const lc, 
			const long kk, 
			const long long nperiods,
			const double mjdref)
{
  if (NULL!=lc->time) {
    // Non-periodic light curve.
    return(lc->time[kk] + lc->timezero + (lc->mjdref-mjdref)*24.*3600.);
  } else {
    // Periodic light curve. 
    return((lc->phase[kk] - lc->phase0 + nperiods)*lc->period +
	   lc->timezero + (lc->mjdref-mjdref)*24.*3600.);    
  }
}


/** Determine the index of the bin in the light curve that corresponds
    to the specified time. */
static long getLCBin(const SimputLC* const lc, 
		     const double time, const double mjdref,
		     long long* nperiods, int* const status)
{
  // Check if the light curve is periodic or not.
  if (NULL!=lc->time) {
    // Non-periodic light curve.

    // Check if the requested time is within the covered interval.
    if ((time<getLCTime(lc, 0, 0, mjdref)) || 
	(time>=getLCTime(lc, lc->nentries-1, 0, mjdref))) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "requested time (%lf MJD) is outside the "
	      "interval covered by the light curve (%lf to %lf MJD)",
	      time/24./3600. + mjdref,
	      getLCTime(lc, 0, 0, 0.)/24./3600., 
	      getLCTime(lc, lc->nentries-1, 0, 0.)/24./3600.);
      SIMPUT_ERROR(msg);
      *status=EXIT_FAILURE;
      return(0);
    }
    
    *nperiods = 0;

  } else {
    // Periodic light curve.
    double dt = time-getLCTime(lc, 0, 0, mjdref);
    if (dt>=0.) {
      *nperiods = (long long)(dt/lc->period);
    } else {
      *nperiods = (long long)(dt/lc->period)-1;
    }      
  }

  // Determine the respective index kk of the light curve (using
  // binary search).
  long lower=0, upper=lc->nentries-2, mid;
  while (upper>lower) {
    mid = (lower+upper)/2;
    if (getLCTime(lc, mid+1, *nperiods, mjdref) < time) {
      lower = mid+1;
    } else {
      upper = mid;
    }
  }
  
  return(lower);
}


SimputMIdpSpec* returnSimputSrcSpec(const SimputSource* const src,
				    const double time, const double mjdref,
				    int* const status)
{
  // Reference to the requested spectrum.
  char fileref[SIMPUT_MAXSTR];
  strcpy(fileref, ""); // Initialize with empty string.

  // Check, whether the source refers to a light curve.
  SimputLC* lc=returnSimputLC(src, time, mjdref, status);
  CHECK_STATUS_RET(*status, NULL);

  // If the source has a light curve, check, whether it contains a
  // SPECTRUM column.
  if (NULL!=lc) {
    if (NULL!=lc->spectrum) {
      // Determine the current light curve bin.
      long long nperiods;
      long bin=getLCBin(lc, time, mjdref, &nperiods, status);
      CHECK_STATUS_RET(*status, NULL);

      // Copy the reference to the spectrum.
      strcpy(fileref, lc->spectrum[bin]);
    }
  }

  // If no spectral information has been found in the light curve,
  // check, whether the source entry in the catalog itself refers 
  // to a spectrum.
  if (0==strlen(fileref)) {
    if (NULL!=src->spectrum) {
      if ((strlen(src->spectrum)>0) && (0!=strcmp(src->spectrum, "NULL"))) {
	strcpy(fileref, src->spectrum);
      }
    }
  }

  // Check if any valid spectrum has been found.
  if (0==strlen(fileref)) {
    SIMPUT_ERROR("source does not refer to a spectrum");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // Check if the filename and filepath properties for the
  // source are properly set.
  CHECK_NULL_RET(src->filepath, *status,
		 "filepath of source catalog is not known", NULL);
  CHECK_NULL_RET(src->filename, *status,
		 "filename of source catalog is not known", NULL);

  // Load the mission-independent spectrum.
  char filename[SIMPUT_MAXSTR];
  if ('['==fileref[0]) {
    strcpy(filename, *src->filepath);
    strcat(filename, *src->filename);
    strcat(filename, fileref);
  } else {
    if ('/'!=fileref[0]) {
      strcpy(filename, *src->filepath);
    } else {
      strcpy(filename, "");
    }
    strcat(filename, fileref);
  }

  return(loadCacheSimputMIdpSpec(filename, status));
}


SimputMIdpSpec* loadCacheSimputMIdpSpec(const char* const filename,
					int* const status)
{
  const int maxspectra=2000; // Maximum number of spectra in storage.
  static int nspectra=0;     // Current number of spectra in storage.
  static int cspectrum=0;    // Index of next position in storage that will be used.

  static SimputMIdpSpec** spectra=NULL;

  // In case there are no spectra available at all, allocate 
  // memory for the array (storage for spectra).
  if (NULL==spectra) {
    spectra = 
      (SimputMIdpSpec**)malloc(maxspectra*sizeof(SimputMIdpSpec*));
    CHECK_NULL_RET(spectra, *status, 
		   "memory allocation for spectra failed", NULL);
  }

  // Search if the required spectrum is available in the storage.
  int ii;
  for (ii=0; ii<nspectra; ii++) {
    // Check if the spectrum is equivalent to the required one.
    if (0==strcmp(spectra[ii]->fileref, filename)) {
      // If yes, return the spectrum.
      return(spectra[ii]);
    }
  }

  // The required spectrum is not contained in the storage.
  // Therefore we must load it from the specified location.

  // Check if there is still space left in the spectral storage buffer.
  if (nspectra<maxspectra) {
    cspectrum=nspectra;
    nspectra++;
  } else {
    cspectrum++;
    if (cspectrum>=maxspectra) {
      cspectrum=0;
    }
    // Release the spectrum that is currently stored at this place in the
    // storage buffer.
    freeSimputMIdpSpec(&spectra[cspectrum]);
  }

  // Load the mission-independent spectrum.
  spectra[cspectrum]=loadSimputMIdpSpec(filename, status);
  CHECK_STATUS_RET(*status, spectra[cspectrum]);

  // Store the file reference to the spectrum for later comparisons.
  spectra[cspectrum]->fileref = 
    (char*)malloc((strlen(filename)+1)*sizeof(char));
  CHECK_NULL_RET(spectra[cspectrum]->fileref, *status, 
		 "memory allocation for file reference failed", 
		 spectra[cspectrum]);
  strcpy(spectra[cspectrum]->fileref, filename);

  return(spectra[cspectrum]);
}


static inline void getSpecEbounds(const SimputMIdpSpec* const spec,
				  const long idx, 
				  float* emin, float* emax)
{
  // Determine the lower boundary.
  if (idx>0) {
    *emin = 0.5*(spec->energy[idx]+spec->energy[idx-1]);
  } else {
    *emin = spec->energy[idx];
  }

  // Determine the upper boundary.
  if (idx<spec->nentries-1) {
    *emax = 0.5*(spec->energy[idx+1]+spec->energy[idx]);
  } else {
    *emax = spec->energy[idx];
  }
}


/** Convolve the given mission-independent spectrum with the
    instrument ARF. The product of this process is the spectral
    probability distribution binned to the energy grid of the ARF. */
static void convSimputMIdpSpecWithARF(SimputMIdpSpec* const spec, 
				      int* const status)
{
  // Check if the ARF is defined.
  CHECK_NULL_VOID(static_arf, *status, "instrument ARF undefined");

  // Allocate memory.
  spec->distribution=
    (double*)malloc(static_arf->NumberEnergyBins*sizeof(double));
  CHECK_NULL_VOID(spec->distribution, *status,
		 "memory allocation for spectral distribution failed");

  // Loop over all bins of the ARF.
  long ii, jj=0;
  int warning_printed=0; // Flag whether warning has been printed.
  for (ii=0; ii<static_arf->NumberEnergyBins; ii++) {
    // Initialize with 0.
    spec->distribution[ii]=0.;

    // Lower boundary of the current bin.
    float lo = static_arf->LowEnergy[ii];

    // Loop over all spectral points within the ARF bin.
    int finished=0;
    do {
      // Determine the next spectral point.
      float spec_emin=0., spec_emax=0.;
      for ( ; jj<spec->nentries; jj++) {
	getSpecEbounds(spec, jj, &spec_emin, &spec_emax);
	if (spec_emax>lo) break;
      }
      
      // Check special cases.
      if ((0==jj) && (spec_emin>static_arf->LowEnergy[ii])) {
	if (0==warning_printed) {
	  printf("*** warning: the spectrum does not cover the "
		 "full energy range of the ARF!\n");
	  warning_printed=1;
	}
	if (spec_emin>static_arf->HighEnergy[ii]) break;

      } else if (jj==spec->nentries) {
	if (0==warning_printed) {
	  printf("*** warning: the spectrum does not cover the "
		 "full energy range of the ARF!\n");
	  warning_printed=1;
	}
	break;
      }
      
      //if (spec_emin>lo) {
      //printf("*** warning: bins do not overlap (delta=%e)!\n",
      //       spec_emin-lo);
      //}

      // Upper boundary of the current bin.
      float hi;
      if (spec_emax<=static_arf->HighEnergy[ii]) {
	hi=spec_emax;
      } else {
	hi=static_arf->HighEnergy[ii];
	finished=1;
      }

      // Add to the spectral probability density.
      spec->distribution[ii]+=
	(hi-lo)*static_arf->EffArea[ii]*spec->pflux[jj];
      
      // Increase the lower boundary.
      lo=hi;

    } while (0==finished);

    // Create the spectral distribution function 
    // normalized to the total photon rate [photon/s]. 
    if (ii>0) {
      spec->distribution[ii]+=spec->distribution[ii-1];
    }

  } // Loop over all ARF bins.
}


/** Determine a random photon energy [keV] according to the specified
    spectral distribution. */
static float getRndPhotonEnergy(SimputMIdpSpec* const spec,
				int* const status) 
{
  // Check if the spectrum has been convolved with the ARF.
  if (NULL==spec->distribution) {
    // Multiply it by the ARF in order to obtain the spectral distribution.
    convSimputMIdpSpecWithARF(spec, status);
    CHECK_STATUS_RET(*status, 0.);
  }

  // Get a random number in the interval [0,1].
  double rnd = static_rndgen();
  assert(rnd>=0.);
  assert(rnd<=1.);

  // Multiply with the total photon rate.
  rnd *= spec->distribution[static_arf->NumberEnergyBins-1];

  // Determine the corresponding point in the spectral 
  // distribution (using binary search).
  long upper=static_arf->NumberEnergyBins-1, lower=0, mid;
  while (upper>lower) {
    mid = (lower+upper)/2;
    if (spec->distribution[mid]<rnd) {
      lower = mid+1;
    } else {
      upper = mid;
    }
  }

  // Return the corresponding photon energy.
  return(static_arf->LowEnergy[lower] + 
	 static_rndgen()*
	 (static_arf->HighEnergy[lower]-static_arf->LowEnergy[lower]));
}


float getSimputPhotonEnergy(const SimputSource* const src,
			    const double time,
			    const double mjdref,
			    int* const status)
{
  // Get the spectrum which is stored in the catalog.
  SimputMIdpSpec* spec=
    returnSimputSrcSpec(src, time, mjdref, status);
  CHECK_STATUS_RET(*status, 0.);

  // Determine a random photon energy from the spectral distribution.
  return(getRndPhotonEnergy(spec, status));
}


float getSimputSpecBandFlux(SimputMIdpSpec* const spec,
			    const float emin, const float emax,
			    int* const status)
{
  if (NULL!=spec->refflux) {
    if ((fabs(emin-spec->refflux->emin)/emin<1.e-6) &&
	(fabs(emax-spec->refflux->emax)/emax<1.e-6)) {
      return(spec->refflux->flux);
    }
  } else {
    spec->refflux=(SimputSpecBandFlux*)malloc(sizeof(SimputSpecBandFlux));
    CHECK_NULL_RET(spec->refflux, *status, 
		 "memory allocation for SimputSpecBandFlux failed", 0.);
  }

  // Conversion factor from [keV] -> [erg].
  const float keV2erg = 1.602e-9;

  // Return value.
  float flux=0.;

  long ii;
  for (ii=0; ii<spec->nentries; ii++) {
    float binmin, binmax;
    getSpecEbounds(spec, ii, &binmin, &binmax);
    if ((emin<binmax) && (emax>binmin)) {
      float min = MAX(binmin, emin);
      float max = MIN(binmax, emax);
      assert(max>min);
      flux += (max-min)*spec->pflux[ii]*spec->energy[ii];
    }
  }

  // Convert units of 'flux' from [keV/s/cm^2] -> [erg/s/cm^2].
  flux *= keV2erg;

  // Store the calculated values in the SimputSpecBandFlux 
  // data structure in order to avoid re-calculation,
  // when this function is called later again.
  spec->refflux->emin = emin;
  spec->refflux->emax = emax;
  spec->refflux->flux = flux;

  return(flux);
}


float getSimputSrcBandFlux(const SimputSource* const src,
			   const double time, const double mjdref,
			   int* const status)
{
  // Determine the spectrum valid for the specified point of time.
  SimputMIdpSpec* spec=
    returnSimputSrcSpec(src, time, mjdref, status);
  CHECK_STATUS_RET(*status, 0.);

  return(getSimputSpecBandFlux(spec, src->e_min, src->e_max, status));
}


float getSimputPhotonRate(const SimputSource* const src,
			  const double time, const double mjdref,
			  int* const status)
{
  SimputMIdpSpec* spec=
    returnSimputSrcSpec(src, time, mjdref, status);
  CHECK_STATUS_RET(*status, 0.);

  // Flux in the reference energy band.
  float refband_flux=
    getSimputSpecBandFlux(spec, src->e_min, src->e_max, status);
  CHECK_STATUS_RET(*status, 0.);

  // Check for ARF convolution. This is required, since spec->distribution 
  // is used in the next step.
  if (NULL==spec->distribution) {
    // Multiply it by the ARF in order to obtain the spectral distribution.
    convSimputMIdpSpecWithARF(spec, status);
    CHECK_STATUS_RET(*status, 0.);
  }

  return(src->eflux / refband_flux * 
	 (float)(spec->distribution[static_arf->NumberEnergyBins-1]));
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
  lc->a       =NULL;
  lc->b       =NULL;
  lc->spectrum=NULL;
  lc->image   =NULL;
  lc->mjdref  =0.;
  lc->timezero=0.;
  lc->phase0  =0.;
  lc->period  =0.;
  lc->fluxscal=0.;
  lc->src_id  =0;
  lc->fileref =NULL;

  return(lc);
}


void freeSimputLC(SimputLC** const lc)
{
  if (NULL!=*lc) {
    if ((*lc)->nentries>0) {
      if (NULL!=(*lc)->spectrum) {
	long ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->spectrum[ii]) {
	    free((*lc)->spectrum[ii]);
	  }
	}
	free((*lc)->spectrum);    
      }
      if (NULL!=(*lc)->image) {
	long ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->image[ii]) {
	    free((*lc)->image[ii]);
	  }
	}
	free((*lc)->image);
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
    if (NULL!=(*lc)->a) {
      free((*lc)->a);
    }
    if (NULL!=(*lc)->b) {
      free((*lc)->b);
    }
    if (NULL!=(*lc)->fileref) {
      free((*lc)->fileref);
    }
    free(*lc);
    *lc=NULL;
  }
}


static int isSimputLC(const char* const filename, 
		      int* const status)
{
  fitsfile* fptr=NULL;
  int ret=0; // Return value.

  do { // Beginning of error handling loop.

    // Open the specified FITS file. The filename must uniquely identify
    // the light curve contained in a binary table via the extended filename 
    // syntax. Therefore we do not have to care about the HDU number.
    fits_open_table(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);
    
    // Read the HDUCLAS1 and HDUCLAS2 header keywords.
    char comment[SIMPUT_MAXSTR];
    char hduclas1[SIMPUT_MAXSTR];
    char hduclas2[SIMPUT_MAXSTR];
    fits_read_key(fptr, TSTRING, "HDUCLAS1", &hduclas1, comment, status);
    fits_read_key(fptr, TSTRING, "HDUCLAS2", &hduclas2, comment, status);
    CHECK_STATUS_BREAK(*status);

    if ((0==strcmp(hduclas1, "SIMPUT")) && (0==strcmp(hduclas2, "LIGHTCUR"))) {
      // This is a SIMPUT light curve.
      ret = 1;
    }

  } while (0); // END of error handling loop.

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, ret);

  return(ret);
}


static void gauss_rndgen(double* const x, double* const y)
{
  double sqrt_2rho = sqrt(-log(static_rndgen())*2.);
  double phi = static_rndgen()*2.*M_PI;

  *x = sqrt_2rho * cos(phi);
  *y = sqrt_2rho * sin(phi);
}


static void setLCAuxValues(SimputLC* const lc, int* const status)
{
  // Memory allocation.
  lc->a       = (float*)malloc(lc->nentries*sizeof(float));
  CHECK_NULL_VOID(lc->a, *status, 
		 "memory allocation for light curve failed");
  lc->b       = (float*)malloc(lc->nentries*sizeof(float));
  CHECK_NULL_VOID(lc->b, *status, 
		  "memory allocation for light curve failed");

  // Determine the auxiliary values for the light curve (including
  // FLUXSCAL).
  long ii;
  for (ii=0; ii<lc->nentries-1; ii++) {
    double dt;
    if (NULL!=lc->time) {
      // Non-periodic light curve.
      dt = lc->time[ii+1]-lc->time[ii];
    } else {
      // Periodic light curve.
      dt = (lc->phase[ii+1]-lc->phase[ii])*lc->period;
    }
    lc->a[ii] = (lc->flux[ii+1]-lc->flux[ii])	
      /dt /lc->fluxscal;
    lc->b[ii] = lc->flux[ii]/lc->fluxscal; 
  }
  lc->a[lc->nentries-1] = 0.;
  lc->b[lc->nentries-1] = lc->flux[lc->nentries-1]/lc->fluxscal;
}


static SimputLC* loadSimputLCfromPSD(const char* const filename, 
				     const double t0,
				     const double mjdref,
				     int* const status)
{
  SimputPSD* psd=NULL;
  SimputLC*  lc =NULL;
  fitsfile* fptr=NULL;

  // Buffer for Fourier transform.
  double* fcomp =NULL;
  float* power  =NULL;

  printf("*** warning: the loadSimputLCfromPSD() routine has to be revised!\n");

  do { // Error handling loop.

    // Load the PSD.
    psd=loadSimputPSD(filename, status);
    CHECK_STATUS_BREAK(*status);

    // Check if the number of bins is a power of 2.
    long nentries=psd->nentries;
    while (0==nentries % 2) {
      nentries = nentries/2;
    }
    if (1!=nentries) {
      SIMPUT_ERROR("PSD length is not a power of 2");
      *status=EXIT_FAILURE;
      break;
    }

    // Get an empty SimputLC data structure.
    lc = getSimputLC(status);
    CHECK_STATUS_BREAK(*status);

    // Set the MJDREF.
    lc->mjdref = mjdref;

    // Allocate memory.
    lc->time    = (double*)malloc(2*psd->nentries*sizeof(double));
    CHECK_NULL_BREAK(lc->time, *status, 
		     "memory allocation for light curve failed");
    lc->flux    = (float*)malloc(2*psd->nentries*sizeof(float));
    CHECK_NULL_BREAK(lc->flux, *status, 
		     "memory allocation for light curve failed");
    lc->nentries=2*psd->nentries;

    // Convert the PSD into a light curve.
    // Set the time bins.
    lc->timezero = t0;
    long ii;
    for (ii=0; ii<lc->nentries; ii++) {
      lc->time[ii] = ii*1./(2.*psd->frequency[psd->nentries-1]);
    }

    // Buffer for PSD. Used to rescale in order to get proper RMS.
    power = (float*)malloc(psd->nentries*sizeof(float));
    CHECK_NULL_BREAK(power, *status, "memory allocation for light "
		     "curve (buffer for PSD transform) failed");

    // The PSD is given in Miyamoto normalization. In order to get the RMS
    // right, we have to multiply each bin with df (delta frequency).
    power[0] = psd->power[0]*psd->frequency[0];
    for (ii=1; ii<psd->nentries; ii++) {
      power[ii] = psd->power[ii]*(psd->frequency[ii]-psd->frequency[ii-1]);
    }

    // Create Fourier components.
    fcomp = (double*)malloc(lc->nentries*sizeof(double));
    CHECK_NULL_BREAK(fcomp, *status, "memory allocation for light "
		     "curve (buffer for Fourier transform) failed");
    
    // Apply the algorithm introduced by Timmer & Koenig (1995).
    double randr, randi;
    lc->fluxscal=1.; // Set Fluxscal to 1.
    gauss_rndgen(&randr, &randi);
    fcomp[0]             = 1.;
    fcomp[psd->nentries] = randi*sqrt(0.5*power[psd->nentries-1]);
    for (ii=1; ii<psd->nentries; ii++) {
      gauss_rndgen(&randr, &randi);
      REAL(fcomp, ii)               = randr*sqrt(0.5*power[ii-1]);
      IMAG(fcomp, ii, lc->nentries) = randi*sqrt(0.5*power[ii-1]);
    }

    // Perform Fourier (back-)transformation.
    gsl_fft_halfcomplex_radix2_backward(fcomp, 1, lc->nentries);

    // Normalization.
    // Calculate the required rms.
    float requ_rms=1.;
    for (ii=0; ii<psd->nentries; ii++) {
      // It is not absolutely clear to me, why we have to divide by 2 here.
      // But the results are only right, if we do that.
      requ_rms += power[ii]/2.;
    }
    requ_rms = sqrt(requ_rms);

    // Calculate the actual rms.
    float act_rms=0.;
    for (ii=0; ii<lc->nentries; ii++) {
      act_rms += (float)pow(fcomp[ii], 2.);
    }
    act_rms = sqrt(act_rms/lc->nentries);

    // Determine the normalized rates from the FFT.
    for (ii=0; ii<lc->nentries; ii++) {
      lc->flux[ii] = (float)fcomp[ii] * requ_rms/act_rms;

      // Avoid negative fluxes (no physical meaning):
      if (lc->flux[ii]<0.) { 
	lc->flux[ii] = 0.; 
      }
    }

    // Determine the auxiliary light curve parameters.
    setLCAuxValues(lc, status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.

  // Release allocated memory.
  if (NULL!=power) free(power);
  if (NULL!=fcomp) free(fcomp);
  freeSimputPSD(&psd);

  // Close the file.
  if (NULL!=fptr) fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, lc);
  
  // TODO For testing we can write the generated light curve to 
  // an default output file.
  //remove("lc.fits");
  //saveSimputLC(lc, "lc.fits", "LC", 1, status);

  return(lc);
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
    lc = getSimputLC(status);
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
      fits_read_key(fptr, TFLOAT,  "PHASE0",   &lc->phase0,   comment, status);
      fits_read_key(fptr, TFLOAT,  "PERIOD",   &lc->period,   comment, status);
    } else {
      lc->phase0 = 0.;
      lc->period = 0.;
    }
    CHECK_STATUS_BREAK(*status);
    // Optional keywords.
    opt_status=EXIT_SUCCESS;
    fits_write_errmark();
    fits_read_key(fptr, TFLOAT,  "FLUXSCAL", &lc->fluxscal, comment, &opt_status);
    if (EXIT_SUCCESS!=opt_status) {
      // FLUXSCAL is not given in the FITS header. We therefore assume
      // that it has a value of 1.
      lc->fluxscal = 1.;
    }
    fits_clear_errmark();


    // Determine the number of rows in the table.
    fits_get_num_rows(fptr, &lc->nentries, status);
    CHECK_STATUS_BREAK(*status);
    printf("light curve '%s' contains %ld data points\n", 
	   filename, lc->nentries);

    // Allocate memory for the arrays.
    if (ctime>0) {
      lc->time  = (double*)malloc(lc->nentries*sizeof(double));
      CHECK_NULL_BREAK(lc->time, *status, 
		       "memory allocation for light curve failed");
    }
    if (cphase>0) {
      lc->phase = (float*)malloc(lc->nentries*sizeof(float));
      CHECK_NULL_BREAK(lc->phase, *status, 
		       "memory allocation for light curve failed");
    }
    lc->flux    = (float*)malloc(lc->nentries*sizeof(float));
    CHECK_NULL_BREAK(lc->flux, *status, 
		     "memory allocation for light curve failed");
    if (cspectrum>0) {
      lc->spectrum = (char**)malloc(lc->nentries*sizeof(char*));
      CHECK_NULL_BREAK(lc->spectrum, *status, 
		       "memory allocation for light curve failed");
      // String buffer.
      spectrum[0]=(char*)malloc(SIMPUT_MAXSTR*sizeof(char));
      CHECK_NULL_BREAK(spectrum[0], *status, 
		       "memory allocation for string buffer failed");
    }
    if (cimage>0) {
      lc->spectrum = (char**)malloc(lc->nentries*sizeof(char*));
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
	lc->time[row] *= ftime;
      }
    }

    // PHASE
    if (cphase>0) {
      fits_read_col(fptr, TFLOAT, cphase, 1, 1, lc->nentries, 
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
	lc->spectrum[row] = 
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
	lc->image[row] = 
	  (char*)malloc((strlen(image[0])+1)*sizeof(char));
	CHECK_NULL_BREAK(lc->image[row], *status,
			 "memory allocation for image string failed");
	strcpy(lc->image[row], image[0]);
      }      
      CHECK_STATUS_BREAK(*status);
    }
    // END of reading the data from the FITS table.


    // Determine the auxiliary values for the light curve (including
    // FLUXSCAL).
    setLCAuxValues(lc, status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // END of error handling loop.
  
  // Release allocated memory.
  if (NULL!=spectrum[0]) free(spectrum[0]);
  if (NULL!=image[0])    free(image[0]);

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
    CHECK_NULL_BREAK(ttype, *status, "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tform, *status, "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tunit, *status, "memory allocation for string buffer failed");
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
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS2", "LIGHTCUR", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0", "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);
    fits_write_key(fptr, TDOUBLE, "MJDREF", &lc->mjdref, "", status);
    fits_write_key(fptr, TDOUBLE, "TIMEZERO", &lc->timezero, "", status);
    fits_write_key(fptr, TFLOAT,  "FLUXSCAL", &lc->fluxscal, "", status);
    int periodic=0;
    if (cphase>0) {
      // Only for periodic light curves.
      periodic=1;
      fits_write_key(fptr, TFLOAT,  "PHASE0", &lc->phase0, "", status);
      fits_write_key(fptr, TFLOAT,  "PERIOD", &lc->period, "", status);
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
      fits_write_col(fptr, TFLOAT, cphase, 1, 1, lc->nentries, 
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
  if (NULL!=spectrum[0]) free(spectrum[0]);
  if (NULL!=image[0])    free(image[0]);

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
    CHECK_NULL_BREAK(ttype, *status, "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tform, *status, "memory allocation for string buffer failed");
    CHECK_NULL_BREAK(tunit, *status, "memory allocation for string buffer failed");
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
      strcpy(ttype[0], "FREQUENC");
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
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS2", "POWSPEC", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0", "", status);
    fits_write_key(fptr, TINT,    "EXTVER", &extver, "", status);

    // Create new rows in the table and store the data of the power spectrum in it.
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

static double rndexp(const double avgdist)
{
  assert(avgdist>0.);

  double rand;
  do {
    rand=static_rndgen();
    assert(rand>=0.);
  } while (rand==0.);

  return(-log(rand)*avgdist);
}


double getSimputPhotonTime(const SimputSource* const src,
			   double prevtime,
			   const double mjdref,
			   int* const failed,
			   int* const status)
{
  // Determine the light curve.
  SimputLC* lc=returnSimputLC(src, prevtime, mjdref, status);
  CHECK_STATUS_RET(*status, 0.);

  // Set the 'failed' flag to its default value.
  *failed=0.;

  // Check, whether the source has constant brightness.
  if (NULL==lc) {
    // The source has a constant brightness.
    float rate=getSimputPhotonRate(src, prevtime, mjdref, status);
    if (0.==rate) {
      *failed=1;
      return(0.);
    } else {
      assert(rate>0.);
      return(prevtime+ rndexp((double)1./rate));
    }
  } else {
    // The source has a time-variable brightness.

    // The LinLightCurve contains data points, so the
    // general algorithm proposed by Klein & Roberts has to 
    // be applied.

    // Apply the individual photon rate of the particular source.
    float avgrate=getSimputPhotonRate(src, prevtime, mjdref, status);
    CHECK_STATUS_RET(*status, 0.);
    assert(avgrate>0.);

    // Step 1 in the algorithm.
    double u = static_rndgen();

    // Determine the respective index kk of the light curve.
    long long nperiods=0;
    long kk=getLCBin(lc, prevtime, mjdref, &nperiods, status);
    CHECK_STATUS_RET(*status, 0.);

    while ((kk<lc->nentries-1)||(lc->src_id>0)) {

      // If the end of the light curve is reached, check if it has
      // been produced from a PSD. In that case one can create new one.
      if ((kk>=lc->nentries-1)&&(lc->src_id>0)) {
	lc=returnSimputLC(src, prevtime, mjdref, status);
	CHECK_STATUS_RET(*status, 0.);
	kk=getLCBin(lc, prevtime, mjdref, &nperiods, status);
	CHECK_STATUS_RET(*status, 0.);
      }

      // Determine the relative time within the kk-th interval, i.e., t=0 lies
      // at the beginning of the kk-th interval.
      double t        =prevtime-(getLCTime(lc, kk, nperiods, mjdref));
      double stepwidth=
	getLCTime(lc, kk+1, nperiods, mjdref)-getLCTime(lc, kk, nperiods, mjdref);

      // Step 2 in the algorithm.
      double uk = 1.-exp((-lc->a[kk]/2.*(pow(stepwidth,2.)-pow(t,2.))
			  -lc->b[kk]*(stepwidth-t))*avgrate);
      // Step 3 in the algorithm.
      if (u <= uk) {
	if ( fabs(lc->a[kk]*stepwidth) > fabs(lc->b[kk]*1.e-6) ) { 
	  // Instead of checking if a_kk = 0. check, whether its product with the 
	  // interval length is a very small number in comparison to b_kk.
	  // If a_kk * stepwidth is much smaller than b_kk, the rate in the interval
	  // can be assumed to be approximately constant.
	  return(getLCTime(lc, kk, nperiods, mjdref) +
		 (-lc->b[kk]+sqrt(pow(lc->b[kk],2.) + pow(lc->a[kk]*t,2.) + 
				  2.*lc->a[kk]*lc->b[kk]*t - 
				  2.*lc->a[kk]*log(1.-u)/avgrate)
		  )/lc->a[kk]);

	} else { // a_kk == 0
	  return(prevtime-log(1.-u)/(lc->b[kk]*avgrate));
	}

      } else {
	// Step 4 (u > u_k).
	u = (u-uk)/(1-uk);
	kk++;
	if ((kk>=lc->nentries-1) && (NULL!=lc->phase)) {
	  kk=0;
	  nperiods++;
	}
	prevtime = getLCTime(lc, kk, nperiods, mjdref);
      }
    }
    // END of while (kk < lc->nentries).
    
    // The range of the light curve has been exceeded.
    *failed=1;
    return(0.);
  }
}


static SimputLC* returnSimputLC(const SimputSource* const src,
				const double time, const double mjdref,
				int* const status)
{
  const int maxlcs=10; // Maximum number of light curves in storage.
  static int nlcs=0;   // Current number of light curves in storage.
  static int clc =0;   // Index of next position in storage that will be used.

  static SimputLC** lcs=NULL;

  // Check, whether the source refers to a light curve.
  if (NULL==src->lightcur) {
    return(NULL);
  }
  if (0==strlen(src->lightcur)) {
    return(NULL);
  }
  if (0==strcmp(src->lightcur, "NULL")) {
    return(NULL);
  }

  // In case there are no light curves available at all, allocate 
  // memory for the array (storage for light curves).
  if (NULL==lcs) {
    lcs = (SimputLC**)malloc(maxlcs*sizeof(SimputLC*));
    CHECK_NULL_RET(lcs, *status, 
		   "memory allocation for light curves failed", NULL);
  }

  // Check if the requested light curve is available in the storage.
  int ii;
  for (ii=0; ii<nlcs; ii++) {
    // Check if the light curve is equivalent to the requested one.
    if (0==strcmp(lcs[ii]->fileref, src->lightcur)) {
      // For a light curve created from a PSD, we also have to check,
      // whether this is the source associated to this light curve.
      if (lcs[ii]->src_id>0) {
	// Check if the SRC_IDs agree.
	if (lcs[ii]->src_id==src->src_id) {
	  // We have a light curve which has been produced from a PSD.
	  // Check if the requested time is covered by the light curve.
	  if (time<getLCTime(lcs[ii], lcs[ii]->nentries-1, 0, mjdref)) {
	    return(lcs[ii]);
	  }
	  // If not, we have to produce a new light curve from the PSD.
	}
      } else {
	// This light curve is loaded from a file and can be re-used
	// for different sources.
	return(lcs[ii]);
      }
    }
  }

  // The requested light curve is not contained in the storage.
  // Therefore we must load it from the specified location.

  // Check if there is still space left in the light curve storage.
  if (nlcs<maxlcs) {
    clc = nlcs;
    nlcs++;
  } else {
    clc++;
    if (clc>=maxlcs) {
      clc=0;
    }
    // Release the light curve that is currently stored at this place
    // in the cache.
    freeSimputLC(&lcs[clc]);
  }

  // Check if the filename and filepath properties for the
  // source are properly set.
  CHECK_NULL_RET(src->filepath, *status,
		 "filepath of source catalog is not known", NULL);
  CHECK_NULL_RET(src->filename, *status,
		 "filename of source catalog is not known", NULL);

  // Load the light curve.
  char filename[SIMPUT_MAXSTR];
  if ('['==src->lightcur[0]) {
    strcpy(filename, *src->filepath);
    strcat(filename, *src->filename);
    strcat(filename, src->lightcur);
  } else {
    if ('/'!=src->lightcur[0]) {
      strcpy(filename, *src->filepath);
    } else {
      strcpy(filename, "");
    }
    strcat(filename, src->lightcur);
  }

  // Check, whether the reference refers to a light curve or to 
  // a PSD.
  int islc = isSimputLC(filename, status);
  CHECK_STATUS_RET(*status, lcs[nlcs]);
  if (1==islc) {
    lcs[clc]=loadSimputLC(filename, status);
    CHECK_STATUS_RET(*status, lcs[clc]);
  } else {
    lcs[clc]=loadSimputLCfromPSD(filename, time, mjdref, status);
    CHECK_STATUS_RET(*status, lcs[clc]);
    // The new LC is assigned to this particular source and 
    // cannot be re-used.
    lcs[clc]->src_id=src->src_id;
  }

  // Store the file reference to the light curve for later comparisons.
  lcs[clc]->fileref = 
    (char*)malloc((strlen(src->lightcur)+1)*sizeof(char));
  CHECK_NULL_RET(lcs[clc]->fileref, *status, 
		 "memory allocation for file reference failed", 
		 lcs[clc]);
  strcpy(lcs[clc]->fileref, src->lightcur);
   
  return(lcs[clc]);
}


SimputImg* getSimputImg(int* const status)
{
  SimputImg* img=(SimputImg*)malloc(sizeof(SimputImg));
  CHECK_NULL_RET(img, *status, 
		 "memory allocation for SimputImg failed", img);

  // Initialize elements.
  img->naxis1  =0;
  img->naxis2  =0;
  img->dist    =NULL;
  img->fluxscal=0.;
  img->fileref =NULL;
  img->wcs     =NULL;

  return(img);
}


void freeSimputImg(SimputImg** const img)
{
  if (NULL!=*img) {
    if (NULL!=(*img)->dist) {
      if ((*img)->naxis1>0) {
	long ii;
	for (ii=0; ii<(*img)->naxis1; ii++) {
	  if (NULL!=(*img)->dist[ii]) {
	    free((*img)->dist[ii]);
	  }
	}
      }
      free((*img)->dist);
    }
    if (NULL!=(*img)->fileref) {
      free((*img)->fileref);
    }
    if (NULL!=(*img)->wcs) {
      wcsfree((*img)->wcs);
    }
    free(*img);
    *img=NULL;
  }
}


SimputImg* loadSimputImg(const char* const filename, int* const status)
{
  // Image input buffer.
  double* image1d=NULL;

  // String buffer for FITS header.
  char* headerstr=NULL;

  // Get an empty SimputImg data structure.
  SimputImg* img = getSimputImg(status);
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
    fits_read_key(fptr, TFLOAT, "FLUXSCAL", &img->fluxscal, comment, &opt_status);
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


    // Allocate memory for the 1-dimensional image buffer (required for
    // output to FITS file).
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
    fits_write_key(fptr, TSTRING, "HDUCLASS", "HEASARC", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS1", "SIMPUT", "", status);
    fits_write_key(fptr, TSTRING, "HDUCLAS2", "IMAGE", "", status);
    fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0", "", status);
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
      strbuffer[80] = '\0';
      fits_write_record(fptr, strbuffer, status);
      CHECK_STATUS_BREAK(*status);
      strptr += 80;
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


static SimputImg* returnSimputImg(const SimputSource* const src,
				  int* const status)
{
  const int maximgs=200;
  static int nimgs =0;
  static SimputImg** imgs=NULL;

  // Check, whether the source refers to an image.
  if (NULL==src->image) {
    return(NULL);
  }
  if (0==strlen(src->image)) {
    return(NULL);
  }
  if (0==strcmp(src->image, "NULL")) {
    return(NULL);
  }

  // In case there are no images available at all, allocate 
  // memory for the array (storage for images).
  if (NULL==imgs) {
    imgs = (SimputImg**)malloc(maximgs*sizeof(SimputImg*));
    CHECK_NULL_RET(imgs, *status, 
		   "memory allocation for images failed", NULL);
  }

  // Search if the requested image is available in the storage.
  int ii;
  for (ii=0; ii<nimgs; ii++) {
    // Check if the image is equivalent to the requested one.
    if (0==strcmp(imgs[ii]->fileref, src->image)) {
      // If yes, return the image.
      return(imgs[ii]);
    }
  }

  // The requested image is not contained in the storage.
  // Therefore we must load it from the specified location.
  if (nimgs>=maximgs) {
    SIMPUT_ERROR("too many images in the internal storage");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // Check if the filename and filepath properties for the
  // source are properly set.
  CHECK_NULL_RET(src->filepath, *status,
		 "filepath of source catalog is not known", NULL);
  CHECK_NULL_RET(src->filename, *status,
		 "filename of source catalog is not known", NULL);

  // Load the image.
  char filename[SIMPUT_MAXSTR];
  if ('['==src->image[0]) {
    strcpy(filename, *src->filepath);
    strcat(filename, *src->filename);
    strcat(filename, src->image);
  } else {
    if ('/'!=src->image[0]) {
      strcpy(filename, *src->filepath);
    } else {
      strcpy(filename, "");
    }
    strcat(filename, src->image);
  }
  imgs[nimgs]=loadSimputImg(filename, status);
  CHECK_STATUS_RET(*status, imgs[nimgs]);
  nimgs++;

  // Store the file reference to the image for later comparisons.
  imgs[nimgs-1]->fileref = 
    (char*)malloc((strlen(src->image)+1)*sizeof(char));
  CHECK_NULL_RET(imgs[nimgs-1]->fileref, *status, 
		 "memory allocation for file reference failed", 
		 imgs[nimgs-1]);
  strcpy(imgs[nimgs-1]->fileref, src->image);
   
  return(imgs[nimgs-1]);
}


static void p2s(struct wcsprm* const wcs,
		const double px, const double py,
		double* sx, double* sy,
		int* const status)
{
  double pixcrd[2]={ px, py };
  double imgcrd[2], world[2];
  double phi, theta;
  
  // If CUNIT is set to 'degree', change this to 'deg'.
  // Otherwise the WCSlib will not work properly.
  if (0==strcmp(wcs->cunit[0], "degree  ")) {
    strcpy(wcs->cunit[0], "deg");
  }
  if (0==strcmp(wcs->cunit[1], "degree  ")) {
    strcpy(wcs->cunit[1], "deg");
  }

  // Perform the transform using WCSlib.
  int retval=wcsp2s(wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, status);
  CHECK_STATUS_VOID(*status);
  if (0!=retval) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("WCS transformation failed");
    return;
  }
  
  // Convert to [rad].
  *sx=world[0]*M_PI/180.;
  *sy=world[1]*M_PI/180.;
}


void getSimputPhotonCoord(const SimputSource* const src,
			  double* const ra, double* const dec,
			  int* const status)
{
  struct wcsprm wcs={ .flag=-1 };

  do { // Error handling loop.

    // Determine the image.
    SimputImg* img=returnSimputImg(src, status);
    CHECK_STATUS_BREAK(*status);

    // Check, whether the source is point-like or spatially extended.
    if (NULL==img) {
      // The source is point-like.
      *ra =src->ra;
      *dec=src->dec;

    } else {
      // The source is spatially extended.

      // Perform a binary search in 2 dimensions.
      double rnd=static_rndgen() * img->dist[img->naxis1-1][img->naxis2-1];

      // Perform a binary search to obtain the x-coordinate.
      long high = img->naxis1-1;
      long xl   = 0;
      long mid;
      long ymax = img->naxis2-1;
      while (high > xl) {
	mid = (xl+high)/2;
	if (img->dist[mid][ymax] < rnd) {
	  xl = mid+1;
	} else {
	  high = mid;
	}
      }

      // Search for the y coordinate.
      high = img->naxis2-1;
      long yl = 0;
      while (high > yl) {
	mid = (yl+high)/2;
	if (img->dist[xl][mid] < rnd) {
	  yl = mid+1;
	} else {
	  high = mid;
	}
      }
      // Now xl and yl have pixel positions [long pixel coordinates].

      // Create a temporary wcsprm data structure, which can be modified
      // to fit this particular source. The wcsprm data structure contained 
      // in the image should not be modified, since it is used for all 
      // sources including the image.
      wcscopy(1, img->wcs, &wcs);

      // Set the position to the origin and assign the correct scaling.
      wcs.crval[0] = src->ra *180./M_PI; // Units (CUNITn) must be [deg]!
      wcs.crval[1] = src->dec*180./M_PI;
      wcs.cdelt[0] *= 1./src->imgscal;
      wcs.cdelt[1] *= 1./src->imgscal;
      wcs.flag = 0;

      // Check that CUNIT is set to "deg". Otherwise there will be a conflict
      // between CRVAL [deg] and CDELT [different unit]. 
      // TODO This is not required by the standard.
      if (((0!=strcmp(wcs.cunit[0], "deg     ")) && 
	   (0!=strcmp(wcs.cunit[0], "degree  "))) || 
	  ((0!=strcmp(wcs.cunit[1], "deg     ")) &&
	   (0!=strcmp(wcs.cunit[1], "degree  ")))) {
	*status=EXIT_FAILURE;
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "units of image coordinates are '%s' and '%s' (must be 'deg')",
		wcs.cunit[0], wcs.cunit[1]);
	SIMPUT_ERROR(msg);
	break;
      }

      // Determine floating point pixel positions shifted by 0.5 in 
      // order to match the FITS conventions and with a randomization
      // over the pixels.
      double xd = (double)xl + 0.5 + static_rndgen();
      double yd = (double)yl + 0.5 + static_rndgen();

      // Rotate the image (pixel coordinates) by IMGROTA around the 
      // reference point.
      double xdrot =  
	(xd-wcs.crpix[0])*cos(src->imgrota) + 
	(yd-wcs.crpix[1])*sin(src->imgrota) + wcs.crpix[0];
      double ydrot = 
	-(xd-wcs.crpix[0])*sin(src->imgrota) + 
	 (yd-wcs.crpix[1])*cos(src->imgrota) + wcs.crpix[1];
      
      // Convert the long-valued pixel coordinates to double values,
      // including a randomization over the pixel and transform from 
      // pixel coordinates to RA and DEC ([rad]) using the  WCS information.
      p2s(&wcs, xdrot, ydrot, ra, dec, status);
      CHECK_STATUS_BREAK(*status);

      // Determine the RA in the interval from [0:2pi).
      while(*ra>=2.*M_PI) {
	*ra -= 2.*M_PI;
      }
      while(*ra<0.) {
	*ra += 2.*M_PI;
      }

    }
  } while(0); // END of error handling loop.

  // Release memory.
  wcsfree(&wcs);
}


float getSimputSourceExtension(const SimputSource* const src,
			       int* const status)
{
  // Return value.
  float extension=0.;

  struct wcsprm wcs = { .flag=-1 };

  do { // Error handling loop.

    // Get the source image for this particular source.
    SimputImg* img=returnSimputImg(src, status);
    CHECK_STATUS_BREAK(*status);

    // Check if it is a point-like or an extended source.
    if (NULL==img) {
      // Point-like source.
      extension=0.;
      break;

    } else {
      // Extended source => determine the maximum extension.
      double maxext=0.;

      // Copy the wcsprm structure and change the size 
      // according to IMGSCAL.
      wcscopy(1, img->wcs, &wcs);

      // Change the scale of the image according to the source specific
      // IMGSCAL property.
      wcs.crval[0] = 0.; 
      wcs.crval[1] = 0.;
      wcs.cdelt[0] *= 1./src->imgscal;
      wcs.cdelt[1] *= 1./src->imgscal;
      wcs.flag = 0;

      // Check lower left corner.
      double px = 0.5;
      double py = 0.5;
      double sx, sy;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      // TODO This has not been tested extensively, yet.
      while(sx>M_PI) {
	sx-=2*M_PI;
      }
      double ext = sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext = ext;
      }

      // Check lower right corner.
      px = img->naxis1*1. + 0.5;
      py = 0.5;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      while(sx>M_PI) {
	sx-=2*M_PI;
      }
      ext = sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext = ext;
      }
      
      // Check upper left corner.
      px = 0.5;
      py = img->naxis2*1. + 0.5;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      while(sx>M_PI) {
	sx-=2*M_PI;
      }
      ext = sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext = ext;
      }
      
      // Check upper right corner.
      px=img->naxis1*1. + 0.5;
      py=img->naxis2*1. + 0.5;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      while(sx>M_PI) {
	sx-=2*M_PI;
      }
      ext=sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext=ext;
      }
      
      extension=(float)maxext;
    }
  } while(0); // END of error handling loop.

  // Release memory.
  wcsfree(&wcs);

  return(extension);
}


SimputPSD* getSimputPSD(int* const status)
{
  SimputPSD* psd=(SimputPSD*)malloc(sizeof(SimputPSD));
  CHECK_NULL_RET(psd, *status, 
		 "memory allocation for SimputPSD failed", psd);

  // Initialize elements.
  psd->nentries =0;
  psd->frequency=NULL;
  psd->power    =NULL;
  psd->fileref  =NULL;

  return(psd);
}


void freeSimputPSD(SimputPSD** const psd)
{
  if (NULL!=*psd) {
    if (NULL!=(*psd)->frequency) {
      free((*psd)->frequency);
    }
    if (NULL!=(*psd)->power) {
      free((*psd)->power);
    }
    if (NULL!=(*psd)->fileref) {
      free((*psd)->fileref);
    }
    free(*psd);
    *psd=NULL;
  }
}


SimputPSD* loadSimputPSD(const char* const filename, int* const status)
{
  // Get an empty SimputPSD data structure.
  SimputPSD* psd = getSimputPSD(status);
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
    fits_get_colnum(fptr, CASEINSEN, "FREQUENC", &cfrequency, status);
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
      SIMPUT_ERROR("unknown units in FREQUENC column");
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
    printf("PSD '%s' contains %ld data points\n", 
	   filename, psd->nentries);

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




