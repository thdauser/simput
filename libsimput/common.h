#ifndef COMMON_H
#define COMMON_H 1


/////////////////////////////////////////////////////////////////
// Includes.
/////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "simput.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/** Common string length. */
#define SIMPUT_MAXSTR (1025)


/** Chatter level:

    0: error messages only

    1: error messages and warnings

    2: error messages, warnings, and informational output. 

    The chatter level can be defined as a compiler flag. If it is not
    set externally, use the default value of 1. */
#ifndef DCHATTY
#define DCHATTY 1
#endif


/** Different types of extensions in SIMPUT FITS files. */
#define EXTTYPE_NONE (0)
#define EXTTYPE_MIDPSPEC (1)
#define EXTTYPE_IMAGE (2)
#define EXTTYPE_PHOLIST (3)
#define EXTTYPE_LC (4)
#define EXTTYPE_PSD (5)


/** Maximum number of mission-independet spectra in the cache. */
#define MAXMIDPSPEC (30000) 


/////////////////////////////////////////////////////////////////
// Macros.
/////////////////////////////////////////////////////////////////


/** Output routine for error messages. */
#define SIMPUT_ERROR(msg) \
  fprintf(stderr, "\nError in %s: %s!\n", __func__, msg)

/** Chatter routine for warnings printed to STDOUT. */
#define SIMPUT_WARNING(msg) \
  if(1<=DCHATTY) {printf("\n*** Warning in %s: %s! ***\n", __func__, msg);}

/** Chatter routine for informational output to STDOUT. */
#define SIMPUT_INFO(msg) \
  if(2==DCHATTY) {printf("%s\n", msg);}


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
// Data structures.
/////////////////////////////////////////////////////////////////


struct SimputSrcBuffer {
  long nsrcs; // Current number of sources in the cache.
  long csrc;  // Index of next position in cache that will be used.
  SimputSrc** srcs; // Cache for the sources.

  // This array contains the row numbers of the sources in the storage
  // given in the same order as the corresponding sources in the storage.
  // The array is used to replace the oldest source in the storage.
  long* rownums;

  // Array with a size corresponding to the number of entries 
  // in the catalog. Each entry in the array refers to the index
  // of the corresponding source in the storage. If the respective
  // source is not contained in the storage, the value in the array
  // is -1.
  long* rowmap;
};


struct SimputExttypeBuffer {
  long nhdus;// Current number of extensions in the cache.
  long chdu; // Index of next position in the cache that will be used.
  int* hdus; // Cache for the extension types.
  char** filenames; // Cache for the corresponding file references.
};


struct SimputSpecBuffer {
  long nspectra;  // Current number of spectra in the cache.
  long cspectrum; // Index of next position in the cache that will be used.
  SimputSpec** spectra; // Cache for the spectra.
};


struct SimputMIdpSpecBuffer {
  long nspectra;  // Current number of spectra in the cache.
  long cspectrum; // Index of next position in the cache that will be used.
  SimputMIdpSpec** spectra; // Cache for the spectra.
};


struct SimputLCBuffer {
  long nlcs; // Current number of light curves in the cache.
  long clc;  // Index of next position in the cache that will be used.
  SimputLC** lcs; // Cache for the light curves.
};


struct SimputKRLCBuffer {
  long nkrlcs; // Current number of K&R light curves in the cache.
  long ckrlc;  // Index of next position in the cache that will be used.
  SimputKRLC** krlcs; // Cache for the K&R light curves.
};


struct SimputPSDBuffer {
  long npsds; // Current number of PSDs in the cache.
  SimputPSD** psds; // Cache for the PSDs.
};


struct SimputImgBuffer {
  long nimgs; // Current number of images in the cache.
  SimputImg** imgs; // Cache for the images.
};


/////////////////////////////////////////////////////////////////
// Functions.
/////////////////////////////////////////////////////////////////


struct SimputSrcBuffer* newSimputSrcBuffer(int* const status);
void freeSimputSrcBuffer(struct SimputSrcBuffer** sb);


struct SimputExttypeBuffer* newSimputExttypeBuffer(int* const status);
void freeSimputExttypeBuffer(struct SimputExttypeBuffer** eb);


struct SimputMIdpSpecBuffer* newSimputMIdpSpecBuffer(int* const status);
void freeSimputMIdpSpecBuffer(struct SimputMIdpSpecBuffer** sb);


SimputSpec* newSimputSpec(int* const status);
void freeSimputSpec(SimputSpec** const spec);


struct SimputSpecBuffer* newSimputSpecBuffer(int* const status);
void freeSimputSpecBuffer(struct SimputSpecBuffer** sb);


struct SimputLCBuffer* newSimputLCBuffer(int* const status);
void freeSimputLCBuffer(struct SimputLCBuffer** sb);


SimputKRLC* newSimputKRLC(int* const status);
void freeSimputKRLC(SimputKRLC** const lc);


struct SimputKRLCBuffer* newSimputKRLCBuffer(int* const status);
void freeSimputKRLCBuffer(struct SimputKRLCBuffer** sb);


struct SimputPSDBuffer* newSimputPSDBuffer(int* const status);
void freeSimputPSDBuffer(struct SimputPSDBuffer** sb);


struct SimputImgBuffer* newSimputImgBuffer(int* const status);
void freeSimputImgBuffer(struct SimputImgBuffer** sb);


/** Determine the extension type of a particular FITS file HDU. */
int getExtType(SimputCtlg* const cat, 
	       char* const filename, 
	       int* const status);


#endif /* COMMON_H */
