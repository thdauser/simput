#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"
#include "heasp.h"

/*
  simput.h

  Version  Date       Author                
  --------------------------------------------------------------------------- 
  0.0.1  2011/04/02   Christian Schmid      initial 
  --------------------------------------------------------------------------- 

   This is the header file of the simput library. The library provides
   basic routines to create a SIMPUT source catalog file according to
   the standard defined in the SIMPUT formats document (reference ???).
*/


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  /** Unique source ID. */
  int src_id;

  /** Source name. */
  char* src_name;

  /** Right ascension of source position [rad]. */
  double ra;
  
  /** Declination of source position [rad]. */
  double dec;

  /** Image and polarization rotation angle [rad]. */
  float imgrota;

  /** Image scaling factor. Default value is 1. May not be 0. */
  float imgscal;

  /** Lower limit of reference energy band [keV]. */
  float e_min;

  /** Upper limit of reference energy band [keV]. */
  float e_max;

  /** Source flux density in the reference energy band
      [erg/s/cm^2]. */
  float flux;

  /** Reference to the storage location of the source spectrum. */
  char* spectrum;

  /** Reference to the storage location of the source image. */
  char* image;

  /** Reference to the storage location of the source light curve. */
  char* lightcur;

} SimputSourceEntry;


typedef struct {
  /** Number of entries in the source catalog. */
  int nentries;

  /** Array of the individual entries in the catalog. */
  SimputSourceEntry** entries;

  // TODO
  /** File name (without path contributions) of the FITS file
      containing the source catalog. */
  char* filename;

  // TODO
  /** Path to the FITS file containing the source catalog. */
  char* path;

} SimputSourceCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for the SimputSourceEntry data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputSourceEntry* getSimputSourceEntry(int* const status);

/** Constructor for the SimputSourceEntry data structure. Allocates
    memory and initializes elements with the given values. */
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
					 int* const status);

/** Destructor for the SimputSourceEntry. Calls destructor routines
    for all contained elements, releases the allocated memory, and
    finally sets the pointer to NULL. */
void freeSimputSourceEntry(SimputSourceEntry** const entry);

/** Constructor for the SimputSourceCatalog data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputSourceCatalog* getSimputSourceCatalog(int* const status);

/** Destructor for the SimputSourceCatalog. Calls destructor routines
    for all contained SimputSourceEntry elements, releases the
    allocated memory, and finally sets the pointer to NULL. */
void freeSimputSourceCatalog(SimputSourceCatalog** const catalog);

/** Load the SIMPUT source catalog from the specified file. */
SimputSourceCatalog* loadSimputSourceCatalog(const char* const filename,
					     int* const status);

/** Store the SimputSourceCatalog in the specified file. */
void saveSimputSourceCatalog(const SimputSourceCatalog* const catalog,
			     const char* const filename,
			     int* const status);


/** Set the telescope ARF containing the effective area. This
    information is required to obtain a mission-specific spectrum from
    the mission-independent format. */
void simputSetARF(struct ARF* const arf);


#endif /* SIMPUT_H */
