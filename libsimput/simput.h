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


/** Single Entry in the SimputSourceCatalog. */
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


/** SIMPUT source catalog. */
typedef struct {
  /** Number of entries in the source catalog. */
  int nentries;

  /** Array of the individual entries in the catalog. */
  SimputSourceEntry** entries;

  /** File name (without path contributions) of the FITS file
      containing the source catalog. */
  char* filename;

  /** Path to the FITS file containing the source catalog. */
  char* filepath;

} SimputSourceCatalog;


/** Mission-independent spectrum. */
typedef struct {
  /** Number of entries in the spectrum. */
  int nentries;
  
  /** Energy values [keV]. */
  float* energy;

  /** Source flux distribution [photons/s/cm**2/keV]. */
  float* flux;

} SimputMissionIndepSpec;


/** Spectral probability distribution. */
typedef struct {
  /** Number of entries in the spectrum. */
  int nentries;
  
  /** Energy values [keV]. */
  float* energy;

  /** Probability distribution normalized to the total photon rate
      [photons/s]. */
  float* distribution;

  /** Reference to the location of the spectrum given by the extended
      filename syntax. */
  char* fileref;

} SimputSpectralDistribution;


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



/** Constructor for the SimputMissionIndepSpec data
    structure. Allocates memory, initializes elements with their
    default values and pointers with NULL. */
SimputMissionIndepSpec* getSimputMissionIndepSpec(int* const status);

/** Destructor for the SimputSourceEntry. Calls destructor routines
    for all contained elements, releases the allocated memory, and
    finally sets the pointer to NULL. */
void freeSimputMissionIndepSpec(SimputMissionIndepSpec** const spec);

/** Load the SimputMissionIndepSpec from the specified file. */
SimputMissionIndepSpec* loadSimputMissionIndepSpec(const char* const filename,
						   int* const status);


/** Set the instrument ARF containing the effective area. This
    information is required to obtain a mission-specific spectrum from
    the mission-independent format. */
void simputSetARF(struct ARF* const arf);


/** Return a random photon energy according to the distribution
    defined by the energy spectrum of the particular source. If the
    spectrum is not stored in memory, its loaded from the location
    specified in the the catalog. */
float getSimputPhotonEnergy(const SimputSourceEntry* const src,
			    const SimputSourceCatalog* const cat,
			    int* const status);

/** Constructor for the SimputSpectralDistribution data
    structure. Allocates memory, initializes elements with their
    default values and pointers with NULL. */
SimputSpectralDistribution* getSimputSpectralDistribution(int* const status);

/** Convolve the given mission-independent spectrum with the
    instrument ARF in order to obtain the spectral probability
    distribution. */
SimputSpectralDistribution* 
convSimputMissionIndepSpecWithARF(const SimputMissionIndepSpec* const indepspec, 
				  int* const status);


#endif /* SIMPUT_H */
