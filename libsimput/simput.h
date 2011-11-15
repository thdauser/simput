#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"
#include <wcslib/wcshdr.h>

// GSL header files
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_halfcomplex.h>

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/*
  simput.h

  Version  Date       Author                
  --------------------------------------------------------------------------- 
  0.0.1  2011/04/02   Christian Schmid      initial 
  0.0.2  2011/09/08   Christian Schmid      updated
  --------------------------------------------------------------------------- 

   This is the header file of the simput library. The library provides
   basic routines to create a SIMPUT source catalog file according to
   the standard defined in the SIMPUT formats document (reference ???).
*/


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Single entry in the SimputSourceCatalog. Requires about 128 bytes,
    depending on the string lengths. */
typedef struct {
  /** Unique source ID. */
  long src_id;

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

  /** Source energy flux density in the reference energy band
      [erg/s/cm^2]. */
  float eflux;

  /** Reference to the storage location of the source spectrum. */
  char* spectrum;

  /** Reference to the storage location of the source image. */
  char* image;

  /** Reference to the storage location of the source light curve. */
  char* lightcur;

  /** Pointer to the filename in the source catalog file. */
  char** filename;

  /** Pointer to the filepath in the source catalog file. */
  char** filepath;

} SimputSource;


typedef struct {
  /** Pointer to the catalog extension. */
  fitsfile* fptr;

  /** Total number of entries in the catalog. */
  long nentries;

  /** Column numbers. */
  int csrc_id, csrc_name, cra, cdec, cimgrota, cimgscal,
    ce_min, ce_max, cflux, cspectrum, cimage, clightcur;

  /** Unit conversion factors. */
  float fra, fdec, fimgrota, fe_min, fe_max, fflux;

  /** File name (without path contributions) of the FITS file
      containing the source catalog. This value is automatically set
      when a catalog is loaded from or saved to a file. This pointer
      should not be modified directly. */
  char* filename;

  /** Path to the FITS file containing the source catalog. This value
      is automatically set when a catalog is loaded from or saved to a
      file. This pointer should not be modified directly. */
  char* filepath;

  /** Buffer for pre-loaded SIMPUT sources. */
  void* srcbuff;

} SimputCatalog;


/** Reference flux of a certain spectrum in a particular energy
    band. */
typedef struct {
  /** Lower Boundary of the energy band [keV]. */
  float emin;

  /** Upper Boundary of the energy bin [keV]. */
  float emax;

  /** Reference flux in this band [erg/s/cm^2]. */
  float flux;

} SimputSpecBandFlux;


/** Mission-independent spectrum. */
typedef struct {
  /** Number of entries in the spectrum. */
  long nentries;
  
  /** Energy values [keV]. */
  float* energy;

  /** Photon flux distribution [photons/s/cm**2/keV]. */
  float* pflux;

  /** Energy flux in a certain reference band. This value is only
      calculated once (if the same reference band is used all the
      time) and stored here in order to avoid recalculation during the
      photon generation. */
  SimputSpecBandFlux* refflux;

  /** Probability distribution normalized to the total photon rate
      [photons/s]. */
  double* distribution;

  /** Unique case-sensitive designator for an individual spectrum. */
  char* name;

  /** Reference to the location of the spectrum given by the extended
      filename syntax. This reference is used to check, whether a
      spectrum is already contained in the internal storage. */
  char* fileref;

} SimputMIdpSpec;


/** SIMPUT light curve. */
typedef struct {
  /** Number of entries in the light curve. */
  long nentries;

  /** Time values [s]. */
  double* time;

  /** Phase values (between 0 and 1). */
  float* phase;

  /** Relative flux values (unitless). */
  float* flux;

  /** Piece-wise linear light curve data points. The value a_k
      represents the gradient of the light curve between the time t_k
      and t_{k+1} (slope, [1/s]). The value b_k represents the
      constant contribution (intercept) at t_k. These values include
      the FLUXSCAL. */
  float *a, *b;

  /** Reference to the storage location of the source spectrum at a
      particular point of time or phase respectively. */
  char** spectrum;

  /** Reference to the storage location of the source image at a
      particular point of time or phase respectively. */
  char** image;

  /** MJD for reference time [d]. */
  double mjdref;

  /** Zero time [s]. */
  double timezero;

  /** Phase of periodic oscillation at timezero. */
  float phase0;

  /** Duration of one oscillation period [s]. */
  float period;
  
  /** Flux scaling factor. */
  float fluxscal;

  /** If the light curve has been produced from a PSD, it is assigned
      to a particular source and cannot be re-used for different
      sources. In that case the SRC_ID of the respective source is
      stored in this variable. Otherwise its value is set to 0. */
  long src_id;

  /** Reference to the location of the light curve given by the
      extended filename syntax. This reference is used to check,
      whether the light curve is already contained in the internal
      storage. */
  char* fileref;

} SimputLC;


/** SIMPUT power spectral density (PSD). */
typedef struct {
  /** Number of entries in the PSD. */
  long nentries;

  /** Frequency [Hz]. */
  float* frequency;

  /** Power spectral density with Miyamoto normalization [Hz^-1]. */
  float* power;

  /** Reference to the location of the light curve given by the
      extended filename syntax. This reference is used to check,
      whether the light curve is already contained in the internal
      storage. */
  char* fileref;

} SimputPSD;


/** SIMPUT source image. */
typedef struct {
  /** Image dimensions. */
  long naxis1, naxis2;

  /** Pixel value distribution function. */
  double** dist;

  /** Flux scaling factor. */
  float fluxscal;

  /** WCS data used by wcslib. */
  struct wcsprm* wcs;

  /** Reference to the location of the source image given by the
      extended filename syntax. This reference is used to check,
      whether the source image is already contained in the internal
      storage. */
  char* fileref;

} SimputImg;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for the SimputCatalog data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputCatalog* getSimputCatalog(int* const status);

/** Destructor for the SimputCatalog. Closes the FITS file,
    releases the allocated memory, and finally sets the pointer to
    NULL. */
void freeSimputCatalog(SimputCatalog** const catalog,
		       int* const status);

/** Open a FITS file with a SIMPUT source catalog. The access mode can
    be either READONLY to open an existing catalog or READWRITE for
    both existing or non-existing files. */
SimputCatalog* openSimputCatalog(const char* const filename,
				 const int mode,
				 int* const status);

/** Constructor for the SimputSource data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputSource* getSimputSource(int* const status);

/** Constructor for the SimputSource data structure. Allocates
    memory and initializes elements with the given values. */
SimputSource* getSimputSourceV(const long src_id, 
			       const char* const src_name,
			       /** ([rad]). */
			       const double ra,
			       /** ([rad]). */
			       const double dec,
			       /** Image rotation angle ([rad]). Only
				   applicable for extended sources
				   with an image. */
			       const float imgrota,
			       const float imgscal,
			       /** Lower boundary of reference energy
				   band ([keV]). */
			       const float e_min,
			       /** Upper boundary of reference energy
				   band ([keV]). */
			       const float e_max,
			       /** Energy flux density in the
				   reference energy band
				   ([erg/s/cm^2]). */
			       const float eflux,
			       const char* const spectrum,
			       const char* const image,
			       const char* const lightcur,
			       int* const status);

/** Destructor for the SimputSource. Calls destructor routines
    for all contained elements, releases the allocated memory, and
    finally sets the pointer to NULL. */
void freeSimputSource(SimputSource** const entry);

/** Load a source from a particular row of the catalog in the FITS
    file. Row numbering starts at 1. The returned pointer to the
    SimputSource must be free'd afterwards in order to avoid a memory
    leak. */
SimputSource* loadSimputSource(SimputCatalog* const cat,
			       const long row,
			       int* const status);

/** Return an entry from a SimputCatalog, which is contained in a
    particular row of the FITS table. According to the FITS
    conventions row numbering starts with 1 for the first line. When
    loading the SimputSource for the first time, the data are stored
    in an static cache, such that they do not have to be loaded again
    on later access. The returned pointer to the SimputSource should
    not be free'd, since the allocated memory is managed by the
    caching mechanism. */
SimputSource* loadCacheSimputSource(SimputCatalog* const cat,
				    const long row,
				    int* const status);

/** Append a SimputSource to an existing catalog. The source is
    inserted at the end of the binary table in the FITS file. */
void appendSimputSource(SimputCatalog* const cat,
			SimputSource* const src,
			int* const status);

/** Append an array of SimputSources to an existing catalog. The
    sources are inserted at the end of the binary table in the FITS
    file. */
void appendSimputSourceBlock(SimputCatalog* const cat,
			     SimputSource** const src,
			     const long nsources,
			     int* const status);

/** Constructor for the SimputMIdpSpec data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputMIdpSpec* getSimputMIdpSpec(int* const status);

/** Destructor for the SimputSource. Calls destructor routines
    for all contained elements, releases the allocated memory, and
    finally sets the pointer to NULL. */
void freeSimputMIdpSpec(SimputMIdpSpec** const spec);

/** Load the SimputMIdpSpec from the specified file. */
SimputMIdpSpec* loadSimputMIdpSpec(const char* const filename,
				   int* const status);

/** Load the requested spectrum. Keeps a certain number of spectra in
    an internal cache. If the requested spectrum is not located in the
    cache, it is loaded from the given filename. The returned pointer
    to the SimputMIdpSpec should not be free'd, since the
    allocated memory is managed by the internal caching mechanism. */
SimputMIdpSpec* 
loadCacheSimputMIdpSpec(const char* const filename,
			int* const status);

/** Return the spectrum of the specified SimputSource for the
    particular point of time. */
SimputMIdpSpec* returnSimputSrcSpec(const SimputSource* const src,
				    const double time, 
				    const double mjdref,
				    int* const status);

/** Save the mission-independent spectrum in the specified extension
    of the given FITS file. If the file does not exist yet, a new file
    is created. If the file exists, but does not contain the specified
    extension, an appropriate HDU is created. If the extension exists,
    the unambiguousness of the spectrum name (contained in the data
    structure) is verified. */
void saveSimputMIdpSpec(SimputMIdpSpec* const spec,
			const char* const filename,
			char* const extname,
			int extver,
			int* const status);

/** Determine the energy and flux values of a particular bin in the
    SimputMIdpSpec. The energy is given in [keV], the flux in
    [photons/s/cm**2/keV]. */
void getSimputSpectrumValue(const SimputMIdpSpec* const spec,
			    const long row,
			    float* const energy, 
			    float* const pflux,
			    int* const status);

/** Set the instrument ARF containing the effective area. This
    information is required to obtain a mission-specific spectrum from
    the mission-independent format. The access to the ARF data
    structure must be guaranteed as long as the SIMPUT library
    routines are used. */
void simputSetARF(struct ARF* const arf);

/** Set the random number generator, which is used by the simput
    library routines. The generator should return double valued,
    uniformly distributed numbers in the interval [0,1). */
void simputSetRndGen(double(*rndgen)(void));

/** Return a randomized photon energy according to the distribution
    defined by the energy spectrum of the particular source. If the
    spectrum is not stored in memory, it is loaded from the location
    specified in the the catalog. If the source has a time-variable
    spectrum, which is defined via the light curve extension, the
    reference time is required. */
float getSimputPhotonEnergy(const SimputSource* const src,
			    const double time,
			    const double mjdref,
			    int* const status);

/** Determine the energy flux in [erg/s/cm**2] within the reference
    energy band of the specified source valid a the requested point of
    time. */
float getSimputSrcBandFlux(const SimputSource* const src,
			   const double time, const double mjdref,
			   int* const status);

/** Determine the energy flux of the spectrum in [erg/s/cm**2] within a
    certain energy band from emin to emax. */
float getSimputSpecBandFlux(SimputMIdpSpec* const spec,
			    const float emin, const float emax,
			    int* const status);

/** Return the photon rate of a particular source. The return value is
    the nominal photon rate given in the source catalog. WARNING: It
    does not contain any light curve or other time-variable
    contributions. Specification of instrument ARF required. */
float getSimputPhotonRate(const SimputSource* const src,
			  const double time, const double mjdref,
			  int* const status);


/** Constructor for the SimputLC data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputLC* getSimputLC(int* const status);

/** Destructor for the SimputLC. Calls destructor routines for all
    contained elements, releases the allocated memory, and finally
    sets the pointer to NULL. */
void freeSimputLC(SimputLC** const lc);

/** Load the SimputLC from the specified file. */
SimputLC* loadSimputLC(const char* const filename, int* const status);

/** Save the light curve in the specified extension of the given FITS
    file. If the file does not exist yet, a new file is created. If
    the file exists, an appropriate HDU is created. */
void saveSimputLC(SimputLC* const lc, const char* const filename,
		  char* const extname, int extver,
		  int* const status);


/** Constructor for the SimputPSD data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputPSD* getSimputPSD(int* const status);

/** Destructor for the SimputPSD. Calls destructor routines for all
    contained elements, releases the allocated memory, and finally
    sets the pointer to NULL. */
void freeSimputPSD(SimputPSD** const lc);

/** Load the SimputPSD from the specified file. */
SimputPSD* loadSimputPSD(const char* const filename, int* const status);


/** Return a randomized photon arrival time at the telescope according
    to the photon rate and possible light curve. As the light curve is
    time-dependet, a reference time defined by the arrival time of the
    previous photon has to be specified. If the source refers to a
    light curve, which is not stored in memory, it is loaded from the
    location specified in the the catalog. If the light curve time is
    exceeded or the source flux is zero, the 'failed' flag is set to
    1. Otherwise its value is 0. */
double getSimputPhotonTime(const SimputSource* const src,
			   double prevtime,
			   const double mjdref,
			   int* const failed,
			   int* const status);


/** Constructor for the SimputImg data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputImg* getSimputImg(int* const status);

/** Destructor for the SimputImg. Calls destructor routines for all
    contained elements, releases the allocated memory, and finally
    sets the pointer to NULL. */
void freeSimputImg(SimputImg** const img);

/** Load a SIMPUT source image from the specified file and store it in
    a SimputImg data structure. */
SimputImg* loadSimputImg(const char* const filename, int* const status);

/** Save the source image in the specified extension of the given FITS
    file. If the file does not exist yet, a new file is created. If
    the file exists, an appropriate HDU is created. */
void saveSimputImg(SimputImg* const img, const char* const filename,
		   char* const extname, int extver,
		   int* const status);

/** Determine the coordinates (RA and Dec) of a new photon emerging
    from a particular source. For point-like sources the coordinates
    are equivalent with the source position. For extended sources a
    random position is determined according to the flux distribution
    defined by the source image. The returned coordinate values are
    given in [rad]. */
void getSimputPhotonCoord(const SimputSource* const src,
			  double* const ra, double* const dec,
			  int* const status);

/** Return the maximum spatial extension of a particular source around
    its reference point in [rad]. */
float getSimputSourceExtension(const SimputSource* const src,
			       int* const status);


#endif /* SIMPUT_H */
