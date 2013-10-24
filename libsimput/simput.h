#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"
#include "wcs.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/*
  simput.h

  Version  Date       Author                
  -------------------------------------------------------------------- 
  0.1  2011/04/02   Christian Schmid      initial 
  0.2  2011/09/08   Christian Schmid      updated

  0.11 2012/01/26   Christian Schmid      updated
  0.12 2012/02/21   Christian Schmid      updated

  0.25 2012/08/31   Christian Schmid      renewed interface
  -------------------------------------------------------------------- 

   This is the header file of the SIMPUT library. The library provides
   basic routines to create a SIMPUT source catalog file according to
   the standard defined in the SIMPUT format document (reference ???).
*/


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Single entry in the SimputCtlg. Requires about 128 bytes,
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

  /** Flux of the source in the reference energy band [erg/s/cm^2]. */
  float eflux;

  /** Photon rate. Determined from the spectrum and the reference flux. */
  float* phrate;

  /** Reference string to the storage location of the spectrum of the
      source. */
  char* spectrum;

  /** Reference string to the storage location of a spatial image of
      the source. */
  char* image;

  /** Reference string to the storage location of a light curve or PSD
      of the source. */
  char* timing;

} SimputSrc;


/** Main data structure providing access to a SIMPUT catalog in a FITS
    file. The data structure contains basic information about the
    catalog as well as buffers for data extensions. */
typedef struct {

  /** Pointer to the catalog extension. */
  fitsfile* fptr;

  /** Total number of entries in the catalog. */
  long nentries;

  /** Column numbers. */
  int csrc_id, csrc_name, cra, cdec, cimgrota, cimgscal,
    ce_min, ce_max, cflux, cspectrum, cimage, ctiming;

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

  /** Buffer for FITS file HDU types. */
  void* extbuff;

  /** Buffer for pre-loaded mission-independent spectra. */
  void* midpspecbuff;
  
  /** Buffer for pre-loaded photon lists. */
  void* phlistbuff;

  /** Buffer for pre-loaded light curves. */
  void* lcbuff;

  /** Buffer for pre-loaded power spectra. */
  void* psdbuff;

  /** Buffer for pre-loaded images. */
  void* imgbuff;

  /** Buffer for pre-loaded spectra. */
  void* specbuff;

  /** Buffer for pre-loaded Klein and Robert light curves. */
  void* krlcbuff;

  /** Instrument ARF. */
  struct ARF* arf;

} SimputCtlg;


/** Mission-independent spectrum. */
typedef struct {

  /** Number of entries in the spectrum. */
  long nentries;
  
  /** Energy values [keV]. */
  float* energy;

  /** Photon flux distribution [photons/cm**2/keV]. */
  float* pflux;

  /** Unique case-sensitive designator for an individual spectrum. */
  char* name;

  /** Reference to the location of the spectrum given by the extended
      filename syntax. This reference is used to check, whether the
      spectrum is already contained in the internal storage. */
  char* fileref;

} SimputMIdpSpec;


/** Spectrum (spectral distribution function). */
typedef struct {

  /** Probability distribution normalized to the total photon number
      [photons]. */
  double* distribution;

  /** Reference to the location of the spectrum given by the extended
      filename syntax. This reference is used to check, whether the
      spectrum is already contained in the internal storage. */
  char* fileref;

} SimputSpec;


/** SIMPUT light curve. */
typedef struct {

  /** Number of entries in the light curve. */
  long nentries;

  /** Time values [s]. */
  double* time;

  /** Phase values (between 0 and 1). */
  double* phase;

  /** Relative flux values (unitless). This value in combination with
      the flux scaling factor determines the variation of the source
      reference flux given in the catalog table. */
  float* flux;

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

  /** Phase of periodic signal at timezero. */
  double phase0;

  /** Duration of one period [s]. */
  double period;

  /** First derivative of period with respect to time [unitless] */
  double dperiod;

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


/** SIMPUT light curve converted to the form needed by the Klein &
    Roberts (1984) algorithm. */
typedef struct {

  /** Number of entries in the light curve. */
  long nentries;

  /** Time values [s]. */
  double* time;

  /** Phase values (between 0 and 1). */
  double* phase;

  /** Piece-wise linear light curve data points. The value b_k
      represents the constant contribution (intercept) at t_k. It
      includes the FLUXSCAL. */
  double *b;

  /** MJD for reference time [d]. */
  double mjdref;

  /** Zero time [s]. */
  double timezero;

  /** Phase of periodic signal at timezero. */
  double phase0;

  /** Duration of one period [s]. */
  double period;
  
  /** First derivative of period with respect to time [unitless] */
  double dperiod;

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

} SimputKRLC;


/** SIMPUT power spectral density (PSD). */
typedef struct {

  /** Number of entries in the PSD. */
  long nentries;

  /** Frequency [Hz]. */
  float* frequency;

  /** Power spectral density with Miyamoto normalization [Hz^-1]. */
  float* power;

  /** Reference to the location of the power spectrum given by the
      extended filename syntax. This reference is used to check,
      whether the power spectrum is already contained in the internal
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


/** SIMPUT photon list. */
typedef struct {

  /** Pointer to the FITS file HDU. */
  fitsfile* fptr;

  /** Column numbers. */
  int cra, cdec, cenergy, ctime;

  /** Unit conversion factors. */
  float fra, fdec, fenergy, ftime;

  /** Total number of photons in the list. */
  long nphs;

  /** MJD for reference time [d]. */
  double mjdref;

  /** Zero time [s]. */
  double timezero;

  /** Start and stop time according to the FITS header keywords. */
  double tstart, tstop;

  /** Reference effective area [cm^2]. This value must be greater or
      equal to the maximum value of the used ARF. */
  float refarea;

  /** Acceptance rate for the photons in the list. This value is
      needed if the time information is used. */
  double accrate;

  /** Current row in the photon list. Numbering starts a 1 for the
      first row. This variable is needed if the time information in
      the photon list is taken into account. */
  long currrow;

  /** Reference to the location of the photon list given by the
      extended filename syntax. This reference is used to check,
      whether the photon list is already contained in the internal
      storage. */
  char* fileref;

} SimputPhList;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for the SimputCtlg data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputCtlg* newSimputCtlg(int* const status);

/** Destructor for the SimputCtlg. Closes the FITS file, releases the
    allocated memory, and finally sets the pointer to NULL. */
void freeSimputCtlg(SimputCtlg** const ctlg, int* const status);

/** Open a FITS file with a SIMPUT source catalog. The access mode can
    be either READONLY to open an existing catalog or READWRITE for
    both existing or non-existing files. The maximum length of strings
    to be inserted in the string columns can be specified. These
    values are only required for the creation of a new catalog. If the
    values are set to 0, variable length string columns will be
    used. Otherwise fixed length string columns are used. */
SimputCtlg* openSimputCtlg(const char* const filename,
			   const int mode, 
			   const int maxstrlen_src_name,
			   const int maxstrlen_spectrum,
			   const int maxstrlen_image,
			   const int maxstrlen_timing,
			   int* const status);

/** Returns the number of sources contained in the specified
    SIMPUTCtlg. */
long getSimputCtlgNSources(const SimputCtlg* const cat);


/** Constructor for the SimputSrc data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputSrc* newSimputSrc(int* const status);

/** Constructor for the SimputSrc data structure. Allocates memory and
    initializes elements with the given values. */
SimputSrc* newSimputSrcV(const long src_id, 
			 const char* const src_name,
			 /** ([rad]). */
			 const double ra,
			 /** ([rad]). */
			 const double dec,
			 /** Image rotation angle ([rad]). Only
			     applicable for extended sources with an
			     image. */
			 const float imgrota,
			 const float imgscal,
			 /** Lower boundary of reference energy band
			     ([keV]). */
			 const float e_min,
			 /** Upper boundary of reference energy band
			     ([keV]). */
			 const float e_max,
			 /** Energy flux density in the reference
			     energy band ([erg/s/cm^2]). */
			 const float eflux,
			 const char* const spectrum,
			 const char* const image,
			 const char* const timing,
			 int* const status);

/** Destructor for the SimputSrc. Calls destructor routines
    for all contained elements, releases the allocated memory, and
    finally sets the pointer to NULL. */
void freeSimputSrc(SimputSrc** const src);

/** Load a source from a particular row of the catalog in the FITS
    file. Row numbering starts at 1. The returned pointer to the
    SimputSrc must be free'd afterwards in order to avoid a memory
    leak. */
SimputSrc* loadSimputSrc(SimputCtlg* const cat,
			 const long row,
			 int* const status);


/** Return an entry from a SimputCtlg, which is contained in a
    particular row of the FITS table. According to the FITS
    conventions row numbering starts with 1 for the first line. When
    loading the SimputSrc for the first time, the data are stored
    in a static cache, such that they do not have to be loaded again
    on later access. The returned pointer to the SimputSrc must
    not be free'd, since the allocated memory is managed by the
    caching mechanism. */
SimputSrc* getSimputSrc(SimputCtlg* const cat,
			const long row,
			int* const status);

/** Append a SimputSrc to an existing catalog. The source is
    inserted at the end of the binary table in the FITS file. */
void appendSimputSrc(SimputCtlg* const cat,
		     SimputSrc* const src,
		     int* const status);

/** Append an array of SimputSrc data structures to an existing
    catalog. The sources are inserted at the end of the binary table
    in the FITS file. */
void appendSimputSrcBlock(SimputCtlg* const cat,
			  SimputSrc** const src,
			  const long nsources,
			  int* const status);

/** Determine the reference of a source to its spectrum. The spectrum
    can be either referred to directly in the source catalog or in the
    optional light curve extension. In the latter case the point of
    time for which the spectrum is requested, needs to be specified. A
    negative value of mjdref indicates that the spectrum referred to
    by the first bin in the light curve should be returned. */
void getSimputSrcSpecRef(SimputCtlg* const cat,
			 const SimputSrc* const src,
			 const double prevtime,
			 const double mjdref,
			 char* const specref,
			 int* const status);

/** Constructor for the SimputMIdpSpec data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputMIdpSpec* newSimputMIdpSpec(int* const status);

/** Destructor for the SimputMIdpSpec. Calls destructor routines for
    all contained elements, releases the allocated memory, and finally
    sets the pointer to NULL. */
void freeSimputMIdpSpec(SimputMIdpSpec** const spec);

/** Load the SimputMIdpSpec from the specified file. */
SimputMIdpSpec* loadSimputMIdpSpec(const char* const filename,
				   int* const status);

/** Loads all spectra from the specified FITS binary table into the
    internal cache. */
void loadCacheAllSimputMIdpSpec(SimputCtlg* const cat,
				const char* const filename,
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

/** Save the mission-independent spectrum in the specified extension
    of the given FITS file. If the file does not exist yet, a new file
    is created. If the file exists, but does not contain the specified
    extension, an appropriate HDU is created. */
void saveSimputMIdpSpecBlock(SimputMIdpSpec** const spec,
			     const long nspec,
			     const char* const filename,
			     char* const extname,
			     int extver,
			     int* const status);

/** Return the spectrum of the specified SimputSrc for the
    particular point of time. If the required spectrum is already
    located in the internal cache, a pointer to this spectrum will be
    returned. Otherwise the spectrum is loaded into the cache. The
    pointer returned by the function must no be free'd. */
SimputMIdpSpec* getSimputSrcMIdpSpec(SimputCtlg* const cat,
				     const SimputSrc* const src,
				     const double prevtime,
				     const double mjdref,
				     int* const status);

/** Determine the energy and flux values of a particular bin in the
    SimputMIdpSpec. The energy is given in [keV], the flux in
    [photons/cm**2/keV]. */
void getSimputMIdpSpecVal(const SimputMIdpSpec* const spec,
			  const long row,
			  float* const energy, 
			  float* const pflux,
			  int* const status);

/** Determine the flux in the specified energy band in
    [erg/cm**2]. */
float getSimputMIdpSpecBandFlux(SimputMIdpSpec* const spec,
				const float emin, 
				const float emax);

/** Set the instrument ARF containing the effective area. This
    information is required to obtain a mission-specific spectrum from
    the mission-independent format. The access to the ARF data
    structure must be guaranteed as long as the SIMPUT library
    routines are used. */
void setSimputARF(SimputCtlg* const cat, struct ARF* const arf);

/** Load and set the instrument ARF containing the effective
    area. This information is required to obtain a mission-specific
    spectrum from the mission-independent format. If an error occurs,
    the status variable is set to EXIT_FAILURE. */
void loadSimputARF(SimputCtlg* const cat, 
		   char* const filename, 
		   int* const status);

/** Set the random number generator, which is used by the simput
    library routines. The generator should return double valued,
    uniformly distributed numbers in the interval [0,1). */
void setSimputRndGen(double(*rndgen)(int* const));

/** Return the photon rate of a particular source. The return value is
    the nominal photon rate for the whole spectrum according to the
    reference flux given in the source catalog. WARNING: It does not
    contain any light curve or other time-variable contributions. A
    specification of an instrument ARF required. */
float getSimputPhotonRate(SimputCtlg* const cat,
			  SimputSrc* const src,
			  const double time, 
			  const double mjdref,
			  int* const status);


/** Constructor for the SimputLC data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputLC* newSimputLC(int* const status);

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
SimputPSD* newSimputPSD(int* const status);

/** Destructor for the SimputPSD. Calls destructor routines for all
    contained elements, releases the allocated memory, and finally
    sets the pointer to NULL. */
void freeSimputPSD(SimputPSD** const lc);

/** Load the SimputPSD from the specified file. */
SimputPSD* loadSimputPSD(const char* const filename, int* const status);

/** Save the PSD in the specified extension of the given FITS
    file. If the file does not exist yet, a new file is created. If
    the file exists, an appropriate HDU is created. */
void saveSimputPSD(SimputPSD* const psd, 
		   const char* const filename,
		   char* const extname, 
		   int extver,
		   int* const status);


/** Constructor for the SimputImg data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputImg* newSimputImg(int* const status);

/** Destructor for the SimputImg data structure. Calls destructor
    routines for all contained elements, releases the allocated
    memory, and finally sets the pointer to NULL. */
void freeSimputImg(SimputImg** const img);

/** Load a SIMPUT source image from the specified file and store it in
    a SimputImg data structure. */
SimputImg* loadSimputImg(const char* const filename, int* const status);

/** Save the source image in the specified extension of the given FITS
    file. If the file does not exist yet, a new file is created. If
    the file exists, an appropriate HDU is created. */
void saveSimputImg(SimputImg* const img, 
		   const char* const filename,
		   char* const extname, 
		   int extver,
		   int* const status);


/** Constructor for the SimputPhList data structure. Allocates memory,
    initializes elements with their default values and pointers with
    NULL. */
SimputPhList* newSimputPhList(int* const status);

/** Destructor for the SimputPhList data structure. Calls destructor
    routines for all contained elements, releases the allocated
    memory, and finally sets the pointer to NULL. */
void freeSimputPhList(SimputPhList** const phl, int* const status);

/** Open a SimputPhList from a FITS file extension. */
SimputPhList* openSimputPhList(const char* const filename,
			       const int mode,
			       int* const status);


/** Return the maximum angular extension (radius) of a particular
    source around its reference point in [rad]. */
float getSimputSrcExt(SimputCtlg* const cat,
		      const SimputSrc* const src,
		      const double prevtime,
		      const double mjdref,
		      int* const status);


/* If no photon can be produced, because there is no light curve
   information available for the specified point of time, the return
   value of the function will be 1. If a photon is successfully
   produced, the return value will be 0. */
int getSimputPhotonTime(SimputCtlg* const cat,
			SimputSrc* const src,
			double prevtime,
			const double mjdref,
			double* const nexttime,
			int* const status);

/** Determine the energy and the direction of origin of a new
    photon. These two calculations have to be combined in one step,
    because, e.g., a photon list contains both information. */
void getSimputPhotonEnergyCoord(SimputCtlg* const cat,
				SimputSrc* const src,
				double currtime,
				const double mjdref,
				float* const energy,
				double* const ra,
				double* const dec,
				int* const status);

/** Produce a photon for a particular source in a SIMPUT catalog. The
    error status variable refers to errors related to the access to
    FITS files. If no photon can be produced, because there is no
    light curve information available for the specified point of time,
    the return value of the function will be 1. If a photon is
    successfully produced, the return value will be 0. */
int getSimputPhoton(SimputCtlg* const cat,
		    SimputSrc* const src,
		    const double prevtime,
		    const double mjdref,
		    /** [s]. */
		    double* const time,
		    /** [keV]. */
		    float* const energy,
		    /** [rad]. */
		    double* const ra,
		    /** [rad]. */
		    double* const dec,
		    int* const status);


#endif /* SIMPUT_H */
