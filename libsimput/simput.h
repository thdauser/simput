#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"

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


// Obsolete:
/** Common string length. */
#define SIMPUT_MAXMSG (1024)


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  /** Unique source ID. */
  unsigned int src_id;

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
  unsigned int nentries;

  /** Array of the individual entries in the catalog. */
  SimputSourceEntry** entries;

} SimputSourceCatalog;



// Obsolete:
/** Time definition keywords. Defined by Angelini et al. (1994). */
struct simput_timing {
  double mjdref;
  double tstart;
  double tstop;
  double timezero;
  char timesys[SIMPUT_MAXMSG];
  char timeunit[SIMPUT_MAXMSG];
  char clockcor[SIMPUT_MAXMSG];
};


/** Polarization of emitted radiation. The Stokes parameters are
    defined by Stokes (1852). */
struct simput_polarization {
  /** 2nd component of normalized Stokes vector. */
  double stokes1;
  /** 3rd component of normalized Stokes vector. */
  double stokes2;
  /** 4th component of normalized Stokes vector. */
  double stokes3;
  /** Reference coordinate system. */
  char radecsys[SIMPUT_MAXMSG];
  char equinox[SIMPUT_MAXMSG];
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for the SimputSourceEntry data structure. Allocates
    memory, initializes elements with their default values and
    pointers with NULL. */
SimputSourceEntry* getSimputSourceEntry(int* const status);

/** Constructor for the SimputSourceEntry data structure. Allocates
    memory and initializes elements with the given values. */
SimputSourceEntry* getSimputSourceEntryV(const unsigned int src_id, 
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


// Obsolete:
/** Add a new X-ray source to the source catalog. Insert a new line
    with with the basic description of the X-ray source in a source
    catalog file. If the file or the source catalog extension does not
    exists yet, a new file with the appropriate extension is
    created. The ID of the source (src_id) must be unique. If a source
    with the same ID already exists in the catalog, an error is
    returned. */
void simput_add_src(/** Name of the FITS file containing the source
			catalog. */
		    const char* const filename,
		    /** Unique identifier of the X-ray source. */
		    long src_id,
		    /** Name of the X-ray source. The string might be
			empty. */
		    char* src_name,
		    /** Right ascension. */
		    float ra,
		    /** Declination. */
		    float dec,
		    /** Source flux density in the specified reference
			energy band [erg/s/cm^2]. */
		    float flux,
		    /** Lower boundary of the reference energy band in
			[keV]. */
		    float e_min,
		    /** Upper boundary of the reference energy band in
			[keV]. */
		    float e_max,
		    /** Extended filename of the FITS file extension
			containing the energy spectrum of the X-ray
			source. A source can be described by several
			spectra. In that case this field must contain
			a reference to a grouping table containing all
			the spectra. */
		    char* spectrum,
		    /** Extended filename of the FITS file extension
			containing an image of the X-ray source. A
			source can be described by several images. In
			that case this field must contain a reference
			to a grouping table containing all the
			spectra. If this field is empty, the X-ray
			source is assumed to be point-like. */
		    char* image,
		    /** Extended filename of the FITS file extension
			containing the light curve of the X-ray
			source. If this field is empty, the X-ray
			source is assumed to have a constant
			brightness. */
		    char* lightcur,
		    /** Error status. */
		    int* const status);


/** Store a flux density spectrum in the specified FITS file. */
void simput_write_spectrum(/** Name of the output file. If the file
			       already exists, a new table is
			       appended. */
			   const char* const filename,
			   /** Extension name of the new HDU. EXTNAM
			       and EXTVER can be used to assign the
			       spectrum to an X-ray source in the
			       catalog via the extended filename
			       syntax. */
			   char* const extname,
			   /** Extension version of the new HDU. This
			       parameter cannot be set by the user,
			       but is a return value of this
			       function. EXTNAM and EXTVER can be used
			       to assign the spectrum to an X-ray
			       source in the catalog via the extended
			       filename syntax. */
			   int* extver,
			   /** Timing parameters. Required if the
			       spectrum is valid for a particular
			       point of time. If not, a NULL pointer
			       should be given here. */
			   struct simput_timing* const timing,
			   /** If the spectrum is part of a set
			       describing spectral variations of an
			       oscillating source, this parameter is
			       used to refer to the appropriate phase
			       the spectrum is valid for. If this
			       parameter is not used, the pointer
			       should be set to NULL. */
			   float* phase,
			   /** The polarization of the emitted
			       radiation can be specified by a set of
			       Stokes parameters. */
			   struct simput_polarization* const polarization,
			   /** Number of bins in the spectrum. */
			   const long nbins,
			   /** Vector containing the lower boundaries
			       of the energy bins in [keV]. */
			   float* const e_min,
			   /** Vector containing the upper boundaries
			       of the energy bins in [keV]. */
			   float* const e_max,
			   /** Vector containing the spectral
			       distribution of the photon flux density
			       in [photons/s/cm^2/keV] for each energy
			       bin. */
			   float* const flux,
			   /** Error status. */
			   int* const status);


/** Store a light curve in the specified file. */
void simput_write_lightcur(/** Name of the output file. If the file
			       already exists, a new table is
			       appended. */
			   const char* const filename,
			   /** Extension name of the new HDU. EXTNAM
			       and EXTVER can be used to assign the
			       light curve to an X-ray source in the
			       catalog via the extended filename
			       syntax. */
			   char* const extname,
			   /** Extension version of the new HDU. This
			       parameter cannot be set by the user,
			       but is a return value of this
			       function. EXTNAM and EXTVER can be used
			       to assign the light curve to an X-ray
			       source in the catalog via the extended
			       filename syntax. */
			   int* extver,
			   /** Timing definition parameters. */
			   struct simput_timing* const timing,
			   /** Lower boundary of the reference energy
			       band in [keV]. */
			   float e_min,
			   /** Upper boundary of the reference energy
			       band in [keV]. */
			   float e_max,
			   /** Flux density in the reference energy
			       band in [erg/s/cm^2]. */
			   float ref_flux,
			   /** Phase of periodic oscillation at
			       TIMEZERO (periodic light curves
			       only). */
			   double* phase0,
			   /** Duration of one oscillation period
			       (periodic light curves only). */
			   double* period,
			   /** Number of bins in the light curve. */
			   const long nbins,
			   /** Vector containing the time of the
			       individual bins of the light curve in
			       [s]. This parameter only applies for
			       non-periodic light curves. Should be
			       NULL for periodic light curves. */
			   double* const time,
			   /** Vector containing the phase of the
			       individual bins of the light
			       curve. This parameter only applies for
			       periodic systems. For non-periodic
			       systems it should be NULL. */
			   float* const phase,
			   /** Vector containing the source flux
			       density in [erg/s/cm^2] in each bin of
			       the light curve, determined for the
			       reference energy band. */			   
			   float* const flux,
			   /** Error status. */
			   int* const status);


/** Insert a reference to a spectrum into a specific line in the
    source catalog table. If the source specification already contains
    a spectrum, combine the references to the multiple spectra
    in a grouping table. */
void simput_add_spectrum(/** File containing the source catalog. */
			 const char* const srcctlg_filename,
			 /** ID of the X-ray source. */
			 const long src_id,
			 /** Extended file name of the FITS file
			     extension containing the spectrum. */
			 const char* const spec_filename,
			 /** Error status. */
			 int* const status);


/** Insert a reference to a light curve into a specific line in the
    source catalog table. */
void simput_add_lightcur(/** File containing the source catalog. */
			 const char* const srcctlg_filename,
			 /** ID of the X-ray source. */
			 const long src_id,
			 /** Extended file name of the FITS file
			     extension containing the light curve. */
			 const char* const lc_filename,
			 /** Error status. */
			 int* const status);


/** Insert a reference to a source image into a specific line in the
    source catalog table. */
void simput_add_image(/** File containing the source catalog. */
		      const char* const srcctlg_filename,
		      /** ID of the X-ray source. */
		      const long src_id,
		      /** Extended file name of the FITS file
			  extension containing the image. */
		      const char* const img_filename,
		      /** Error status. */
		      int* const status);


#endif /* SIMPUT_H */
