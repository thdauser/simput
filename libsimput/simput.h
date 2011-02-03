#ifndef SIMPUT_H
#define SIMPUT_H (1)

#include "fitsio.h"

/*
  simput.h

  Version  Date       Author                
  --------------------------------------------------------------------------- 
  0.0.1  2011/02/02   Christian Schmid      initial 
  --------------------------------------------------------------------------- 

   This is the header file of the simput library. The library provides
   basic routines to create a SIMPUT source catalog file according to
   the standard defined in ???.
*/

/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

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
void simput_store_spectrum(/** Name of the output file. If the file
			       already exists, a new table is
			       appended. */
			   const char* const filename,
			   /** Extension name of the new HDU. EXTNAM
			       and EXTVER can be used to assign the
			       spectrum to an X-ray source in the
			       catalog via the extended filename
			       syntax. */
			   char* const extname,
			   /** Number of bins in the spectrum. */
			   const long nbins,
			   /** Vector containing the lower boundaries
			       of the energy bins in [keV]. */
			   float* const e_min,
			   /** Vector containing the upper boundaries
			       of the energy bins in [keV]. */
			   float* const e_max,
			   /** Vector containing the source flux
			       density in [photons/s/cm^2/keV] in each
			       bin. */
			   float* const flux,
			   /** If the spectrum is part of a set
			       describing spectral variations of an
			       oscillating source, this parameter is
			       used to refer to the appropriate phase
			       the spectrum is valid for. If this
			       parameter is not used, it should be set
			       to 0. */
			   float phase,
			   /** Extension version of the new HDU. This
			       parameter cannot be set by the user,
			       but is a return value of this
			       function. EXTNAM and EXTVER can be used
			       to assign the spectrum to an X-ray
			       source in the catalog via the extended
			       filename syntax. */
			   int* extver,
			   /** Error status. */
			   int* const status);


/** Store a light curve in the specified file. */
void simput_store_lightcur(/** Name of the output file. If the file
			       already exists, a new table is
			       appended. */
			   const char* const filename,
			   /** Extension name of the new HDU. EXTNAM
			       and EXTVER can be used to assign the
			       light curve to an X-ray source in the
			       catalog via the extended filename
			       syntax. */
			   char* const extname,
			   /** Number of bins in the light curve. */
			   const long nbins,
			   /** Vector containing the time of the
			       individual bins of the light curve in
			       [s]. */
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
			   /** Vector containing the polarization
			       fraction in each bin. */
			   float* const pol_frac,
			   /** Vector containing the polarization
			       direction in [deg] in each bin. */
			   float* const pol_dir,
			   /** Lower boundary of the reference energy
			       band. */
			   float e_min,
			   /** Upper boundary of the reference energy
			       band. */
			   float e_max,
			   /** Extension version of the new HDU. This
			       parameter cannot be set by the user,
			       but is a return value of this
			       function. EXTNAM and EXTVER can be used
			       to assign the light curve to an X-ray
			       source in the catalog via the extended
			       filename syntax. */
			   int* extver,
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
