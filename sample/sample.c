#include <stdio.h>
#include <stdlib.h>

#include "simput.h"

#define MAXMSG (1024)

#define CHECK_STATUS(a) if (EXIT_SUCCESS!=a) break

int main()
{
  int status=EXIT_SUCCESS;

  do { // Error handling loop.

    // Create a source catalog with 2 sources.
    const char filename[] = "simput.fits";
    remove(filename);
    simput_add_src(filename, 1, "myPOINTSOURCE",
    		   0., 0., 1.e-10, 1., 10.,
    		   "", "", "", &status);
    simput_add_src(filename, 2, "2ndPOINTSOURCE",
    		   1., 1.5, 2.5e-9, 1.5, 9.5,
    		   "", "", "", &status);
    CHECK_STATUS(status);


    // Create a spectrum and add it to the source catalog file.
    float spec1_e_min[] = { 0., 2., 5. };
    float spec1_e_max[] = { 2., 5., 10. };
    float spec1_flux[]  = { 10., 5., 2. };
    int spec1_extver;
    simput_store_spectrum(filename, "SPECTRUM", 3,
			  spec1_e_min, spec1_e_max, spec1_flux,
			  0., &spec1_extver, &status);
    CHECK_STATUS(status);
    // Create a 2nd spectrum and add it to the source catalog file.
    float spec2_e_min[] = { 0., 3., 6. };
    float spec2_e_max[] = { 3., 6., 10. };
    float spec2_flux[]  = { 1., 4., 3. };
    int spec2_extver;
    simput_store_spectrum(filename, "SPECTRUM", 3,
			  spec2_e_min, spec2_e_max, spec2_flux,
			  0., &spec2_extver, &status);
    CHECK_STATUS(status);
    // Create a 3rd spectrum and add it to a separate file.
    float spec3_e_min[] = { 0., 1.5, 6.5 };
    float spec3_e_max[] = { 1.5, 6.5, 12. };
    float spec3_flux[]  = { 0.8, 1.2, 0.6 };
    const char spec3_filename[] = "spectrum.fits";
    remove(spec3_filename);
    int spec3_extver=0;
    simput_store_spectrum(spec3_filename, "SPECTRUM", 3,
			  spec3_e_min, spec3_e_max, spec3_flux,
			  0., &spec3_extver, &status);
    CHECK_STATUS(status);
    // Assign the 3 spectra to the 1st source.
    char spec_filename[MAXMSG];
    sprintf(spec_filename, "%s[%s, %d]", filename, "SPECTRUM", spec1_extver);
    simput_add_spectrum(filename, 1, spec_filename, &status);
    CHECK_STATUS(status);
    sprintf(spec_filename, "%s[%s, %d]", filename, "SPECTRUM", spec2_extver);
    simput_add_spectrum(filename, 1, spec_filename, &status);
    CHECK_STATUS(status);
    sprintf(spec_filename, "%s[%s, %d]", spec3_filename, "SPECTRUM", spec3_extver);
    simput_add_spectrum(filename, 1, spec_filename, &status);
    CHECK_STATUS(status);
    // Assign the 3rd spectrum also to the 2nd source.
    sprintf(spec_filename, "%s[%s, %d]", spec3_filename, "SPECTRUM", spec3_extver);
    simput_add_spectrum(filename, 2, spec_filename, &status);
    CHECK_STATUS(status);


    // Create a lightcurve and add it to the file.
    double time[] = { 0., 1., 2., 3., 4. };
    float flux_lc[] = { 0.8, 0.7, 1.0, 1.2, 1.0 };
    int lc_extver;
    simput_store_lightcur(filename, "LIGHTCUR", 5,
			  time, NULL, flux_lc, NULL, NULL,
			  1., 10., &lc_extver, &status);
    CHECK_STATUS(status);
    // Assign the light curve to the 1st source.
    char lc_filename[MAXMSG];
    sprintf(lc_filename, "%s[%s, %d]", filename, "LIGHTCUR", lc_extver);
    simput_add_lightcur(filename, 1, lc_filename, &status);
    CHECK_STATUS(status);

    
    // Assign source images to the 2nd source.
    simput_add_image(filename, 2, "image1.fits", &status);
    CHECK_STATUS(status);
    simput_add_image(filename, 2, "image2.fits", &status);
    CHECK_STATUS(status);

    // TODO Create a second source.

  } while(0); // END of error handling loop.

  fits_report_error(stderr, status);
  return(status);
}
