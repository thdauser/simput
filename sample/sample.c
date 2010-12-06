#include <stdio.h>
#include <stdlib.h>

#include "simput.h"

#define MAXMSG (1024)

#define CHECK_STATUS(a) if (EXIT_SUCCESS!=a) break

int main()
{
  const char filename[] = "simput.fits";
  int status=EXIT_SUCCESS;

  do { // Error handling loop.

    // Create a source catalog with a single source linked to the
    // spectrum and the light curve.
    simput_add_src(filename, 1, "myPOINTSOURCE",
    		   0., 0., 1.e-10, 1., 10.,
    		   "", "", "", &status);
    CHECK_STATUS(status);


    // Create a spectrum and add it to the file.
    float spec1_e_min[] = { 0., 2., 5. };
    float spec1_e_max[] = { 2., 5., 10. };
    float spec1_flux[]  = { 10., 5., 2. };
    int spec1_extver;
    simput_store_spectrum(filename, "SPECTRUM", 3,
			  spec1_e_min, spec1_e_max, spec1_flux,
			  0., &spec1_extver, &status);
    CHECK_STATUS(status);
    // Create a 2nd spectrum and add it to the file.
    float spec2_e_min[] = { 0., 3., 6. };
    float spec2_e_max[] = { 3., 6., 10. };
    float spec2_flux[]  = { 1., 4., 3. };
    int spec2_extver;
    simput_store_spectrum(filename, "SPECTRUM", 3,
			  spec2_e_min, spec2_e_max, spec2_flux,
			  0., &spec2_extver, &status);
    CHECK_STATUS(status);
    // Assign the 2 spectra to the source.
    char spec1_filename[MAXMSG];
    sprintf(spec1_filename, "%s[%s, %d]", filename, "SPECTRUM", spec1_extver);
    simput_add_spectrum(filename, 1, spec1_filename, &status);
    CHECK_STATUS(status);
    char spec2_filename[MAXMSG];
    sprintf(spec2_filename, "%s[%s, %d]", filename, "SPECTRUM", spec2_extver);
    simput_add_spectrum(filename, 1, spec2_filename, &status);
    CHECK_STATUS(status);


    // Create a lightcurve and add it to the file.
    double time[] = { 0., 1., 2., 3., 4. };
    float flux_lc[] = { 0.8, 0.7, 1.0, 1.2, 1.0 };
    int hdunum_lc;
    simput_store_lightcur(filename, 5,
			  time, NULL, flux_lc, NULL, NULL,
			  1., 10.,
			  &hdunum_lc, &status);
    CHECK_STATUS(status);

  } while(0); // END of error handling loop.

  fits_report_error(stderr, status);
  return(status);
}
