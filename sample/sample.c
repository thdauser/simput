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
    float e_min[] = { 0., 2., 5. };
    float e_max[] = { 2., 5., 10. };
    float flux_spec[]  = { 10., 5., 2. };
    int extver_spec;
    simput_store_spectrum(filename, "SPECTRUM", 3,
			  e_min, e_max, flux_spec,
			  0., &extver_spec, &status);
    CHECK_STATUS(status);
    char spec_filename[MAXMSG];
    sprintf(spec_filename, "%s[%s, %d]", filename, "SPECTRUM", extver_spec);
    simput_add_spectrum(filename, 1, spec_filename, &status);
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
