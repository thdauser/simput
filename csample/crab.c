#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simput.h"

#define MAXSTR (1024)

#define CHECK_STATUS_BREAK(a) if (EXIT_SUCCESS!=a) break

#define CHECK_NULL_BREAK(a,status,msg) \
  if (NULL==a) { \
    printf(msg); \
    status=EXIT_FAILURE; \
    break;\
  }


int main( )
{
  const char filename[] = "crab.fits";

  SimputCatalog* cat=NULL;
  SimputMIdpSpec* spec=NULL;
  SimputLC* lc=NULL;
  int status=EXIT_SUCCESS;

  do { // Error handling loop.

    // Create a source catalog with a single source.
    cat = openSimputCatalog(filename, READWRITE, &status);
    CHECK_STATUS_BREAK(status);
    SimputSource* src = 
      getSimputSourceV(1, "PSRB0531+21", 
		       83.633125*M_PI/180., 22.014472*M_PI/180., 
		       0., 1., 
		       1., 10., 2.79674e-8, 
		       "[SPECTRUM,1]", "", "[LIGHTCUR,1]", 
		       &status);
    CHECK_STATUS_BREAK(status);
    
    appendSimputSource(cat, src, &status);
    CHECK_STATUS_BREAK(status);
    freeSimputCatalog(&cat, &status);
    CHECK_STATUS_BREAK(status);
    // END of catalog creation.


    // Create a spectrum and append it to the file with the source catalog.
    // TODO
    // END of spectrum creation.


    // Create a light curve and append it to the file with the source catalog.
    // TODO
    // END of light curve creation.

  } while(0); // END of error handling loop.
  
  freeSimputLC(&lc);
  freeSimputMIdpSpec(&spec);
  freeSimputCatalog(&cat, &status);
  
  fits_report_error(stderr, status);
  return(status);
}
