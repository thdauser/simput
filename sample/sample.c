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


int main(int argc, char **argv)
{
  const char filename[] = "simput.fits";

  SimputSourceCatalog* catalog=NULL;
  int status=EXIT_SUCCESS;

  do { // Error handling loop.

    if ((argc>1) && (0==strcmp(argv[1],"r"))) {
      // Read the source catalog from the file.
      catalog = loadSimputSourceCatalog(filename, &status);
      CHECK_STATUS_BREAK(status);

      printf("catalog with %d entries successfully loaded\n", catalog->nentries);

    } else {
      // Create a source catalog with 2 sources.
      catalog = getSimputSourceCatalog(&status);
      CHECK_STATUS_BREAK(status);
      catalog->entries = (SimputSourceEntry**)malloc(2*sizeof(SimputSourceEntry*));
      CHECK_NULL_BREAK(catalog->entries, status, "Error: memory allocation failed!\n");
      catalog->entries[0] = 
	getSimputSourceEntryV(1, "", -1.*M_PI/180., 2.5*M_PI/180., 0., 1., 
			      1., 10., 5.e-12, "", "", "", &status);
      CHECK_STATUS_BREAK(status);
      catalog->entries[1] = 
	getSimputSourceEntryV(2, "", -1.8*M_PI/180., 2.3*M_PI/180., 30.*M_PI/180., 1.2, 
			      0.5, 15., 8.e-13, "", "", "", &status);
      CHECK_STATUS_BREAK(status);
      catalog->nentries = 2;

      //    remove(filename);
      saveSimputSourceCatalog(catalog, filename, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of create a source catalog.

  } while(0); // END of error handling loop.

  freeSimputSourceCatalog(&catalog);

  fits_report_error(stderr, status);
  return(status);
}
