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


void printSimputSourceEntry(SimputSourceEntry* sse) {
  printf("Simput Source Entry:\n");
  printf(" SRC_ID:\t%ld\n", sse->src_id);
  printf(" SRC_NAME:\t'%s'\n", sse->src_name);
  printf(" RA:\t\t%lf\n", sse->ra);
  printf(" DEC:\t\t%lf\n", sse->dec);
  printf(" IMGROTA:\t%f\n", sse->imgrota);
  printf(" IMGSCAL:\t%f\n", sse->imgscal);
  printf(" E_MIN:\t\t%f\n", sse->e_min);
  printf(" E_MAX:\t\t%f\n", sse->e_max);
  printf(" FLUX:\t\t%f\n", sse->flux);
  printf(" SPECTRUM:\t'%s'\n", sse->spectrum);
  printf(" IMAGE:\t\t'%s'\n", sse->image);
  printf(" LIGHTCUR:\t'%s'\n", sse->lightcur);
}


int main(int argc, char **argv)
{
  const char filename[] = "simput.fits";

  SimputSourceCatalog* catalog =NULL;
  SimputMissionIndepSpec* spec =NULL;
  SimputMissionIndepSpec* spec2=NULL;
  SimputImg* img=NULL;
  int status=EXIT_SUCCESS;

  do { // Error handling loop.

    if ((argc>1) && (0==strcmp(argv[1],"r"))) {
      // Read the source catalog from the file.
      catalog = loadSimputSourceCatalog(filename, &status);
      CHECK_STATUS_BREAK(status);

      printf("catalog with %ld entries successfully loaded\n", catalog->nentries);
      printSimputSourceEntry(catalog->entries[0]);
      printSimputSourceEntry(catalog->entries[1]);
    } else {

      remove(filename);

      // Create a source catalog with 2 sources.
      catalog = getSimputSourceCatalog(&status);
      CHECK_STATUS_BREAK(status);
      catalog->entries = (SimputSourceEntry**)malloc(2*sizeof(SimputSourceEntry*));
      CHECK_NULL_BREAK(catalog->entries, status, "memory allocation failed!\n");
      catalog->entries[0] = 
	getSimputSourceEntryV(1, "", -1.*M_PI/180., 2.5*M_PI/180., 0., 1., 
			      1., 10., 5.e-12, "", "", "", &status);
      CHECK_STATUS_BREAK(status);
      catalog->entries[1] = 
	getSimputSourceEntryV(2, "", -1.8*M_PI/180., 2.3*M_PI/180., 30.*M_PI/180., 1.2, 
			      0.5, 15., 8.e-13, "", "", "", &status);
      CHECK_STATUS_BREAK(status);
      catalog->nentries = 2;

      saveSimputSourceCatalog(catalog, filename, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a source catalog.

      // Create a spectrum and append it to the file with the source catalog.
      spec = getSimputMissionIndepSpec(&status);
      CHECK_STATUS_BREAK(status);

      spec->nentries=2;
      spec->energy  = (float*)malloc(2*sizeof(float));
      CHECK_NULL_BREAK(spec->energy, status, "memory allocation failed!\n");
      spec->flux    = (float*)malloc(2*sizeof(float));
      CHECK_NULL_BREAK(spec->flux, status, "memory allocation failed!\n");
      spec->name    = (char*)malloc(1024*sizeof(char));
      CHECK_NULL_BREAK(spec->name, status, "memory allocation failed!\n");

      spec->energy[0] = 1.;
      spec->flux[0]   = 0.8;
      spec->energy[1] = 2.;
      spec->flux[1]   = 0.6;
      strcpy(spec->name, "spec1");

      saveSimputMissionIndepSpec(spec, filename, "SPEC", 1, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a spectrum.

      // Create a 2nd spectrum and append it to the file with the source catalog.
      spec2 = getSimputMissionIndepSpec(&status);
      CHECK_STATUS_BREAK(status);

      spec2->nentries=3;
      spec2->energy  = (float*)malloc(3*sizeof(float));
      CHECK_NULL_BREAK(spec2->energy, status, "memory allocation failed!\n");
      spec2->flux    = (float*)malloc(3*sizeof(float));
      CHECK_NULL_BREAK(spec2->flux, status, "memory allocation failed!\n");
      spec2->name    = (char*)malloc(1024*sizeof(char));
      CHECK_NULL_BREAK(spec2->name, status, "memory allocation failed!\n");

      spec2->energy[0] = 1.5;
      spec2->flux[0]   = 0.75;
      spec2->energy[1] = 2.5;
      spec2->flux[1]   = 0.65;
      spec2->energy[2] = 3.5;
      spec2->flux[2]   = 0.55;
      strcpy(spec2->name, "spec2");

      saveSimputMissionIndepSpec(spec2, filename, "SPEC", 1, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a spectrum.

      // Create a source image and append it to the file with the source catalog.
      img = getSimputImg(&status);
      CHECK_STATUS_BREAK(status);
      img->dist = (double**)malloc(3*sizeof(double*));
      CHECK_NULL_BREAK(img->dist, status, "memory allocation failed!\n");
      long ii;
      for (ii=0; ii<3; ii++) {
	img->dist[ii] = (double*)malloc(sizeof(double));
	CHECK_NULL_BREAK(img->dist[ii], status, "memory allocation failed!\n");
      }
      CHECK_STATUS_BREAK(status);
      img->naxis1   = 3;
      img->naxis2   = 1;
      img->dist[0][0] = 0.6;
      img->dist[1][0] = 1.6;
      img->dist[2][0] = 2.0;
      img->fluxscal = 2.0;

      //      saveSimputImg(img, filename, "IMG1", 0, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a source image.
    }
    // END of creating a new file.

  } while(0); // END of error handling loop.

  freeSimputImg(&img);
  freeSimputMissionIndepSpec(&spec);
  freeSimputMissionIndepSpec(&spec2);
  freeSimputSourceCatalog(&catalog);
  
  fits_report_error(stderr, status);
  return(status);
}
