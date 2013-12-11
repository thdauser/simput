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


void printSimputSource(SimputSrc* sse) {
  printf("Simput Source :\n");
  printf(" SRC_ID:\t%ld\n", sse->src_id);
  printf(" SRC_NAME:\t'%s'\n", sse->src_name);
  printf(" RA:\t\t%lf\n", sse->ra);
  printf(" DEC:\t\t%lf\n", sse->dec);
  printf(" IMGROTA:\t%f\n", sse->imgrota);
  printf(" IMGSCAL:\t%f\n", sse->imgscal);
  printf(" E_MIN:\t\t%f\n", sse->e_min);
  printf(" E_MAX:\t\t%f\n", sse->e_max);
  printf(" FLUX:\t\t%f\n", sse->eflux);
  printf(" SPECTRUM:\t'%s'\n", sse->spectrum);
  printf(" IMAGE:\t\t'%s'\n", sse->image);
  printf(" TIMING:\t'%s'\n", sse->timing);
}


int main(int argc, char **argv)
{
  const char filename[]="simput.fits";

  SimputCtlg* cat=NULL;
  SimputMIdpSpec* spec =NULL;
  SimputMIdpSpec* spec2=NULL;
  SimputImg* img=NULL;
  int status=EXIT_SUCCESS;

  do { // Error handling loop.

    if ((argc>1) && (0==strcmp(argv[1],"r"))) {

      // Read the source catalog from the file.
      cat=openSimputCtlg(filename, READONLY, 0, 0, 0, 0, &status);
      CHECK_STATUS_BREAK(status);

      printf("catalog contains %ld entries\n", cat->nentries);

      SimputSrc* src=getSimputSrc(cat, 1, &status);
      CHECK_STATUS_BREAK(status);
      printSimputSource(src);

      src=getSimputSrc(cat, 2, &status);
      CHECK_STATUS_BREAK(status);
      printSimputSource(src);

    } else {

      remove(filename);

      // Create a source catalog with 2 sources.
      cat=openSimputCtlg(filename, READWRITE, 0, 0, 0, 0, &status);
      CHECK_STATUS_BREAK(status);

      SimputSrc* src=
	newSimputSrcV(1, "", -1.*M_PI/180., 2.5*M_PI/180., 0., 1., 
		      1., 10., 5.e-12, "[SPEC,1][#row==1]", "[IMG1,1]", "", &status);
      CHECK_STATUS_BREAK(status);
      appendSimputSrc(cat, src, &status);
      CHECK_STATUS_BREAK(status);

      src=newSimputSrcV(2, "", -1.8*M_PI/180., 2.3*M_PI/180., 30.*M_PI/180., 1.2,
			0.5, 15., 8.e-13, "[SPEC,1][#row==2]", "", "", &status);
      CHECK_STATUS_BREAK(status);
      appendSimputSrc(cat, src, &status);
      CHECK_STATUS_BREAK(status);

      freeSimputCtlg(&cat, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a source catalog.

      // Create a spectrum and append it to the file with the source catalog.
      spec=newSimputMIdpSpec(&status);
      CHECK_STATUS_BREAK(status);

      spec->nentries=2;
      spec->energy=(float*)malloc(2*sizeof(float));
      CHECK_NULL_BREAK(spec->energy, status, "memory allocation failed!\n");
      spec->fluxdensity=(float*)malloc(2*sizeof(float));
      CHECK_NULL_BREAK(spec->fluxdensity, status, "memory allocation failed!\n");
      spec->name=(char*)malloc(1024*sizeof(char));
      CHECK_NULL_BREAK(spec->name, status, "memory allocation failed!\n");

      spec->energy[0]=1.0;
      spec->fluxdensity[0]=0.8;
      spec->energy[1]=2.0;
      spec->fluxdensity[1]=0.6;
      strcpy(spec->name, "spec1");

      saveSimputMIdpSpec(spec, filename, "SPEC", 1, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a spectrum.

      // Create a 2nd spectrum and append it to the file with the source catalog.
      spec2=newSimputMIdpSpec(&status);
      CHECK_STATUS_BREAK(status);

      spec2->nentries=3;
      spec2->energy=(float*)malloc(3*sizeof(float));
      CHECK_NULL_BREAK(spec2->energy, status, "memory allocation failed!\n");
      spec2->fluxdensity=(float*)malloc(3*sizeof(float));
      CHECK_NULL_BREAK(spec2->fluxdensity, status, "memory allocation failed!\n");
      spec2->name=(char*)malloc(1024*sizeof(char));
      CHECK_NULL_BREAK(spec2->name, status, "memory allocation failed!\n");

      spec2->energy[0]=1.5;
      spec2->fluxdensity[0]=0.75;
      spec2->energy[1]=2.5;
      spec2->fluxdensity[1]=0.65;
      spec2->energy[2]=3.5;
      spec2->fluxdensity[2]=0.55;
      strcpy(spec2->name, "spec2");

      saveSimputMIdpSpec(spec2, filename, "SPEC", 1, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a spectrum.

      // Create a source image and append it to the file with the source catalog.
      img=newSimputImg(&status);
      CHECK_STATUS_BREAK(status);
      img->dist=(double**)malloc(3*sizeof(double*));
      CHECK_NULL_BREAK(img->dist, status, "memory allocation failed!\n");
      long ii;
      for (ii=0; ii<3; ii++) {
	img->dist[ii]=(double*)malloc(sizeof(double));
	CHECK_NULL_BREAK(img->dist[ii], status, "memory allocation failed!\n");
      }
      CHECK_STATUS_BREAK(status);
      img->naxis1=3;
      img->naxis2=1;
      img->dist[0][0]=0.6;
      img->dist[1][0]=1.6;
      img->dist[2][0]=2.0;

      img->wcs=(struct wcsprm*)malloc(sizeof(struct wcsprm));
      if (NULL==img->wcs) {
	printf("error: memory allocation for WCS data structure failed!\n");
	status=EXIT_FAILURE;
	break;
      }
      img->wcs->flag=-1;
      if (0!=wcsini(1, 2, img->wcs)) {
	printf("error: initalization of WCS data structure failed!\n");
	status=EXIT_FAILURE;
	break;
      }
      img->wcs->naxis=2;
      strcpy(img->wcs->ctype[0], "RA---TAN");
      strcpy(img->wcs->ctype[1], "DEC--TAN");
      strcpy(img->wcs->cunit[0], "deg");
      strcpy(img->wcs->cunit[1], "deg");
      img->wcs->crval[0]=0.0;
      img->wcs->crval[1]=0.0;
      img->wcs->crpix[0]=2.0;
      img->wcs->crpix[1]=1.0;
      img->wcs->cdelt[0]=0.01;
      img->wcs->cdelt[1]=0.01;

      saveSimputImg(img, filename, "IMG1", 0, &status);
      CHECK_STATUS_BREAK(status);
      // END of create a source image.
    }
    // END of creating a new file.

  } while(0); // END of error handling loop.

  freeSimputImg(&img);
  freeSimputMIdpSpec(&spec);
  freeSimputMIdpSpec(&spec2);
  freeSimputCtlg(&cat, &status);
  
  fits_report_error(stderr, status);
  return(status);
}
