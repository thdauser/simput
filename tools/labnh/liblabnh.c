#include "labnh.h"
#include "posstring.h"

#include <math.h>
#include <err.h>
#include <stdio.h>

#include <fitsio.h>
#include <wcshdr.h>
#include <wcs.h>
#include <wcsfix.h>

#include <atFunctions.h>

typedef struct {
  int nwcs;
  struct wcsprm *wcs;
  long naxes[3];
  long imgsiz;
  double vmin;
  double dv;
  double *nhval;
  double *nhcache;
} nhinfo;

static nhinfo *nhmap=NULL;

static int verbosity=0;
static int sampling=1;


//
// Helper Functions
//

double get21cm(long x,long y,long z);
double angdist(double l1, double b1, double l2, double b2);

// return one value in the full blown (velocity resolved)
// LAB map
// for speed reasons: noout of bounds checking and the like...
// -> this is an internal function!
double get21cm(long x,long y,long z) {
    long ind=x+nhmap->naxes[0]*y+nhmap->imgsiz*z;
    double xx=nhmap->nhval[ind];
    if (xx<=-11.02) {
      xx=0.;
    }
    return xx;
}

double angdist(double l1, double b1, double l2, double b2) {
  //
  // calculate the angular distance between points l1/b1, l2/b2
  // using the Vincenty formula (numerically stable for all cases)
  //

  // all input and output is in DEGREES

  double deg2rad=M_PI/180.;

  l1*=deg2rad; b1*=deg2rad;
  l2*=deg2rad; b2*=deg2rad;

  #ifndef _GNU_SOURCE
  double sinb1=sin(b1);
  double cosb1=cos(b1);

  double sinb2=sin(b2);
  double cosb2=cos(b2);

  double sindl=sin(l2-l1);
  double cosdl=cos(l2-l1);
  #else
  #if defined(__APPLE__) && defined(__MACH__)
  double sinb1,cosb1,sinb2,cosb2,sindl,cosdl;
  __sincos(b1,&sinb1,&cosb1);
  __sincos(b2,&sinb2,&cosb2);
  __sincos(l2-l1,&sindl,&cosdl);
  #else
  double sinb1,cosb1,sinb2,cosb2,sindl,cosdl;
  sincos(b1,&sinb1,&cosb1);
  sincos(b2,&sinb2,&cosb2);
  sincos(l2-l1,&sindl,&cosdl);
  #endif
  #endif

  double nom=sqrt(
		  pow(cosb1*sindl,2)+
		  pow(cosb2*sinb1-sinb2*cosb1*cosdl,2)
		  );
  double denom=sinb2*sinb1+cosb2*cosb1*cosdl;

  double dist=atan2(nom,denom);

  return dist/deg2rad; 
}


//
// Interface routines (exported)
//

void nh_verbosity(int verb) {
  verbosity=verb;
}

void nh_sampling(int sampl) {
  sampling=fabs(sampl); // must be positive
}


// initializer function. 
int nhinit() {
  if (verbosity>2) {
    warnx("Initializing LABNH function");
  }

  char *map=getenv("LABNH");

  if (map==NULL) {
    warnx("Environment variable LABNH not set");
    return 100;
  }
  int ierr=0;

  fitsfile *fptr;
  char *header;
  int nkeys,nreject,nwcs;
  int ifix,iwcs;
  int stat[NWCSFIX];
  struct wcsprm *wcs;

  nhmap=(nhinfo *) malloc(sizeof(nhinfo));

  // open map file
  if (ffopen(&fptr,map,0,&ierr) ) {
    fits_report_error(stderr,ierr);
    return 1;
  };

  // get header (without comments)
  if ( fits_hdr2str(fptr,1,NULL,0,&header,&nkeys,&ierr) ) {
    fits_report_error(stderr,ierr);
    ffclos(fptr,&ierr);
    return 1;
  };

  // extract get WCS keywords
  if ( (ierr=wcsbth(header,nkeys,WCSHDR_all,2,0,NULL,&nreject,&nwcs,&wcs)) ) {
    warnx("wcsbth ERROR %d: %s.\n", ierr, wcs_errmsg[ierr]);
    return 2;
  }

  // repair nonstandard keywords
  for (iwcs=0; iwcs<nwcs; iwcs++ ) {
    if ((ierr=wcsfix(7,0,wcs+iwcs,stat) )) {
      fprintf(stderr,"wcsfix ERROR, status returns: (");
      for (ifix = 0; ifix < NWCSFIX; ifix++) {
        fprintf(stderr,ifix ? ", %d" : "%d", stat[ifix]);
      }
      warnx(")\n");
    }
    // initialize conversions
    if ( (ierr=wcsset(wcs)) ) {
      warnx("Error setting up coordinate conversions: Error %i\n",ierr);
      return 3;
    };
  }

  nhmap->wcs=wcs;
  nhmap->nwcs=nwcs;

  //
  // WCS has worked, now read the full image extension
  //
  int bitpix;
  int naxis;
  if (fits_get_img_param(fptr,3,&bitpix,&naxis,nhmap->naxes,&ierr)) {
    fits_report_error(stderr,ierr);
    return 4;
  };

  if (naxis!=3) {
    warnx("FITS file has %i axes, but exactly three are required.\n",naxis);
    return 5;
  }
  
  // ... and suck it in
  LONGLONG numpxl=(nhmap->naxes[0])*(nhmap->naxes[1])*(nhmap->naxes[2]);
  long *fpixel=malloc(3*sizeof(long));
  fpixel[0]=1 ; fpixel[1]=1 ; fpixel[2]=1;
  double nan=FP_NAN;
  int anynul;
  nhmap->nhval=malloc(numpxl*sizeof(double));
  if (nhmap->nhval == 0 ) {
    warnx("Cannot allocate memory for nhmap\n");
    return 6;
  }

  if (fits_read_pix(fptr,TDOUBLE,fpixel,numpxl,&nan,
		    nhmap->nhval,&anynul,&ierr)) {
    fits_report_error(stderr,ierr);
    return 7;
  }
  nhmap->imgsiz=(nhmap->naxes[0])*(nhmap->naxes[1]);

  //
  // parameters for velocity grid
  //
  if (fits_read_key(fptr,TDOUBLE,"CRVAL3",&(nhmap->vmin),NULL,&ierr)){
    fits_report_error(stderr,ierr);
    return 8;
  }
  if (fits_read_key(fptr,TDOUBLE,"CDELT3",&(nhmap->dv),NULL,&ierr)){
    fits_report_error(stderr,ierr);
    return 9;
  }

  // allocate cache for integrated Nh values and
  // initialize it to 0 (all points on the sky have
  // NH>0, so this is safe)
  nhmap->nhcache=calloc(nhmap->imgsiz,sizeof(double));
  if (nhmap->nhcache == NULL ) {
    warnx("Cannot allocate memory for nhcache\n");
    return 10;
  }

  return 0;
}

void nhfinish() {
  if (verbosity>2) {
    warnx("Freeing NH map memory");
  }
  wcsvfree(&(nhmap->nwcs),&(nhmap->wcs));
  free(nhmap);
  nhmap=NULL;
}

double nh_equ(double ra, double dec) {
  double lii,bii;

  if (nhmap==NULL && nhinit(NULL)!=0) {
    return -1.;
  }

  // convert ra,dec to Galactic
  atJ2000toGal(ra,dec,&lii,&bii);

  if (verbosity>0) {
    fprintf(stderr," l/b: %8.3f -%8.3f\n",lii,bii);
  }

  return nh_gal(lii,bii);
}

double nh_gal(double lii, double bii) {
  double xpix,ypix;
  double world[3];
  double phi[1],theta[1],imgcrd[3],pixcrd[3];
  int status[3];

  // initialize if necessary
  if (nhmap==NULL && nhinit(NULL)!=0) {
    return -1.;
  }

  // Integrate over -400...+400 km/s
  double vstart=-400000.;
  double vstop=+400000.;

  int zmin=(int) ( (vstart-nhmap->vmin)/nhmap->dv );
  if (zmin<0) { zmin=0;}
  if (zmin>=nhmap->naxes[2]) { zmin=nhmap->naxes[2]-1;}

  int zmax=(int) ( (vstop-nhmap->vmin)/nhmap->dv );
  if (zmax<0) { zmax=0;}
  if (zmax>=nhmap->naxes[2]) { zmax=nhmap->naxes[2]-1;}

  if ( verbosity>0 ) {
    fprintf(stderr," v-range: %8.3f - %8.3f km/s\n",vstart/1000.,vstop/1000.);
  }
  if ( verbosity>1) {
    fprintf(stderr," z-range: %4i -%4i\n",zmin,zmax);
  }

  //
  // get nearest map pixel
  //
  world[nhmap->wcs->lng]=lii;
  world[nhmap->wcs->lat]=bii;
  world[2]=0;

  int ret=wcss2p(nhmap->wcs,1,3,world,phi,theta,imgcrd,pixcrd,status);

  if (ret != 0 ) {
    warnx("problem with wcss2p: %i",ret);
    return -1.;
  }

  // pixel position, incl. correction that FORTRAN array starts at 1!
  xpix=pixcrd[0]-1;
  ypix=pixcrd[1]-1;

  if (verbosity>0) {
    fflush(stdout);
  }

  if (verbosity>1) {
    fprintf(stderr," nearest map pixel: %3li/%3li\n",(long) xpix,(long) ypix);
  }

  // 
  // nearest neighbor weighing
  // (not really optimal as the LAB survey does not come
  // in an equal area map, but good enough for government
  // use [i.e., that's what they do in the nh tool])
  //

  double nhsum=0.;
  double weightsum=0.;
  for (long dx=-sampling;dx<=sampling;dx++) {
    for (long dy=-sampling;dy<=sampling;dy++) {

      long xx=xpix+dx;
      long yy=ypix+dy;

      // we avoid wrap around in y, but allow it in x
      if (yy>=0 && yy<nhmap->naxes[1]) {
	if (xx<0) { xx=nhmap->naxes[0]+xx; }
	if (xx>=nhmap->naxes[0]) {xx=xx-nhmap->naxes[0];}

	if (verbosity>1) {
	  fprintf(stderr," %04li/%04li:",xx,yy);
	}

	// get position of pixel
	pixcrd[0]=xx+1;
	pixcrd[1]=yy+1;
	pixcrd[2]=0;
	int ret=wcsp2s(nhmap->wcs,1,3,pixcrd,imgcrd,phi,theta,world,status);
	if (ret != 0 ) {
	  warnx("problem with wcsp2s: %i",ret);
	  return -1.;
	}
	double pxllong=world[nhmap->wcs->lng];
	double pxllat=world[nhmap->wcs->lat];
	
	if (verbosity>0){

	  double pxllongrad=pxllong;
	  double pxllatrad=pxllat;

	  fprintf(stderr," %7.2f/%7.2f ",pxllong,pxllat);

	  // convert to equatorial
	  double rra,rrd;
	  atGaltoJ2000(pxllongrad,pxllatrad,&rra,&rrd);

	  char *rrst=posstring(rra,HOURPOS);
	  char *ddst=posstring(rrd,LATPOS);
	  fprintf(stderr," %s/%s ",rrst,ddst);
	  free(rrst);
	  free(ddst);
	}

	// angular distance
	double dist;
	if (sampling>0 || verbosity>0) {
	  dist=angdist(lii,bii,pxllong,pxllat);

	  if (verbosity>0) {
	    fprintf(stderr," %6.4f",dist);
	  }
	}

	// only calculate NH if not in cache
	long cachpos=xx+nhmap->naxes[0]*yy;
	if (  nhmap->nhcache[cachpos] == 0.) {
	  double nh=0.;
	  for (int z=zmin; z<zmax; z++) {
	    nh+=get21cm(xx,yy,z);
	  }
	  // conversion to NH
	  nh=1.82e15*nh*nhmap->dv;

	  // put in cache
	  nhmap->nhcache[cachpos]=nh;
	}

	// perform weighted sum
	double weight=1.;

	if ( sampling>0 ) {
	  double sigma=0.3;
	  weight=exp(-dist*dist/(2.*sigma*sigma));
	}

	nhsum+=weight*nhmap->nhcache[cachpos];
	weightsum+=weight;

	if (verbosity>0) {
	  fprintf(stderr," %9.2e %7.3f \n",
		  nhmap->nhcache[cachpos],weight);
	}

      }
    }
  }

  double nh=nhsum/weightsum;

  if (verbosity>0) {
    fprintf(stderr," N_H: %e\n",nh);
  }

  return nh;
}
