/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include <float.h>
#include "common.h"


/** Random number generator. */
static double(*static_rndgen)(int* const)=NULL;


void setSimputARF(SimputCtlg* const cat, struct ARF* const arf)
{
  cat->arf=arf;
}


void setSimputConstantARF(SimputCtlg* const cat, double low_energy, double high_energy, double eff_area, char* telescope, int* const status)
{
  cat->arf=getconstantARF(low_energy, high_energy, eff_area, telescope, status);
}



void loadSimputARF(SimputCtlg* const cat, char* const filename, int* const status)
{
  cat->arf=loadARF(filename, status);
}


/** Use the C rand() function to determine a random number between 0
    and 1. */
static double getCRand(int* const status) 
{
  double r=(double)rand()/((double)RAND_MAX+1.0);
  assert(r<1.0);
  return(r);

  // Status variable is not needed.
  (void)(*status);
}


void setSimputRndGen(double(*rndgen)(int* const))
{
  static_rndgen=rndgen;
}


/** Determine a random number between 0 and 1 with the specified
    random number generator. */
static inline double getRndNum(int* const status)
{
  // Check if a random number generator has been set.
  if (NULL==static_rndgen) {
    // If not use the C rand() generator as default.
    SIMPUT_WARNING("use C rand() as default since no random number generator "
		   "is specified");

    // Initialize random seed.
    srand(time(NULL));

    setSimputRndGen(getCRand);
  }

  // Obtain a random number.
  return(static_rndgen(status));
}


SimputSrc* getSimputSrc(SimputCtlg* const cf,
			const long row,
			int* const status)
{
  // Maximum number of sources in the cache.
  const long maxsrcs=1000000; 

  // Check if the source catalog contains a source buffer.
  if (NULL==cf->srcbuff) {
    cf->srcbuff=newSimputSrcBuffer(status);
    CHECK_STATUS_RET(*status, NULL);
  }

  // Convert the void* pointer to the source buffer into the right
  // format.
  struct SimputSrcBuffer* sb=(struct SimputSrcBuffer*)cf->srcbuff;

  // Allocate memory for the cache.
  if (NULL==sb->rowmap) {
    sb->rowmap=(long*)malloc(cf->nentries*sizeof(long));
    CHECK_NULL_RET(sb->rowmap, *status, 
		   "memory allocation for row map failed", NULL);
    long jj;
    for (jj=0; jj<cf->nentries; jj++) {
      sb->rowmap[jj]=-1;
    }
  }

  // Check if the cache already exists or whether this routine is 
  // called for the first time.
  if (NULL==sb->srcs) {
    sb->srcs=(SimputSrc**)malloc(maxsrcs*sizeof(SimputSrc*));
    CHECK_NULL_RET(sb->srcs, *status,
		   "memory allocation for source cache failed", NULL);
  }
  if (NULL==sb->rownums) {
    sb->rownums=(long*)malloc(maxsrcs*sizeof(long));
    CHECK_NULL_RET(sb->rownums, *status,
		   "memory allocation for array of row numbers failed", NULL);
  }

  // Check if the specified source index (row number) is valid.
  if ((row<=0) || (row>cf->nentries)) {
    SIMPUT_ERROR("invalid row number");
    return(NULL);
  }

  // Check if the requested row is already available in the cache.
  if (sb->rowmap[row-1]>0) {
    return(sb->srcs[sb->rowmap[row-1]]);
  }

  // The requested source is not contained in the cache.
  // Therefore we must load it from the FITS file.
  
  // Check if the cache is already full.
  if (sb->nsrcs<maxsrcs) {
    sb->csrc = sb->nsrcs;
    sb->nsrcs++;
  } else {
    sb->csrc++;
    if (sb->csrc>=maxsrcs) {
      sb->csrc=0;
    }
    // Destroy the source that is currently stored at this place 
    // in the cache.
    freeSimputSrc(&(sb->srcs[sb->csrc]));
    sb->rowmap[sb->rownums[sb->csrc]-1] = -1;
    sb->rownums[sb->csrc] = 0;
  }

  // Load the source from the FITS file.
  sb->srcs[sb->csrc]=loadSimputSrc(cf, row, status);
  CHECK_STATUS_RET(*status, sb->srcs[sb->csrc]);
  sb->rownums[sb->csrc]=row;
  sb->rowmap[row-1]=sb->csrc;

  return(sb->srcs[sb->csrc]);
}


static void getSrcTimeRef(SimputCtlg* const cat,
			  const SimputSrc* const src,
			  char* const timeref)
{
  // Initialize with an empty string.
  strcpy(timeref, "");

  // Determine the reference to the timing extension.
  if (NULL==src->timing) {
    strcpy(timeref, "");
  } else if ((0==strlen(src->timing)) || 
	     (0==strcmp(src->timing, "NULL")) ||
	     (0==strcmp(src->timing, " "))) {
    strcpy(timeref, "");
  } else {
    if ('['==src->timing[0]) {
      strcpy(timeref, cat->filepath);
      strcat(timeref, cat->filename);
    } else {
      if ('/'!=src->timing[0]) {
	strcpy(timeref, cat->filepath);
      } else {
	strcpy(timeref, "");
      }
    }
    strcat(timeref, src->timing);
  }
}


static void gauss_rndgen(double* const x, double* const y, int* const status)
{
  double sqrt_2rho=sqrt(-log(getRndNum(status))*2.);
  CHECK_STATUS_VOID(*status);
  double phi=getRndNum(status)*2.*M_PI;
  CHECK_STATUS_VOID(*status);

  *x=sqrt_2rho * cos(phi);
  *y=sqrt_2rho * sin(phi);
}


static SimputPSD* getSimputPSD(SimputCtlg* const cat,
			       char* const filename,
			       int* const status)
{
  const long maxpsds=100;

  // Check if the source catalog contains a PSD buffer.
  if (NULL==cat->psdbuff) {
    cat->psdbuff=newSimputPSDBuffer(status);
    CHECK_STATUS_RET(*status, NULL);
  }

  // Convert the void* pointer to the PSD buffer into the right
  // format.
  struct SimputPSDBuffer* sb=(struct SimputPSDBuffer*)cat->psdbuff;

  // In case there are no PSDs available at all, allocate 
  // memory for the array (storage for PSDs).
  if (NULL==sb->psds) {
    sb->psds=(SimputPSD**)malloc(maxpsds*sizeof(SimputPSD*));
    CHECK_NULL_RET(sb->psds, *status, 
		   "memory allocation for PSD buffer failed", NULL);
  }

  // Search if the requested PSD is available in the storage.
  long ii;
  for (ii=0; ii<sb->npsds; ii++) {
    // Check if the PSD is equivalent to the requested one.
    if (0==strcmp(sb->psds[ii]->fileref, filename)) {
      // If yes, return the PSD.
      return(sb->psds[ii]);
    }
  }

  // The requested PSD is not contained in the storage.
  // Therefore we must load it from the specified location.
  if (sb->npsds>=maxpsds) {
    SIMPUT_ERROR("too many PSDs in the internal storage");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // Load the PSD from the file.
  sb->psds[sb->npsds]=loadSimputPSD(filename, status);
  CHECK_STATUS_RET(*status, sb->psds[sb->npsds]);
  sb->npsds++;

  return(sb->psds[sb->npsds-1]);
}


static inline double getLCTime(const SimputLC* const lc,
			       const long kk, 
			       const long long nperiods,
			       const double mjdref)
{
  if (NULL!=lc->time) {
    // Non-periodic light curve.
    return(lc->time[kk]+lc->timezero+(lc->mjdref-mjdref)*24.*3600.);
  } else {
    // Periodic light curve.
    long double phase=lc->phase[kk]-lc->phase0+nperiods;
    if (fabs(lc->dperiod)<1.e-20) {
      return(phase*lc->period);
    } else {
      return((double)((expm1l(phase*lc->dperiod))*lc->period/lc->dperiod
		      +lc->timezero+(lc->mjdref-mjdref)*24.*3600.));
    }
  }
}


/** Determine the index of the bin of the light curve that corresponds
    to the specified time. For periodic light curves the function
    stores the number of periods since the specified origin of the
    light curve in the parameter nperiods. */
static inline long getLCBin(const SimputLC* const lc, 
			    const double time, 
			    const double mjdref,
			    long long* nperiods, 
			    int* const status)
{
  // Check if the value of MJDREF is negative. This indicates
  // that the spectrum referred to in the first bin of the light curve
  // is requested in order to calculate the source count rate.
  if (mjdref<0.) {
    return(0);
  }

  // Check if the light curve is periodic or not.
  if (NULL!=lc->time) {
    // Non-periodic light curve.
    *nperiods=0;

    // Check if the requested time is within the covered interval.
    double t0=getLCTime(lc, 0, 0, mjdref);
    double t1=getLCTime(lc, lc->nentries-1, 0, mjdref);
    if ((time<t0) || (time>=t1)) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "requested time (%lf MJD) is outside the "
	      "interval covered by the light curve '%s' (%lf to %lf MJD)",
	      time/24./3600.+mjdref, lc->fileref, t0/24./3600., t1/24./3600.);
      SIMPUT_ERROR(msg);
      *status=EXIT_FAILURE;
      return(0);
    }

  } else {
    // Periodic light curve.
    // Make a first guess on the number of passed periods.
    double dt=time-getLCTime(lc, 0, 0, mjdref);
    double phase;
    if (fabs(lc->dperiod)<1.e-20) {
      phase=lc->phase0+dt/lc->period;
    } else {
      phase=lc->phase0+log(1.+dt*lc->dperiod/lc->period)/lc->dperiod;
    }
    *nperiods=(long long)phase;

    // Correct the first guess such that the requested time lies within the
    // covered period. Deviations with respect to the first guess can 
    // introduced by a \dot{P} (DPERIOD).
    while (getLCTime(lc, 0, (*nperiods)+1, mjdref) <= time) {
      (*nperiods)++;
    }
    while (getLCTime(lc, 0, *nperiods, mjdref) > time) {
      (*nperiods)--;
    }
  }

  // Determine the index of the light curve bin (using binary search).
  long lower=0, upper=lc->nentries-2, mid;
  while (upper>lower) {
    mid=(lower+upper)/2;
    if (getLCTime(lc, mid+1, *nperiods, mjdref) < time) {
      lower=mid+1;
    } else {
      upper=mid;
    }
  }

  return(lower);
}


static SimputLC* getSimputLC(SimputCtlg* const cat, 
			     const SimputSrc* const src,
			     char* const filename, 
			     const double prevtime, 
			     const double mjdref, 
			     int* const status)
{
  SimputLC* lc=NULL;

  // Check if the SimputLC is contained in the internal cache.
  const long maxlcs=1;

  // Check if the source catalog contains a light curve buffer.
  if (NULL==cat->lcbuff) {
    cat->lcbuff=newSimputLCBuffer(status);
    CHECK_STATUS_RET(*status, NULL);
  }

  // Convert the void* pointer to the light curve buffer into 
  // the right format.
  struct SimputLCBuffer* lb=(struct SimputLCBuffer*)cat->lcbuff;

  // In case there are no light curves available at all, allocate 
  // memory for the array (storage for images).
  if (NULL==lb->lcs) {
    lb->lcs=(SimputLC**)malloc(maxlcs*sizeof(SimputLC*));
    CHECK_NULL_RET(lb->lcs, *status, 
		   "memory allocation for light curves failed", NULL);
  }

  // Search if the requested light curve is available in the storage.
  long ii;
  for (ii=0; ii<lb->nlcs; ii++) {
    // Check if the light curve is equivalent to the requested one.
    if (0==strcmp(lb->lcs[ii]->fileref, filename)) {
      // For a light curve created from a PSD, we also have to check,
      // whether this is the source associated to this light curve.
      if (lb->lcs[ii]->src_id>0) {
	// Check if the SRC_IDs agree.
	if (lb->lcs[ii]->src_id==src->src_id) {
	  // We have a light curve which has been produced from a PSD.
	  // Check if the requested time is covered by the light curve.
	  if (prevtime<getLCTime(lb->lcs[ii], lb->lcs[ii]->nentries-1, 
				 0, mjdref)) {
	    return(lb->lcs[ii]);
	  }
	  // If not, we have to produce a new light curve from the PSD.
	}
      } else {
	// This light curve is loaded from a file and can be re-used
	// for different sources.
	return(lb->lcs[ii]);
      }
    }
  }


  // If the LC is not contained in the cache, load it either from 
  // a file or create it from a SimputPSD.
  int timetype=getSimputExtType(cat, filename, status);
  CHECK_STATUS_RET(*status, lc);

  if (EXTTYPE_LC==timetype) {
    // Load directly from file.
    lc=loadSimputLC(filename, status);
    CHECK_STATUS_RET(*status, lc);

  } else {
    // Create the SimputLC from a SimputPSD.
    
    // Buffer for Fourier transform.
    float* power=NULL;
    double *fftw_in=NULL, *fftw_out=NULL;

    // Length of the PSD for the FFT, which is obtained by 
    // interpolation of the input PSD.
    const long psdlen=100000000;

    do { // Error handling loop.
  
      SimputPSD* psd=getSimputPSD(cat, filename, status);
      CHECK_STATUS_BREAK(*status);
      
      // Get an empty SimputLC data structure.
      lc=newSimputLC(status);
      CHECK_STATUS_BREAK(*status);
      
      // Set the MJDREF.
      lc->mjdref=mjdref;

      // Allocate memory for the light curve.
      lc->nentries=2*psdlen;
      lc->time    =(double*)malloc(lc->nentries*sizeof(double));
      CHECK_NULL_BREAK(lc->time, *status, 
		       "memory allocation for K&R light curve failed");
      lc->flux    =(float*)malloc(lc->nentries*sizeof(float));
      CHECK_NULL_BREAK(lc->flux, *status, 
		       "memory allocation for K&R light curve failed");

      // Set the time bins of the light curve.
      lc->timezero=prevtime;
      long ii;
      double dt=1./(2.*psd->frequency[psd->nentries-1]);
      for (ii=0; ii<lc->nentries; ii++) {
	lc->time[ii]=ii*dt;
      }

      // Interpolate the PSD to a uniform frequency grid.
      // The PSD is given in Miyamoto normalization. In order to get the RMS
      // right, we have to multiply each bin with df (delta frequency).
      power=(float*)malloc(psdlen*sizeof(float));
      CHECK_NULL_BREAK(power, *status, 
		       "memory allocation for PSD buffer failed");

      // Check that the PSD is positive.
      long jj;
      for (jj=0; jj<psd->nentries; jj++) {
	if (psd->power[jj]<=0.0) {
	  SIMPUT_ERROR("PSD may only have positive entries");
	  *status=EXIT_FAILURE;
	  break;
	}
      }
      CHECK_STATUS_RET(*status, lc);

      float delta_f=psd->frequency[psd->nentries-1]/psdlen;
      jj=0;
      for (ii=0; ii<psdlen; ii++) {
	float frequency=(ii+1)*delta_f;
	while((frequency>psd->frequency[jj]) &&
	      (jj<psd->nentries-1)) {
	  jj++;
	}
	if (jj==0) {
	  power[ii]=0.;
	  /* frequency/psd->frequency[jj]* 
	     psd->power[jj]* 
	     delta_f; */
	} else {
	  power[ii]=
	    (psd->power[jj-1]+
	     (frequency-psd->frequency[jj-1])/
	     (psd->frequency[jj]-psd->frequency[jj-1])*
	     (psd->power[jj]-psd->power[jj-1]))*
	    delta_f;
	}
      }
    
      // Allocate the data structures required by the fftw routines.
      fftw_in =(double*)fftw_malloc(sizeof(double)*lc->nentries);
      CHECK_NULL_BREAK(fftw_in, *status, "memory allocation for fftw "
		       "data structure (fftw_in) failed");
      fftw_out=(double*)fftw_malloc(sizeof(double)*lc->nentries);
      CHECK_NULL_BREAK(fftw_out, *status, "memory allocation for fftw "
		       "data structure (fftw_out) failed");

      // Apply the algorithm introduced by Timmer & Koenig (1995).
      double randr, randi;
      lc->fluxscal=1.; // Set Fluxscal to 1.
      gauss_rndgen(&randr, &randi, status);
      CHECK_STATUS_RET(*status, lc);
      fftw_in[0]     =1.;
      fftw_in[psdlen]=randi*sqrt(power[psdlen-1]);
      for (ii=1; ii<psdlen; ii++) {
	gauss_rndgen(&randr, &randi, status);
	CHECK_STATUS_RET(*status, lc);
	REAL(fftw_in, ii)              =randr*0.5*sqrt(power[ii-1]);
	IMAG(fftw_in, ii, lc->nentries)=randi*0.5*sqrt(power[ii-1]);
      }

      // Perform the inverse Fourier transformation.
      fftw_plan iplan=fftw_plan_r2r_1d(lc->nentries, fftw_in, fftw_out, 
				       FFTW_HC2R, FFTW_ESTIMATE);
      fftw_execute(iplan);
      fftw_destroy_plan(iplan);
      
      // Determine the normalized rates from the FFT.
      for (ii=0; ii<lc->nentries; ii++) {
	lc->flux[ii]=(float)fftw_out[ii] /* *requ_rms/act_rms */;

	// Avoid negative fluxes (no physical meaning):
	if (lc->flux[ii]<0.) { 
	  lc->flux[ii]=0.; 
	}
      }

      // This light curve has been generated for this
      // particular source.
      lc->src_id=src->src_id;

      // Store the file reference to the timing extension for later 
      // comparisons.
      lc->fileref=
	(char*)malloc((strlen(filename)+1)*sizeof(char));
      CHECK_NULL_RET(lc->fileref, *status, 
		     "memory allocation for file reference failed", 
		     lc);
      strcpy(lc->fileref, filename);

    } while(0); // END of error handling loop.
    
    // Release allocated memory.
    if (NULL!=power) free(power);
    if (NULL!=fftw_in) fftw_free(fftw_in);
    if (NULL!=fftw_out) fftw_free(fftw_out);
    
  }

  // Check if there is still space left in the internal cache.
  if (lb->nlcs<maxlcs) {
    lb->clc=lb->nlcs;
    lb->nlcs++;
  } else {
    lb->clc++;
    if (lb->clc>=maxlcs) {
      lb->clc=0;
    }
    // Release the SimputLC that is currently stored at this place 
    // in the cache.
    freeSimputLC(&(lb->lcs[lb->clc]));
  }

  // Store the SimputLC in the internal cache.
  lb->lcs[lb->clc]=lc;

  return(lb->lcs[lb->clc]);
}


void getSimputSrcSpecRef(SimputCtlg* const cat,
			 const SimputSrc* const src,
			 const double prevtime,
			 const double mjdref,
			 char* const specref,
			 int* const status)
{
  // Initialize with an empty string.
  strcpy(specref, "");

  // Determine the timing extension.
  char timeref[SIMPUT_MAXSTR];
  getSrcTimeRef(cat, src, timeref);
  CHECK_STATUS_VOID(*status);

  int timetype=getSimputExtType(cat, timeref, status);
  CHECK_STATUS_VOID(*status);

  if (EXTTYPE_LC==timetype) {
    // Get the respective light curve.
    SimputLC* lc=getSimputLC(cat, src, timeref, prevtime, mjdref, status);
    CHECK_STATUS_VOID(*status);

    // Check if there is a spectrum column in the light curve.
    if (NULL!=lc->spectrum) {
      // Determine the current light curve bin.
      long long nperiods;
      long bin=getLCBin(lc, prevtime, mjdref, &nperiods, status);
      CHECK_STATUS_VOID(*status);

      // Check if a spectrum is defined in this light curve bin.
      if ((0==strcmp(lc->spectrum[bin], "NULL")) || 
	  (0==strcmp(lc->spectrum[bin], " ")) ||
	  (0==strlen(lc->spectrum[bin]))) {
	SIMPUT_ERROR("in the current implementation light curves "
		     "must not contain blank entries in a given "
		     "spectrum column");
	*status=EXIT_FAILURE;
	return;
      }

      // Copy the reference to the spectrum.
      strcpy(specref, lc->spectrum[bin]);

      // Determine the relative location with respect to the 
      // location of the light curve.
      if ('['==specref[0]) {
	char* firstbrack=strchr(timeref, '[');
	strcpy(firstbrack, specref);
	strcpy(specref, timeref);
      } else if (('/'!=specref[0])&&(NULL!=strchr(timeref, '/'))) {
	char* lastslash=strrchr(timeref, '/');
	lastslash++;
	strcpy(lastslash, specref);
	strcpy(specref, timeref);	
      }

      return;
    }
  }

  // If no light curve extension with a spectrum column is given, 
  // determine the spectrum reference directly from the source 
  // description.
  if (NULL!=src->spectrum) {   
    // Check if this is a valid HDU reference.
    if ((0==strcmp(src->spectrum, "NULL")) || 
	(0==strcmp(src->spectrum, " ")) ||
	(0==strlen(src->spectrum))) {
      return;
    }

    strcpy(specref, src->spectrum);
    
    // Set path and file name if missing.
    if ('['==specref[0]) {
      char buffer[SIMPUT_MAXSTR];
      strcpy(buffer, cat->filepath);
      strcat(buffer, cat->filename);
      strcat(buffer, specref);
      strcpy(specref, buffer);
    } else {
      if ('/'!=specref[0]) {
	char buffer[SIMPUT_MAXSTR];
	strcpy(buffer, cat->filepath);
	strcat(buffer, specref);
	strcpy(specref, buffer);
      }
    }
  }
}


static void getSrcImagRef(SimputCtlg* const cat,
			  const SimputSrc* const src,
			  const double prevtime,
			  const double mjdref,
			  char* const imagref,
			  int* const status)
{
  // Initialize with an empty string.
  strcpy(imagref, "");

  // Determine the timing extension.
  char timeref[SIMPUT_MAXSTR];
  getSrcTimeRef(cat, src, timeref);
  CHECK_STATUS_VOID(*status);

  int timetype=getSimputExtType(cat, timeref, status);
  CHECK_STATUS_VOID(*status);

  if (EXTTYPE_LC==timetype) {
    // Get the respective light curve.
    SimputLC* lc=getSimputLC(cat, src, timeref, prevtime, mjdref, status);
    CHECK_STATUS_VOID(*status);

    // Check if there is an image column in the light curve.
    if (NULL!=lc->image) {
      // Determine the current light curve bin.
      long long nperiods;
      long bin=getLCBin(lc, prevtime, mjdref, &nperiods, status);
      CHECK_STATUS_VOID(*status);

      // Check if an image is defined in this light curve bin.
      if ((0==strcmp(lc->image[bin], "NULL")) ||
	  (0==strcmp(lc->image[bin], " ")) ||
	  (0==strlen(lc->image[bin]))) {
	SIMPUT_ERROR("in the current implementation light curves "
		     "must not contain blank entries in a given "
		     "image column");
	*status=EXIT_FAILURE;
	return;
      }

      // Copy the reference to the image.
      strcpy(imagref, lc->image[bin]);

      // Determine the relative location with respect to the 
      // location of the light curve.
      if ('['==imagref[0]) {
	char* firstbrack=strchr(timeref, '[');
	strcpy(firstbrack, imagref);
	strcpy(imagref, timeref);
      } else if (('/'!=imagref[0])&&(NULL!=strchr(timeref, '/'))) {
	char* lastslash=strrchr(timeref, '/');
	lastslash++;
	strcpy(lastslash, imagref);
	strcpy(imagref, timeref);
      }

      return;
    }
  }

  // If no light curve extension with an image column is given, 
  // determine the spectrum reference directly from the source 
  // description.
  if (NULL!=src->image) {
    // Check if this is a valid HDU reference.
    if ((0==strcmp(src->image, "NULL")) ||
	(0==strcmp(src->image, " ")) ||
	(0==strlen(src->image))) {
      return;
    }
    
    strcpy(imagref, src->image);

    // Set path and file name if missing.
    if ('['==imagref[0]) {
      char buffer[SIMPUT_MAXSTR];
      strcpy(buffer, cat->filepath);
      strcat(buffer, cat->filename);
      strcat(buffer, imagref);
      strcpy(imagref, buffer);
    } else {
      if ('/'!=imagref[0]) {
	char buffer[SIMPUT_MAXSTR];
	strcpy(buffer, cat->filepath);
	strcat(buffer, imagref);
	strcpy(imagref, buffer);
      }
    }
  }
}


static SimputMIdpSpec* getSimputMIdpSpec(SimputCtlg* const cat,
					 const char* const filename,
					 int* const status)
{
  // Search if the spectrum is available in the buffer.
  SimputMIdpSpec* spec=
    searchSimputMIdpSpecBuffer(cat->midpspecbuff, filename);
  if (NULL!=spec) {
    return(spec);
  }

  // The required spectrum is not contained in the buffer.
  // Therefore it must be loaded from the specified location.

  // Load the mission-independent spectrum.
  spec=loadSimputMIdpSpec(filename, status);
  CHECK_STATUS_RET(*status, spec);

  // Insert the spectrum into the buffer.
  insertSimputMIdpSpecBuffer(&(cat->midpspecbuff), spec, status);
  CHECK_STATUS_RET(*status, spec);

  return(spec);
}


SimputMIdpSpec* getSimputSrcMIdpSpec(SimputCtlg* const cat,
				     const SimputSrc* const src,
				     const double prevtime,
				     const double mjdref,
				     int* const status)
{
  char specref[SIMPUT_MAXSTR];
  getSimputSrcSpecRef(cat, src, prevtime, mjdref, specref, status);
  CHECK_STATUS_RET(*status, 0);

  int spectype=getSimputExtType(cat, specref, status);
  CHECK_STATUS_RET(*status, 0);

  if (EXTTYPE_MIDPSPEC==spectype) {
    // Determine the spectrum.
    SimputMIdpSpec* spec=getSimputMIdpSpec(cat, specref, status);
    CHECK_STATUS_RET(*status, 0);

    return(spec);

  } else {
    SIMPUT_ERROR("source does not refer to a mission-independent spectrum");
    *status=EXIT_FAILURE;
    return(NULL);
  }
}


static inline void getMIdpSpecEbounds(const SimputMIdpSpec* const spec,
				      const long idx,
				      float* const emin, 
				      float* const emax)
{
  // Determine the lower boundary.
  if (idx>0) {
    *emin=0.5*(spec->energy[idx]+spec->energy[idx-1]);
  } else {
    *emin=spec->energy[idx];
  }

  // Determine the upper boundary.
  if (idx<spec->nentries-1) {
    *emax=0.5*(spec->energy[idx+1]+spec->energy[idx]);
  } else {
    *emax=spec->energy[idx];
  }
}


/** Convolve the given mission-independent spectrum with the
    instrument ARF. The product of this process is the spectral
    probability distribution binned to the energy grid of the ARF. */
static SimputSpec* convSimputMIdpSpecWithARF(SimputCtlg* const cat, 
					     SimputMIdpSpec* const midpspec, 
					     int* const status)
{
  SimputSpec* spec=NULL;

  // Check if the ARF is defined.
  CHECK_NULL_RET(cat->arf, *status, "instrument ARF undefined", spec);

  // Allocate memory.
  spec=newSimputSpec(status);
  CHECK_STATUS_RET(*status, spec);
  spec->distribution=
    (double*)malloc(cat->arf->NumberEnergyBins*sizeof(double));
  CHECK_NULL_RET(spec->distribution, *status,
		 "memory allocation for spectral distribution failed", spec);

  // Loop over all bins of the ARF.
  long ii, jj=0;
  int warning_printed=0; // Flag whether warning has been printed.
  for (ii=0; ii<cat->arf->NumberEnergyBins; ii++) {
    // Initialize with 0.
    spec->distribution[ii]=0.;

    // Lower boundary of the current bin.
    float lo=cat->arf->LowEnergy[ii];

    // Loop over all spectral points within the ARF bin.
    int finished=0;
    do {
      // Determine the next spectral point.
      float spec_emin=0., spec_emax=0.;
      for ( ; jj<midpspec->nentries; jj++) {
	getMIdpSpecEbounds(midpspec, jj, &spec_emin, &spec_emax);
	if (spec_emax>lo) break;
      }
      
      // Check special cases.
      if ((0==jj) && (spec_emin>cat->arf->LowEnergy[ii])) {
	if (0==warning_printed) {
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "the spectrum '%s' does not cover the "
		  "full energy range of the ARF", midpspec->fileref);
	  SIMPUT_WARNING(msg);
	  warning_printed=1;
	}
	if (spec_emin>cat->arf->HighEnergy[ii]) break;

      } else if (jj==midpspec->nentries) {
	if (0==warning_printed) {
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "the spectrum '%s' does not cover the "
		  "full energy range of the ARF", midpspec->fileref);
	  SIMPUT_WARNING(msg);
	  warning_printed=1;
	}
	break;
      }
      
      // Upper boundary of the current bin.
      float hi;
      if (spec_emax<=cat->arf->HighEnergy[ii]) {
	hi=spec_emax;
      } else {
	hi=cat->arf->HighEnergy[ii];
	finished=1;
      }

      // Add to the spectral probability density.
      spec->distribution[ii]+=
	(hi-lo)*cat->arf->EffArea[ii]*midpspec->fluxdensity[jj];
      
      // Increase the lower boundary.
      lo=hi;

    } while (0==finished);

    // Create the spectral distribution function 
    // normalized to the total photon number [photons]. 
    if (ii>0) {
      spec->distribution[ii]+=spec->distribution[ii-1];
    }
  } // Loop over all ARF bins.


  // Copy the file reference to the spectrum for later comparisons.
  spec->fileref=
    (char*)malloc((strlen(midpspec->fileref)+1)*sizeof(char));
  CHECK_NULL_RET(spec->fileref, *status, 
		 "memory allocation for file reference failed", 
		 spec);
  strcpy(spec->fileref, midpspec->fileref);

  return(spec);
}


static SimputSpec* getSimputSpec(SimputCtlg* const cat,
				 const char* const filename,
				 int* const status)
{
  // Search if the spectrum is available in the buffer.
  SimputSpec* spec=searchSimputSpecBuffer(cat->specbuff, filename);
  if (NULL!=spec) {
    return(spec);
  }

  // The required spectrum is not contained in the buffer.
  // Therefore we must determine it from the referred mission-
  // independent spectrum and store it in the buffer.

  // Obtain the mission-independent spectrum.
  SimputMIdpSpec* midpspec=getSimputMIdpSpec(cat, filename, status);
  CHECK_STATUS_RET(*status, NULL);

  // Convolve it with the ARF.
  spec=convSimputMIdpSpecWithARF(cat, midpspec, status);
  CHECK_STATUS_RET(*status, spec);

  // Insert the spectrum into the buffer.
  insertSimputSpecBuffer(&(cat->specbuff), spec, status);
  CHECK_STATUS_RET(*status, spec);

  return(spec);
}


static inline double rndexp(const double avgdist, int* const status)
{
  assert(avgdist>0.);

  double rand;
  do {
    rand=getRndNum(status);
    CHECK_STATUS_RET(*status, 0.);
    assert(rand>=0.);
  } while (rand==0.);

  return(-log(rand)*avgdist);
}


/** Return the requested image. Keeps a certain number of images in an
    internal storage. If the requested image is not located in the
    internal storage, it is loaded from the reference given in the
    source catalog. */
static SimputImg* getSimputImg(SimputCtlg* const cat,
			       char* const filename,
			       int* const status)
{
  const long maximgs=200;

  // Check if the source catalog contains an image buffer.
  if (NULL==cat->imgbuff) {
    cat->imgbuff=newSimputImgBuffer(status);
    CHECK_STATUS_RET(*status, NULL);
  }

  // Convert the void* pointer to the image buffer into the right
  // format.
  struct SimputImgBuffer* sb=(struct SimputImgBuffer*)cat->imgbuff;

  // In case there are no images available at all, allocate 
  // memory for the array (storage for images).
  if (NULL==sb->imgs) {
    sb->imgs=(SimputImg**)malloc(maximgs*sizeof(SimputImg*));
    CHECK_NULL_RET(sb->imgs, *status, 
		   "memory allocation for images failed", NULL);
  }

  // Search if the requested image is available in the storage.
  long ii;
  for (ii=0; ii<sb->nimgs; ii++) {
    // Check if the image is equivalent to the requested one.
    if (0==strcmp(sb->imgs[ii]->fileref, filename)) {
      // If yes, return the image.
      return(sb->imgs[ii]);
    }
  }

  // The requested image is not contained in the storage.
  // Therefore we must load it from the specified location.
  if (sb->nimgs>=maximgs) {
    SIMPUT_ERROR("too many images in the internal storage");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // Load the image from the file.
  sb->imgs[sb->nimgs]=loadSimputImg(filename, status);
  CHECK_STATUS_RET(*status, sb->imgs[sb->nimgs]);
  sb->nimgs++;
  
  return(sb->imgs[sb->nimgs-1]);
}


void simput_s2p(struct wcsprm* const wcs,
		double* px,
		double* py,
		const double sx,
		const double sy,
		int* const status)
{
	// convert to degree
  double world[2] = { sx*180./M_PI, sy*180./M_PI };
  double pixcrd[2], imgcrd[2];
  double phi, theta;

  // If CUNIT is set to 'degree', change this to 'deg'.
  // Otherwise the WCSlib will not work properly.
  if (0==strcmp(wcs->cunit[0], "degree  ")) {
    strcpy(wcs->cunit[0], "deg");
  }
  if (0==strcmp(wcs->cunit[1], "degree  ")) {
    strcpy(wcs->cunit[1], "deg");
  }

// Perform the transform using WCSlib.
  int retval=wcss2p(wcs, 1, 2, world, &phi, &theta,  imgcrd, pixcrd, status);
  CHECK_STATUS_VOID(*status);
  if (0!=retval) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("WCS transformation failed");
    return;
  }

  *px=pixcrd[0];
  *py=pixcrd[1];
}


inline static void p2s(struct wcsprm* const wcs,
		const double px, 
		const double py,
		double* sx, 
		double* sy,
		int* const status)
{
  double pixcrd[2]={ px, py };
  double imgcrd[2], world[2];
  double phi, theta;
  
  // If CUNIT is set to 'degree', change this to 'deg'.
  // Otherwise the WCSlib will not work properly.
  if (0==strcmp(wcs->cunit[0], "degree  ")) {
    strcpy(wcs->cunit[0], "deg");
  }
  if (0==strcmp(wcs->cunit[1], "degree  ")) {
    strcpy(wcs->cunit[1], "deg");
  }

  // Perform the transform using WCSlib.
  int retval=wcsp2s(wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, status);
  CHECK_STATUS_VOID(*status);
  if (0!=retval) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("WCS transformation failed");
    return;
  }
  
  // Convert to [rad].
  *sx=world[0]*M_PI/180.;
  *sy=world[1]*M_PI/180.;
}

// make the p2s function publicly available
void simput_p2s(struct wcsprm* const wcs,
		const double px,
		const double py,
		double* sx,
		double* sy,
		int* const status)
{
	p2s(wcs,px,py,sx,sy,status);
}

/** Return the requested photon lists. Keeps a certain number of
    photon lists in an internal storage. If the requested photon list
    is not located in the internal storage, it is obtained from the
    reference given in the source catalog. */
static SimputPhList* getSimputPhList(SimputCtlg* const cat,
				     char* const filename,
				     int* const status)
{
  const long maxphls=100;

  // Check if the source catalog contains a photon list buffer.
  if (NULL==cat->phlistbuff) {
    cat->phlistbuff=newSimputPhListBuffer(status);
    CHECK_STATUS_RET(*status, NULL);
  }

  // Convert the void* pointer to the photon list buffer into the right
  // format.
  struct SimputPhListBuffer* pb=(struct SimputPhListBuffer*)cat->phlistbuff;

  // In case there are no photon lists available at all, allocate 
  // memory for the array (storage for photon lists).
  if (NULL==pb->phls) {
    pb->phls=(SimputPhList**)malloc(maxphls*sizeof(SimputPhList*));
    CHECK_NULL_RET(pb->phls, *status, 
		   "memory allocation for photon lists failed", NULL);
  }

  // Search if the requested photon list is available in the storage.
  long ii;
  for (ii=0; ii<pb->nphls; ii++) {
    // Check if the photon list is equivalent to the requested one.
    if (0==strcmp(pb->phls[ii]->fileref, filename)) {
      // If yes, return the photon list.
      return(pb->phls[ii]);
    }
  }

  // The requested photon list is not contained in the storage.
  // Therefore we must open it from the specified location.
  if (pb->nphls>=maxphls) {
    SIMPUT_ERROR("too many photon lists in the internal storage");
    *status=EXIT_FAILURE;
    return(NULL);
  }

  // Open the photon list from the file.
  pb->phls[pb->nphls]=openSimputPhList(filename, READONLY, status);
  CHECK_STATUS_RET(*status, pb->phls[pb->nphls]);
  pb->nphls++;
   
  return(pb->phls[pb->nphls-1]);
}


float getSimputMIdpSpecBandFlux(SimputMIdpSpec* const spec,
				const float emin, 
				const float emax)
{
  // Return value.
  float flux=0.;

  long ii;
  for (ii=0; ii<spec->nentries; ii++) {
    float binmin, binmax;
    getMIdpSpecEbounds(spec, ii, &binmin, &binmax);
    if ((emin<binmax) && (emax>binmin)) {
      float min=MAX(binmin, emin);
      float max=MIN(binmax, emax);
      assert(max>min);
      flux+=(max-min)*spec->fluxdensity[ii]*spec->energy[ii];
    }
  }

  // Convert units of 'flux' from [keV/cm^2]->[erg/cm^2].
  flux*=keV2erg;

  return(flux);
}


float getSimputPhotonRate(SimputCtlg* const cat,
			  SimputSrc* const src,
			  const double prevtime, 
			  const double mjdref,
			  int* const status)
{
  // Check if the photon rate has already been determined before.
  if (NULL==src->phrate) {
    // Obtain the spectrum.
    char specref[SIMPUT_MAXSTR];
    getSimputSrcSpecRef(cat, src, prevtime, mjdref, specref, status);
    CHECK_STATUS_RET(*status, 0.);

    int spectype=getSimputExtType(cat, specref, status);
    CHECK_STATUS_RET(*status, 0.);

    // Check if the ARF is defined.
    if (NULL==cat->arf) {
      *status=EXIT_FAILURE;
      SIMPUT_ERROR("ARF not found");
      return(0.);
    }

    if (EXTTYPE_MIDPSPEC==spectype) {
      SimputMIdpSpec* midpspec=getSimputMIdpSpec(cat, specref, status);
      CHECK_STATUS_RET(*status, 0.);

      // Flux in the reference energy band.
      float refband_flux=
	getSimputMIdpSpecBandFlux(midpspec, src->e_min, src->e_max);

      SimputSpec* spec=getSimputSpec(cat, specref, status);
      CHECK_STATUS_RET(*status, 0.);
      
      // Store the determined photon rate in the source data structure
      // for later use.
      // Allocate memory.
      src->phrate=(float*)malloc(sizeof(float));
      CHECK_NULL_RET(src->phrate, *status,
		     "memory allocation for photon rate buffer failed", 0.);
      *(src->phrate)=
	src->eflux / refband_flux * 
	(float)(spec->distribution[cat->arf->NumberEnergyBins-1]);
      
    } else if (EXTTYPE_PHLIST==spectype) {

      // Get the photon list.
      SimputPhList* phl=getSimputPhList(cat, specref, status);
      CHECK_STATUS_RET(*status, 0.);
      
      // Determine the flux in the reference energy band
      // and the reference number of photons after weighing
      // with the instrument ARF.
      double refband_flux=0.; // [erg]
      double refnumber=0.; // [photons*cm^2]
      const long buffsize=10000;
      long ii;
      for (ii=0; ii*buffsize<phl->nphs; ii++) {
	// Read a block of photons.
	int anynul=0;
	long nphs=MIN(buffsize, phl->nphs-(ii*buffsize));
	float buffer[buffsize];
	fits_read_col(phl->fptr, TFLOAT, phl->cenergy, ii*buffsize+1, 
		      1, nphs, NULL, buffer, &anynul, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading unit of energy column in photon list");
	  return(0.);
	}

	// Determine the illuminated energy in the
	// reference energy band.
	long jj;
	for (jj=0; jj<nphs; jj++) {
	  float energy=buffer[jj]*phl->fenergy;
	  if ((energy>=src->e_min)&&(energy<=src->e_max)) {
	    refband_flux+=energy*keV2erg;
	  }

	  // Determine the ARF value for this energy.
	  long upper=cat->arf->NumberEnergyBins-1, lower=0, mid;
	  while (upper>lower) {
	    mid=(lower+upper)/2;
	    if (cat->arf->HighEnergy[mid]<energy) {
	      lower=mid+1;
	    } else {
	      upper=mid;
	    }
	  }
	  refnumber+=cat->arf->EffArea[lower];	  
	}
	// END of loop over all photons in the buffer.
      }

      // Store the determined photon rate in the source data structure
      // for later use.
      // Allocate memory.
      src->phrate=(float*)malloc(sizeof(float));
      CHECK_NULL_RET(src->phrate, *status,
		     "memory allocation for photon rate buffer failed", 0.);
      *(src->phrate)=src->eflux / refband_flux * refnumber;

    } else {
      SIMPUT_ERROR("could not find valid spectrum extension");
      *status=EXIT_FAILURE;
      return(0.);
    }
  }
  
  return(*(src->phrate));
}


int getSimputPhotonTime(SimputCtlg* const cat,
			SimputSrc* const src,
			double prevtime,
			const double mjdref,
			double* const nexttime,
			int* const status)
{
  // Determine the time of the next photon.

  // Determine the reference to the timing extension.
  char timeref[SIMPUT_MAXSTR];
  getSrcTimeRef(cat, src, timeref);
  CHECK_STATUS_RET(*status, 0);
    
  // Check if a timing extension has been specified.
  if (0==strlen(timeref)) {
    // The source has a constant brightness.

    // Determine the average photon rate.
    float avgrate=getSimputPhotonRate(cat, src, prevtime, mjdref, status);
    CHECK_STATUS_RET(*status, 0);
    
    // Check if the rate is 0.
    if (0.==avgrate) {
      return(1);
    }
    assert(avgrate>0.);

    // Time intervals between subsequent photons are exponentially 
    // distributed.
    *nexttime=prevtime+rndexp((double)1./avgrate, status);
    CHECK_STATUS_RET(*status, 0);

    // Successfully produced a photon.
    return(0);

  } else {
    // The source has a time-variable brightness.
    int timetype=getSimputExtType(cat, timeref, status);
    CHECK_STATUS_RET(*status, 0);

    // Check if the extension type of the timing reference.
    if (EXTTYPE_PHLIST==timetype) {
      // The timing reference points to a photon list.

      // Get the photon list.
      SimputPhList* phl=getSimputPhList(cat, timeref, status);
      CHECK_STATUS_RET(*status, 0);

      // Determine the acceptance rate for photons.
      if (0.==phl->accrate) {
	// Determine the average photon rate.
	float avgrate=getSimputPhotonRate(cat, src, prevtime, mjdref, status);
	CHECK_STATUS_RET(*status, 0);
	
	// Check if the rate is 0.
	if (0.==avgrate) {
	  return(1);
	}
	assert(avgrate>0.);

	// Acceptance rate.
	assert(phl->nphs>0);
	phl->accrate=avgrate *(phl->tstop-phl->tstart)/phl->nphs;
	assert(phl->accrate>0.);

	// Check if the rate greater than 1. In that case the rate of
	// photons in the list is too low to fulfill the requested
	// source flux.
	if (phl->accrate>1.0) {
	  *status=EXIT_FAILURE;
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "required photon rate (%e) higher than provided "
		  "by photon list (%e)", 
		  avgrate, phl->nphs/(phl->tstop-phl->tstart));
	  SIMPUT_ERROR(msg);
	  return(0);
	}
      }

      // Verify that the time column is present.
      if (0==phl->ctime) {
	*status=EXIT_FAILURE;
	SIMPUT_ERROR("photon list does not contain a time column");
	return(0);
      }

      // Select a photon.
      double rand=0.0;
      double newtime;
      do {
	// Move one row further.
	phl->currrow++;

	// Check if the end of the list has been reached.
	if (phl->currrow>phl->nphs) {
	  // No valid photon could be selected.
	  return(1);
	}

	// Read the time from the file.
	int anynul=0;
	fits_read_col(phl->fptr, TDOUBLE, phl->ctime, phl->currrow, 1, 1, 
		      NULL, &newtime, &anynul, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading time from photon list");
	  return(0);
	}
	newtime*=phl->ftime;
	
	// Check if the time lies within the requested interval.
	if (newtime+phl->timezero+(phl->mjdref-mjdref)*24.*3600.<prevtime) {
	  continue;
	}

	// Determine a random number in order to apply the acceptance rate.
	rand=getRndNum(status);
	CHECK_STATUS_RET(*status, 0);

      } while (rand>=phl->accrate);

      // Successfully produced a photon.
      *nexttime=newtime;
      return(0);

    } else if ((EXTTYPE_LC==timetype) || (EXTTYPE_PSD==timetype)) {
      // The timing reference points either to a light curve
      // or a PSD.

      // Determine the average photon rate.
      float avgrate=getSimputPhotonRate(cat, src, prevtime, mjdref, status);
      CHECK_STATUS_RET(*status, 0);
    
      // Check if the rate is 0.
      if (0.==avgrate) {
	return(1);
      }
      assert(avgrate>0.);

      // Get the light curve.
      SimputLC* lc=getSimputLC(cat, src, timeref, prevtime, mjdref, status);
      CHECK_STATUS_RET(*status, 0);
      
      // Determine the photon time according to the light curve.
      // The light curve is a piece-wise constant function, so the
      // general algorithm proposed by Klein & Roberts has to 
      // be applied.
      // Step 1 in the algorithm.
      double u=getRndNum(status);
      CHECK_STATUS_RET(*status, 0);

      // Determine the respective index kk of the light curve.
      long long nperiods=0;
      long kk=getLCBin(lc, prevtime, mjdref, &nperiods, status);
      CHECK_STATUS_RET(*status, 0);
      
      while ((kk<lc->nentries-1)||(lc->src_id>0)) {
	
	// If the end of the light curve is reached, check if it has
	// been produced from a PSD. In that case one can create new one.
	if ((kk>=lc->nentries-1)&&(lc->src_id>0)) {
	  lc=getSimputLC(cat, src, timeref, prevtime, mjdref, status);
	  CHECK_STATUS_RET(*status, 0);
	  kk=getLCBin(lc, prevtime, mjdref, &nperiods, status);
	  CHECK_STATUS_RET(*status, 0);
	}

	// Determine the relative time within the kk-th interval 
	// (i.e., t=0 lies at the beginning of the kk-th interval).
	double tk=getLCTime(lc, kk, nperiods, mjdref);
	double t =prevtime-tk;
	double stepwidth=
	  getLCTime(lc, kk+1, nperiods, mjdref)-tk;

	// Make sure that FLUXSCAL and stepwidth are positive.
	assert(lc->fluxscal>0.0);
	assert(stepwidth>0.0);
	double ak=(lc->flux[kk+1]-lc->flux[kk])/lc->fluxscal/stepwidth;
	double bk=lc->flux[kk]/lc->fluxscal;

	// Make sure that stepwidth is positive.
	if (stepwidth<=0.0) {
	  *status=EXIT_FAILURE;
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "encountered nonpositive step width (%es) in light curve '%s'",
		  stepwidth, lc->fileref);
	  SIMPUT_ERROR(msg);
	  return(0);
	}

	// Step 2 in the algorithm.
	double uk=1.-exp((-ak/2.*(pow(stepwidth,2.)-pow(t,2.))
			  -bk*(stepwidth-t))*avgrate);

	// Step 3 in the algorithm.
	if (u<=uk) {
	  if (fabs(ak*stepwidth)>fabs(bk*1.e-6)) { 
	    // Instead of checking if ak = 0., check, whether its product 
	    // with the interval length is a very small number in comparison 
	    // to b_kk. If ak * stepwidth is much smaller than b_kk, the 
	    // rate in the interval can be assumed to be approximately constant.
	    *nexttime=tk+
	      (-bk+sqrt(pow(bk,2.)+pow(ak*t,2.)+2.*ak*t*bk-2.*ak*log(1.-u)/avgrate))/ak;
	    return(0);
	    
	  } else { // ak == 0
	    *nexttime=prevtime-log(1.-u)/(avgrate*bk);
	    return(0);
	  }
	  
	} else {
	  // Step 4 (u > u_k).
	  u=(u-uk)/(1-uk);
	  kk++;
	  if ((kk>=lc->nentries-1)&&(NULL!=lc->phase)) {
	    kk=0;
	    nperiods++;
	  }
	  prevtime=getLCTime(lc, kk, nperiods, mjdref);
	}
      }
      // END of while (kk < lc->nentries).
    
      // The range of the light curve has been exceeded.
      // So the routine has failed to determine a photon time.
      return(1);

    } else {
      // The timing reference does not point to any of the 
      // above extension types.
      *status=EXIT_FAILURE;
      SIMPUT_ERROR("invalid timing extension");
      return(0);
    }
  }
}


static void getSimputPhFromPhList(const SimputCtlg* const cat,
				  SimputPhList* const phl, 
				  float* const energy, 
				  double* const ra, 
				  double* const dec,
				  int* const status)
{
  // Check if we have to read from a particular row in the FITS file,
  // or if we need to return a randomly selected photon.
  if (phl->currrow>0) {
    // Read a photon from a particular row.
    int anynul=0;
    fits_read_col(phl->fptr, TFLOAT, phl->cenergy, phl->currrow, 1, 1,
		  NULL, energy, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading energy from photon list");
      return;
    }
    *energy *=phl->fenergy;
    
    fits_read_col(phl->fptr, TDOUBLE, phl->cra, phl->currrow, 1, 1, 
		  NULL, ra, &anynul, status);      
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading right ascension from photon list");
      return;
    }
    *ra *=phl->fra;
	
    fits_read_col(phl->fptr, TDOUBLE, phl->cdec, phl->currrow, 1, 1, 
		  NULL, dec, &anynul, status);
    if (EXIT_SUCCESS!=*status) {
      SIMPUT_ERROR("failed reading declination from photon list");
      return;
    }
    *dec *=phl->fdec;

    return;

  } else {
    // Randomly select a photon from the file.
    
    // Determine the maximum value of the instrument ARF.
    if (0.==phl->refarea) {
      long kk;
      for (kk=0; kk<cat->arf->NumberEnergyBins; kk++) {
	if (cat->arf->EffArea[kk]>phl->refarea) {
	  phl->refarea=cat->arf->EffArea[kk];
	}
      }
    }
    
    while(1) {
      // Determine a random photon within the list.
      long ii=(long)(getRndNum(status)*phl->nphs);
      CHECK_STATUS_VOID(*status);
      
      // Read the photon energy.
      int anynul=0;
      fits_read_col(phl->fptr, TFLOAT, phl->cenergy, ii+1, 1, 1, 
		    NULL, energy, &anynul, status);
      if (EXIT_SUCCESS!=*status) {
	SIMPUT_ERROR("failed reading energy from photon list");
	return;
      }
      *energy *=phl->fenergy;
      
      // Determine the ARF value for the photon energy.
      long upper=cat->arf->NumberEnergyBins-1, lower=0, mid;
      while (upper>lower) {
	mid=(lower+upper)/2;
	if (cat->arf->HighEnergy[mid]<*energy) {
	  lower=mid+1;
	} else {
	  upper=mid;
	}
      }
      
      // Randomly determine according to the effective area
      // of the instrument, whether this photon is seen or not.
      double r=getRndNum(status);
      CHECK_STATUS_VOID(*status);
      if (r<cat->arf->EffArea[lower]/phl->refarea) {
	// Read the position of the photon.
	fits_read_col(phl->fptr, TDOUBLE, phl->cra, ii+1, 1, 1, 
		      NULL, ra, &anynul, status);      
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading right ascension from photon list");
	  return;
	}
	*ra *=phl->fra;
	
	fits_read_col(phl->fptr, TDOUBLE, phl->cdec, ii+1, 1, 1, 
		      NULL, dec, &anynul, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading declination from photon list");
	  return;
	}
	*dec *=phl->fdec;
     
	// Increase the counter of the number of returned photons and
	// check if it exceeds one fifth of the total number of available
	// photons.
	if (0==(++phl->nrphs) % (phl->nphs/5)) {
	  char msg[SIMPUT_MAXSTR];
	  float ratio=phl->nrphs*1./phl->nphs;
	  sprintf(msg, "ratio of the number of randomly drawn photons (%ld) "
		  "versus total number of photons (%ld) exceeds %.0lf%%! ",
		  phl->nrphs, phl->nphs, ratio*100.);
	  if (ratio<1.) {
	    strcat(msg, "Individual photons might be used multiple times");
	  } else {
	    strcat(msg, "Individual photons are used multiple times");
	  }
	  SIMPUT_WARNING(msg);
	}
      
	return;
      }
    }
  }
}


void getSimputPhotonEnergyCoord(SimputCtlg* const cat,
				SimputSrc* const src,
				double currtime,
				const double mjdref,
				float* const energy,
				double* const ra,
				double* const dec,
				int* const status)
{
  // Determine the references to the spectrum and 
  // image for the updated photon arrival time.
  char specref[SIMPUT_MAXSTR];
  getSimputSrcSpecRef(cat, src, currtime, mjdref, specref, status);
  CHECK_STATUS_VOID(*status);
  char imagref[SIMPUT_MAXSTR];
  getSrcImagRef(cat, src, currtime, mjdref, imagref, status);
  CHECK_STATUS_VOID(*status);

  
  // Determine the extension type of the spectrum and the image
  // reference.
  int spectype=getSimputExtType(cat, specref, status);
  CHECK_STATUS_VOID(*status);
  int imagtype=getSimputExtType(cat, imagref, status);
  CHECK_STATUS_VOID(*status);


  // If the spectrum or the image reference point to a photon list,
  // determine simultaneously the energy and spatial information.
  SimputPhList* phl=NULL;
  if (EXTTYPE_PHLIST==spectype) {
    phl=getSimputPhList(cat, specref, status);
    CHECK_STATUS_VOID(*status);
  } else if (EXTTYPE_PHLIST==imagtype) {
    phl=getSimputPhList(cat, imagref, status);
    CHECK_STATUS_VOID(*status);
  }
  if (NULL!=phl) {
    float b_energy;
    double b_ra, b_dec;
    getSimputPhFromPhList(cat, phl, &b_energy, &b_ra, &b_dec, status);
    CHECK_STATUS_VOID(*status);

    if (EXTTYPE_PHLIST==spectype) {
      *energy=b_energy;
    }
    if (EXTTYPE_PHLIST==imagtype) {
      // Shift the photon position according to the 
      // RA,Dec values defined for this source in the catalog.

      // Apply IMGSCAL.
      b_ra *=1./src->imgscal*cos(b_dec)/cos(b_dec/src->imgscal);
      b_dec*=1./src->imgscal;

      // Get a Carteesian coordinate vector for the photon location.
      Vector p=unit_vector(b_ra, b_dec);

      // Apply IMGROTA by rotation around the x-axis.
      double cosimgrota=cos(src->imgrota);
      double sinimgrota=sin(src->imgrota);
      Vector r;
      r.x= p.x;
      r.y= cosimgrota*p.y + sinimgrota*p.z;
      r.z=-sinimgrota*p.y + cosimgrota*p.z;

      // Rotate the vector towards the source position.
      double cosra=cos(src->ra);
      double sinra=sin(src->ra);
      double cosdec=cos(src->dec);
      double sindec=sin(src->dec);
      Vector f;
      f.x=r.x*cosra*cosdec - r.y*sinra - r.z*cosra*sindec;
      f.y=r.x*sinra*cosdec + r.y*cosra - r.z*sinra*sindec;
      f.z=r.x      *sindec +     0.0   + r.z      *cosdec;

      // Determine RA and Dec of the photon.
      calculate_ra_dec(f, ra, dec); 
    }
  }


  // ---
  // If the spectrum does NOT refer to a photon list, determine the 
  // photon energy from a spectrum.
  if (EXTTYPE_MIDPSPEC==spectype) {
    // Determine the spectrum.
    SimputSpec* spec=getSimputSpec(cat, specref, status);
    CHECK_STATUS_VOID(*status);

    // Get a random number in the interval [0,1].
    double rnd=getRndNum(status);
    CHECK_STATUS_VOID(*status);
    assert(rnd>=0.);
    assert(rnd<=1.);

    // Multiply with the total photon rate (i.e. the spectrum
    // does not have to be normalized).
    rnd*=spec->distribution[cat->arf->NumberEnergyBins-1];

    // Determine the corresponding point in the spectral 
    // distribution (using binary search).
    long upper=cat->arf->NumberEnergyBins-1, lower=0, mid;
    while (upper>lower) {
      mid=(lower+upper)/2;
      if (spec->distribution[mid]<rnd) {
	lower=mid+1;
      } else {
	upper=mid;
      }
    }

    // Return the corresponding photon energy.
    *energy=
      cat->arf->LowEnergy[lower] + 
      getRndNum(status)*
      (cat->arf->HighEnergy[lower]-cat->arf->LowEnergy[lower]);
    CHECK_STATUS_VOID(*status);
  }
  // END of determine the photon energy.
  // --- 


  // ---
  // If the image does NOT refer to a photon list, determine the
  // spatial information.

  // Point-like sources.
  if (EXTTYPE_NONE==imagtype) {
    *ra =src->ra;
    *dec=src->dec;
  }

  // Spatially extended sources.
  else if (EXTTYPE_IMAGE==imagtype) {
    // Determine the photon direction from an image.
    struct wcsprm wcs={ .flag=-1 };

    do { // Error handling loop.

      // Determine the image.
      SimputImg* img=getSimputImg(cat, imagref, status);
      CHECK_STATUS_BREAK(*status);

      // Perform a binary search in 2 dimensions.
      double rnd=getRndNum(status)*img->dist[img->naxis1-1][img->naxis2-1];
      CHECK_STATUS_BREAK(*status);

      // Perform a binary search to obtain the x-coordinate.
      long high=img->naxis1-1;
      long xl  =0;
      long mid;
      long ymax=img->naxis2-1;
      while (high > xl) {
	mid=(xl+high)/2;
	if (img->dist[mid][ymax] < rnd) {
	  xl=mid+1;
	} else {
	  high=mid;
	}
      }

      // Search for the y coordinate.
      high=img->naxis2-1;
      long yl=0;
      while (high > yl) {
	mid=(yl+high)/2;
	if (img->dist[xl][mid] < rnd) {
	  yl=mid+1;
	} else {
	  high=mid;
	}
      }
      // Now xl and yl have pixel positions [long pixel coordinates].

      // Create a temporary wcsprm data structure, which can be modified
      // to fit this particular source. The wcsprm data structure contained 
      // in the image should not be modified, since it is used for all 
      // sources including the image.
      wcscopy(1, img->wcs, &wcs);

      // Set the position to the origin and assign the correct scaling.
      // TODO: This assumes that the image WCS is equivalent to the 
      // coordinate system used in the catalog!!
      wcs.crval[0] =src->ra *180./M_PI;
      wcs.crval[1] =src->dec*180./M_PI;
      wcs.cdelt[0]*=1./src->imgscal;
      wcs.cdelt[1]*=1./src->imgscal;
      wcs.flag=0;

      // Check that CUNIT is set to "deg". Otherwise there will be a conflict
      // between CRVAL [deg] and CDELT [different unit]. 
      // TODO This is not required by the standard.
      if (((0!=strcmp(wcs.cunit[0], "deg     ")) &&
    		  (0!=strcmp(wcs.cunit[0], "degree  ")) &&
    		  (0!=strcmp(wcs.cunit[0], "deg")) &&
    		  (0!=strcmp(wcs.cunit[0], "degree"))) ||
    		  ((0!=strcmp(wcs.cunit[1], "deg     ")) &&
    				  (0!=strcmp(wcs.cunit[1], "degree  ")) &&
    				  (0!=strcmp(wcs.cunit[1], "deg")) &&
    				  (0!=strcmp(wcs.cunit[1], "degree")))) {
    	  *status=EXIT_FAILURE;
    	  char msg[SIMPUT_MAXSTR];
    	  sprintf(msg, "units of image coordinates are '%s' and '%s' "
    			  "(must be 'deg')", wcs.cunit[0], wcs.cunit[1]);
    	  SIMPUT_ERROR(msg);
    	  break;
      }

      // Determine floating point pixel positions shifted by 0.5 in 
      // order to match the FITS conventions and with a randomization
      // over the pixels.
      double xd=(double)xl + 0.5 + getRndNum(status);
      CHECK_STATUS_BREAK(*status);
      double yd=(double)yl + 0.5 + getRndNum(status);
      CHECK_STATUS_BREAK(*status);

      // Rotate the image (pixel coordinates) by IMGROTA around the 
      // reference point.
      double xdrot=
	(xd-wcs.crpix[0])*cos(src->imgrota) +
	(yd-wcs.crpix[1])*sin(src->imgrota) + wcs.crpix[0];
      double ydrot=
	-(xd-wcs.crpix[0])*sin(src->imgrota) +
	 (yd-wcs.crpix[1])*cos(src->imgrota) + wcs.crpix[1];
      
      // Convert the long-valued pixel coordinates to double values,
      // including a randomization over the pixel and transform from 
      // pixel coordinates to RA and DEC ([rad]) using the  WCS information.
      p2s(&wcs, xdrot, ydrot, ra, dec, status);
      CHECK_STATUS_BREAK(*status);

      // Determine the RA in the interval from [0:2pi).
      while(*ra>=2.*M_PI) {
	*ra-=2.*M_PI;
      }
      while(*ra<0.) {
	*ra+=2.*M_PI;
      }
      
    } while(0); // END of error handling loop.

    // Release memory.
    wcsfree(&wcs);
  }
  // END of determine the photon direction.
  // ---

  return;
}


int getSimputPhoton(SimputCtlg* const cat,
		    SimputSrc* const src,
		    double prevtime,
		    const double mjdref,
		    double* const nexttime,
		    float* const energy,
		    double* const ra,
		    double* const dec,
		    int* const status)
{
  // Determine the time of the next photon.
  int failed=getSimputPhotonTime(cat, src, prevtime, mjdref, nexttime, status);
  CHECK_STATUS_RET(*status, failed);
  if (failed>0) return(failed);

  // Determine the energy and the direction of origin of the photon.
  getSimputPhotonEnergyCoord(cat, src, *nexttime, mjdref,
			     energy, ra, dec, status);
  CHECK_STATUS_RET(*status, 0);

  return(0);
}

int precompute_photon (SimputCtlg *cat, long sourcenumber,
		       double mjdref, double prevtime,
		       SimputPhoton *next_photons,
		       int* const status)
{
  int lightcurve_status;
  double time;
  float energy;
  double ra;
  double dec;
  SimputSrc* src = NULL;

  src = getSimputSrc(cat, sourcenumber, status);
  CHECK_STATUS_RET(*status, 0);

  lightcurve_status = getSimputPhoton(cat, src,
				      prevtime, mjdref,
				      &time, &energy, &ra, &dec,
				      status);
  CHECK_STATUS_RET(*status, 0);

  // sourcenumber starts at index 1, but C array start at index 0
  next_photons[sourcenumber-1].time = time;
  next_photons[sourcenumber-1].energy = energy;
  next_photons[sourcenumber-1].ra = ra;
  next_photons[sourcenumber-1].dec = dec;
  next_photons[sourcenumber-1].lightcurve_status = lightcurve_status;
  next_photons[sourcenumber-1].polarization = 0;
  return 0;
}


SimputPhoton* startSimputPhotonAnySource(SimputCtlg* const cat,
				const double mjdref,
				int* const status)
{
  long ii, n_sources;
  n_sources = getSimputCtlgNSources(cat);

  SimputPhoton* next_photons = malloc(n_sources * sizeof(*next_photons));
  CHECK_NULL_RET(next_photons, *status,
		 "memory allocation for photon cache failed", NULL);
  for (ii=1; ii < (n_sources+1); ii++) {
    precompute_photon(cat, ii, mjdref, 0., next_photons, status);
    CHECK_STATUS_RET(*status, NULL);
  }
  return next_photons;
}


int getSimputPhotonAnySource(SimputCtlg* const cat,
			     SimputPhoton *next_photons,
			     const double mjdref,
			     double* const time,
			     float* const energy,
			     double* const ra,
			     double* const dec,
			     double* const polarization,
			     long* const source_index,
			     int* const status)
{
  long ii, n_sources, next_index, number_valid_photons;
  double min_time = DBL_MAX;
  n_sources = getSimputCtlgNSources(cat);

  CHECK_NULL_RET(next_photons, *status,
		   "next_photons has not been initialized", 1);

  // From the list of pre-generated photons, find the next one
  for (ii=0; ii < n_sources; ii++) {
    if ((next_photons[ii].time < min_time) && (next_photons[ii].lightcurve_status == 0)){
      min_time = next_photons[ii].time;
      next_index = ii;
      number_valid_photons++;
    }
  }
  if (number_valid_photons == 0){
    return 1;
  }

  *time = next_photons[next_index].time;
  *energy = next_photons[next_index].energy;
  *ra = next_photons[next_index].ra;
  *dec = next_photons[next_index].dec;
  *polarization = next_photons[next_index].polarization;
  *source_index = next_index + 1; //FTIS file index starts at 1

  // Replace the photon in the precomputed list with the next one.
  precompute_photon(cat, *source_index, mjdref, *time, next_photons, status);

  return 0;
}

void closeSimputPhotonAnySource(SimputPhoton *next_photons)
{
  if (next_photons!=NULL){
    free(next_photons);
    next_photons = NULL;
  }
}

float getSimputSrcExt(SimputCtlg* const cat,
		      const SimputSrc* const src,
		      const double prevtime,
		      const double mjdref,
		      int* const status)
{
  // Return value [rad].
  float extension=0.;

  struct wcsprm wcs={ .flag=-1 };

  do { // Error handling loop.

    // Get the source image for this particular source.
    char imagref[SIMPUT_MAXSTR];
    getSrcImagRef(cat, src, prevtime, mjdref, imagref, status);
    CHECK_STATUS_BREAK(*status);
    int imagtype=getSimputExtType(cat, imagref, status);
    CHECK_STATUS_RET(*status, 0);

    // Check if it is a point-like or an extended source.
    if (EXTTYPE_NONE==imagtype) {
      // Point-like source.
      extension=0.;
      break;
      
    } else if (EXTTYPE_IMAGE==imagtype) {
      // Extended source => determine the maximum extension.
      double maxext=0.;

      SimputImg* img=getSimputImg(cat, imagref, status);
      CHECK_STATUS_BREAK(*status);

      // Copy the wcsprm structure and change the size 
      // according to IMGSCAL.
      wcscopy(1, img->wcs, &wcs);

      // Change the scale of the image according to the source specific
      // IMGSCAL property.
      wcs.crval[0] = 0.; 
      wcs.crval[1] = 0.;
      wcs.cdelt[0]*= 1./src->imgscal;
      wcs.cdelt[1]*= 1./src->imgscal;
      wcs.flag = 0;

      // Check lower left corner.
      double px=0.5;
      double py=0.5;
      double sx, sy;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      // TODO This has not been tested extensively, yet.
      while(sx>M_PI) {
	sx-=2.0*M_PI;
      }
      double ext=sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext=ext;
      }

      // Check lower right corner.
      px = img->naxis1*1. + 0.5;
      py = 0.5;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      while(sx>M_PI) {
	sx-=2.0*M_PI;
      }
      ext = sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext = ext;
      }
      
      // Check upper left corner.
      px=0.5;
      py=img->naxis2*1. + 0.5;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      while(sx>M_PI) {
	sx-=2.0*M_PI;
      }
      ext=sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext=ext;
      }
      
      // Check upper right corner.
      px=img->naxis1*1. + 0.5;
      py=img->naxis2*1. + 0.5;
      p2s(&wcs, px, py, &sx, &sy, status);
      CHECK_STATUS_BREAK(*status);
      while(sx>M_PI) {
	sx-=2.0*M_PI;
      }
      ext=sqrt(pow(sx,2.)+pow(sy,2.));
      if (ext>maxext) {
	maxext=ext;
      }
      
      extension=(float)maxext;

    } else if (EXTTYPE_PHLIST==imagtype) {
      // Get the photon list.
      SimputPhList* phl=getSimputPhList(cat, imagref, status);
      CHECK_STATUS_BREAK(*status);
      
      // Determine the maximum extension by going through
      // the whole list of photons.
      double maxext=0.;
      const long buffsize=10000;
      double rabuffer[buffsize];
      double decbuffer[buffsize];
      long ii;
      for (ii=0; ii*buffsize<phl->nphs; ii++) {
	// Read a block of photons.
	int anynul=0;
	long nphs=MIN(buffsize, phl->nphs-(ii*buffsize));
	fits_read_col(phl->fptr, TDOUBLE, phl->cra, ii*buffsize+1, 
		      1, nphs, NULL, rabuffer, &anynul, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading right ascension from photon list");
	  return(0.);
	}
	fits_read_col(phl->fptr, TDOUBLE, phl->cdec, ii*buffsize+1, 
		      1, nphs, NULL, decbuffer, &anynul, status);
	if (EXIT_SUCCESS!=*status) {
	  SIMPUT_ERROR("failed reading declination from photon list");
	  return(0.);
	}

	// Determine the maximum extension.
	long jj;
	for (jj=0; jj<nphs; jj++) {
	  double ext=
	    sqrt(pow(rabuffer[jj],2.0)+pow(decbuffer[jj],2.0))*M_PI/180.;
	  if (ext>maxext) {
	    maxext=ext;
	  }
	}
      }

      extension=(float)maxext;
    }
  } while(0); // END of error handling loop.

  // Release memory.
  wcsfree(&wcs);

  return(extension);
}

