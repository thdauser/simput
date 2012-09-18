#include "common.h"


SimputSrc* getSimputSrc(int* const status)
{
  SimputSrc* entry=(SimputSrc*)malloc(sizeof(SimputSrc));
  CHECK_NULL_RET(entry, *status, 
		 "memory allocation for SimputSrc failed", entry);

  // Initialize elements.
  entry->src_id  =0;
  entry->src_name=NULL;
  entry->ra      =0.;
  entry->dec     =0.;
  entry->imgrota =0.;
  entry->imgscal =1.;
  entry->e_min   =0.;
  entry->e_max   =0.;
  entry->eflux   =0.;
  entry->spectrum=NULL;
  entry->image   =NULL;
  entry->timing  =NULL;

  return(entry);
}


SimputSrc* getSimputSrcV(const long src_id, 
			 const char* const src_name,
			 const double ra,
			 const double dec,
			 const float imgrota,
			 const float imgscal,
			 const float e_min,
			 const float e_max,
			 const float eflux,
			 const char* const spectrum,
			 const char* const image,
			 const char* const timing,
			 int* const status)
{
  SimputSrc* src=getSimputSrc(status);
  CHECK_STATUS_RET(*status, src);

  // Initialize with the given values.
  if (src_id<0) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "SRC_ID (%ld) must have a positive value", src_id);
    SIMPUT_ERROR(msg);
    *status=EXIT_FAILURE;
    return(src);
  }
  src->src_id = src_id;

  src->src_name=(char*)malloc((strlen(src_name)+1)*sizeof(char));
  CHECK_NULL_RET(src->src_name, *status,
		 "memory allocation for source name failed", src);
  strcpy(src->src_name, src_name);

  src->ra      = ra;
  src->dec     = dec;
  src->imgrota = imgrota;
  src->imgscal = imgscal;
  src->e_min   = e_min;
  src->e_max   = e_max;
  src->eflux   = eflux;

  src->spectrum=(char*)malloc((strlen(spectrum)+1)*sizeof(char));
  CHECK_NULL_RET(src->spectrum, *status,
		 "memory allocation for reference to spectrum failed", src);
  strcpy(src->spectrum, spectrum);

  if (NULL!=image) {
    src->image=(char*)malloc((strlen(image)+1)*sizeof(char));
    CHECK_NULL_RET(src->image, *status,
		   "memory allocation for reference to image extension failed", 
		   src);
    strcpy(src->image, image);
  }

  if (NULL!=timing) {
    src->timing=(char*)malloc((strlen(timing)+1)*sizeof(char));
    CHECK_NULL_RET(src->timing, *status,
		   "memory allocation for reference to timing extension failed", 
		   src);
    strcpy(src->timing, timing);
  }

  return(src);
}


void freeSimputSrc(SimputSrc** const src)
{
  if (NULL!=*src) {
    if (NULL!=(*src)->src_name) {
      free((*src)->src_name);
    }
    if (NULL!=(*src)->spectrum) {
      free((*src)->spectrum);
    }
    if (NULL!=(*src)->image) {
      free((*src)->image);
    }
    if (NULL!=(*src)->timing) {
      free((*src)->timing);
    }
    free(*src);
    *src=NULL;
  }
}


SimputCtlg* getSimputCtlg(int* const status)
{
  SimputCtlg* cat=
    (SimputCtlg*)malloc(sizeof(SimputCtlg));
  CHECK_NULL_RET(cat, *status, 
		 "memory allocation for SimputCtlg failed", cat);

  // Initialize elements.
  cat->fptr     =NULL;
  cat->nentries =0;
  cat->csrc_id  =0;
  cat->csrc_name=0;
  cat->cra      =0;
  cat->cdec     =0;
  cat->cimgrota =0;
  cat->cimgscal =0;
  cat->ce_min   =0;
  cat->ce_max   =0;
  cat->cflux    =0;
  cat->cspectrum=0;
  cat->cimage   =0;
  cat->ctiming  =0;
  cat->fra      =0.;
  cat->fdec     =0.;
  cat->fimgrota =0.;
  cat->fflux    =0.;
  cat->fe_min   =0.;
  cat->fe_max   =0.;
  cat->filepath =NULL;
  cat->filename =NULL;
  cat->srcbuff  =NULL;
  cat->midpspecbuff=NULL;
  cat->phlistbuff  =NULL;
  cat->lcbuff   =NULL;
  cat->psdbuff  =NULL;
  cat->imgbuff  =NULL;
  cat->specbuff =NULL;
  cat->krlcbuff =NULL;
  cat->arf      =NULL;

  return(cat);
}


void freeSimputCtlg(SimputCtlg** const cat,
		    int* const status)
{
  if (NULL!=*cat) {
    if (NULL!=(*cat)->fptr) {
      fits_close_file((*cat)->fptr, status);
    }
    if (NULL!=(*cat)->filepath) {
      free((*cat)->filepath);
    }
    if (NULL!=(*cat)->filename) {
      free((*cat)->filename);
    }
    if (NULL!=(*cat)->srcbuff) {
      freeSimputSourceBuffer((struct SimputSourceBuffer**)&((*cat)->srcbuff));
    }
    // TODO midpspecbuff
    // TODO phlistbuff
    if (NULL!=(*cat)->lcbuff) {
      freeSimputLCBuffer((struct SimputLCBuffer**)&((*cat)->lcbuff));
    }
    // TODO psdbuff
    if (NULL!=(*cat)->imgbuff) {
      freeSimputImgBuffer((struct SimputImgBuffer**)&((*cat)->imgbuff));
    }
    if (NULL!=(*cat)->specbuff) {
      freeSimputSpectrumBuffer((struct SimputSpectrumBuffer**)&((*cat)->specbuff));
    }
    // TODO krlcbuff
    free(*cat);
    *cat=NULL;
  }
}


SimputMIdpSpec* getSimputMIdpSpec(int* const status)
{
  SimputMIdpSpec* spec=
    (SimputMIdpSpec*)malloc(sizeof(SimputMIdpSpec));
  CHECK_NULL_RET(spec, *status, 
		 "memory allocation for SimputMIdpSpec failed", spec);

  // Initialize elements.
  spec->nentries=0;
  spec->energy  =NULL;
  spec->pflux   =NULL;
  spec->name    =NULL;
  spec->fileref =NULL;

  return(spec);  
}


void freeSimputMIdpSpec(SimputMIdpSpec** const spec)
{
  if (NULL!=*spec) {
    if (NULL!=(*spec)->energy) {
      free((*spec)->energy);
    }
    if (NULL!=(*spec)->pflux) {
      free((*spec)->pflux);
    }
    if (NULL!=(*spec)->name) {
      free((*spec)->name);
    }
    if (NULL!=(*spec)->fileref) {
      free((*spec)->fileref);
    }
    free(*spec);
    *spec=NULL;
  }
}


SimputSpec* getSimputSpec(int* const status)
{
  SimputSpec* spec=(SimputSpec*)malloc(sizeof(SimputSpec));
  CHECK_NULL_RET(spec, *status, 
		 "memory allocation for SimputSpec failed", spec);

  // Initialize elements.
  spec->distribution=NULL;
  spec->fileref     =NULL;

  return(spec);  
}


void freeSimputSpec(SimputSpec** const spec)
{
  if (NULL!=*spec) {
    if (NULL!=(*spec)->distribution) {
      free((*spec)->distribution);
    }
    if (NULL!=(*spec)->fileref) {
      free((*spec)->fileref);
    }
    free(*spec);
    *spec=NULL;
  }
}


SimputLC* getSimputLC(int* const status)
{
  SimputLC* lc=(SimputLC*)malloc(sizeof(SimputLC));
  CHECK_NULL_RET(lc, *status, 
		 "memory allocation for SimputLC failed", lc);

  // Initialize elements.
  lc->nentries=0;
  lc->time    =NULL;
  lc->phase   =NULL;
  lc->flux    =NULL;
  lc->spectrum=NULL;
  lc->image   =NULL;
  lc->mjdref  =0.;
  lc->timezero=0.;
  lc->phase0  =0.;
  lc->period  =0.;
  lc->fluxscal=0.;
  lc->src_id  =0;
  lc->fileref =NULL;

  return(lc);
}


void freeSimputLC(SimputLC** const lc)
{
  if (NULL!=*lc) {
    if ((*lc)->nentries>0) {
      if (NULL!=(*lc)->spectrum) {
	long ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->spectrum[ii]) {
	    free((*lc)->spectrum[ii]);
	  }
	}
	free((*lc)->spectrum);    
      }
      if (NULL!=(*lc)->image) {
	long ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->image[ii]) {
	    free((*lc)->image[ii]);
	  }
	}
	free((*lc)->image);
      }
    }
    if (NULL!=(*lc)->time) {
      free((*lc)->time);
    }
    if (NULL!=(*lc)->phase) {
      free((*lc)->phase);
    }
    if (NULL!=(*lc)->flux) {
      free((*lc)->flux);
    }
    if (NULL!=(*lc)->fileref) {
      free((*lc)->fileref);
    }
    free(*lc);
    *lc=NULL;
  }
}


SimputKRLC* getSimputKRLC(int* const status)
{
  SimputKRLC* lc=(SimputKRLC*)malloc(sizeof(SimputKRLC));
  CHECK_NULL_RET(lc, *status, 
		 "memory allocation for SimputKRLC failed", lc);

  // Initialize elements.
  lc->nentries=0;
  lc->time    =NULL;
  lc->phase   =NULL;
  lc->a       =NULL;
  lc->b       =NULL;
  lc->spectrum=NULL;
  lc->image   =NULL;
  lc->mjdref  =0.;
  lc->timezero=0.;
  lc->phase0  =0.;
  lc->period  =0.;
  lc->fluxscal=0.;
  lc->src_id  =0;
  lc->fileref =NULL;

  return(lc);
}


void freeSimputKRLC(SimputKRLC** const lc)
{
  if (NULL!=*lc) {
    if ((*lc)->nentries>0) {
      if (NULL!=(*lc)->spectrum) {
	long ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->spectrum[ii]) {
	    free((*lc)->spectrum[ii]);
	  }
	}
	free((*lc)->spectrum);    
      }
      if (NULL!=(*lc)->image) {
	long ii;
	for (ii=0; ii<(*lc)->nentries; ii++) {
	  if (NULL!=(*lc)->image[ii]) {
	    free((*lc)->image[ii]);
	  }
	}
	free((*lc)->image);
      }
    }
    if (NULL!=(*lc)->time) {
      free((*lc)->time);
    }
    if (NULL!=(*lc)->phase) {
      free((*lc)->phase);
    }
    if (NULL!=(*lc)->a) {
      free((*lc)->a);
    }
    if (NULL!=(*lc)->b) {
      free((*lc)->b);
    }
    if (NULL!=(*lc)->fileref) {
      free((*lc)->fileref);
    }
    free(*lc);
    *lc=NULL;
  }
}


SimputPSD* getSimputPSD(int* const status)
{
  SimputPSD* psd=(SimputPSD*)malloc(sizeof(SimputPSD));
  CHECK_NULL_RET(psd, *status, 
		 "memory allocation for SimputPSD failed", psd);

  // Initialize elements.
  psd->nentries =0;
  psd->frequency=NULL;
  psd->power    =NULL;
  psd->fileref  =NULL;

  return(psd);
}


void freeSimputPSD(SimputPSD** const psd)
{
  if (NULL!=*psd) {
    if (NULL!=(*psd)->frequency) {
      free((*psd)->frequency);
    }
    if (NULL!=(*psd)->power) {
      free((*psd)->power);
    }
    if (NULL!=(*psd)->fileref) {
      free((*psd)->fileref);
    }
    free(*psd);
    *psd=NULL;
  }
}


SimputImg* getSimputImg(int* const status)
{
  SimputImg* img=(SimputImg*)malloc(sizeof(SimputImg));
  CHECK_NULL_RET(img, *status, 
		 "memory allocation for SimputImg failed", img);

  // Initialize elements.
  img->naxis1  =0;
  img->naxis2  =0;
  img->dist    =NULL;
  img->fluxscal=0.;
  img->fileref =NULL;
  img->wcs     =NULL;

  return(img);
}


void freeSimputImg(SimputImg** const img)
{
  if (NULL!=*img) {
    if (NULL!=(*img)->dist) {
      if ((*img)->naxis1>0) {
	long ii;
	for (ii=0; ii<(*img)->naxis1; ii++) {
	  if (NULL!=(*img)->dist[ii]) {
	    free((*img)->dist[ii]);
	  }
	}
      }
      free((*img)->dist);
    }
    if (NULL!=(*img)->fileref) {
      free((*img)->fileref);
    }
    if (NULL!=(*img)->wcs) {
      wcsfree((*img)->wcs);
    }
    free(*img);
    *img=NULL;
  }
}


