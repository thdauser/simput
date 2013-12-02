#include "common.h"


SimputSrc* newSimputSrc(int* const status)
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
  entry->phrate  =NULL;
  entry->spectrum=NULL;
  entry->image   =NULL;
  entry->timing  =NULL;

  return(entry);
}


SimputSrc* newSimputSrcV(const long src_id, 
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
  SimputSrc* src=newSimputSrc(status);
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
    if (NULL!=(*src)->phrate) {
      free((*src)->phrate);
    }
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


struct SimputSrcBuffer* newSimputSrcBuffer(int* const status)
{
  struct SimputSrcBuffer *srcbuff = 
    (struct SimputSrcBuffer*)malloc(sizeof(struct SimputSrcBuffer));

  CHECK_NULL_RET(srcbuff, *status, 
		 "memory allocation for SimputSrcBuffer failed", srcbuff);
    
  srcbuff->nsrcs  =0;
  srcbuff->csrc   =0;
  srcbuff->srcs   =NULL;
  srcbuff->rownums=NULL;
  srcbuff->rowmap =NULL;

  return(srcbuff);
}


void freeSimputSrcBuffer(struct SimputSrcBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->srcs) {
      long ii;
      for (ii=0; ii<(*sb)->nsrcs; ii++) {
	freeSimputSrc(&((*sb)->srcs[ii]));
      }
      free((*sb)->srcs);
    }
    if (NULL!=(*sb)->rownums) {
      free((*sb)->rownums);
    }
    if (NULL!=(*sb)->rowmap) {
      free((*sb)->rowmap);
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputCtlg* newSimputCtlg(int* const status)
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
  cat->extbuff  =NULL;
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
      freeSimputSrcBuffer((struct SimputSrcBuffer**)&((*cat)->srcbuff));
    }
    if (NULL!=(*cat)->midpspecbuff) {
      freeSimputMIdpSpecBuffer((struct SimputMIdpSpecBuffer**)&((*cat)->midpspecbuff));
    }
    if (NULL!=(*cat)->phlistbuff) {
      freeSimputPhListBuffer((struct SimputPhListBuffer**)&((*cat)->phlistbuff), status);
    }
    if (NULL!=(*cat)->lcbuff) {
      freeSimputLCBuffer((struct SimputLCBuffer**)&((*cat)->lcbuff));
    }
    if (NULL!=(*cat)->psdbuff) {
      freeSimputPSDBuffer((struct SimputPSDBuffer**)&((*cat)->psdbuff));
    }
    if (NULL!=(*cat)->imgbuff) {
      freeSimputImgBuffer((struct SimputImgBuffer**)&((*cat)->imgbuff));
    }
    if (NULL!=(*cat)->specbuff) {
      freeSimputSpecBuffer((struct SimputSpecBuffer**)&((*cat)->specbuff));
    }
    if (NULL!=(*cat)->extbuff) {
      freeSimputExttypeBuffer((struct SimputExttypeBuffer**)&((*cat)->extbuff));
    }
    free(*cat);
    *cat=NULL;
  }
}


long getSimputCtlgNSources(const SimputCtlg* const cat)
{
  return(cat->nentries);
}


struct SimputExttypeBuffer* newSimputExttypeBuffer(int* const status)
{
  struct SimputExttypeBuffer *extbuff= 
    (struct SimputExttypeBuffer*)malloc(sizeof(struct SimputExttypeBuffer));

  CHECK_NULL_RET(extbuff, *status, 
		 "memory allocation for SimputExttypeBuffer failed", extbuff);
    
  extbuff->type   =EXTTYPE_NONE;
  extbuff->fileref=NULL;
  extbuff->left   =NULL;
  extbuff->right  =NULL;

  return(extbuff);
}


void freeSimputExttypeBuffer(struct SimputExttypeBuffer** eb)
{
  if (NULL!=*eb) {
    if (NULL!=(*eb)->fileref) {
      free((*eb)->fileref);
    }
    if (NULL!=(*eb)->left) {
      freeSimputExttypeBuffer(&((*eb)->left));
    }
    if (NULL!=(*eb)->right) {
      freeSimputExttypeBuffer(&((*eb)->right));
    }
    free(*eb);
    *eb=NULL;
  }
}


int searchSimputExttypeBuffer(void* buffer, const char* const filename)
{
  struct SimputExttypeBuffer* eb=(struct SimputExttypeBuffer*)buffer;

  // Check if this is a leaf.
  if (NULL==eb) {
    return(EXTTYPE_NONE);
  }
  
  // Compare with the current type.
  int cmp=strcmp(filename, eb->fileref);
  if (0==cmp) {
    return(eb->type);
  } else if (cmp<0) {
    return(searchSimputExttypeBuffer(eb->left, filename));
  } else {
    return(searchSimputExttypeBuffer(eb->right, filename));
  }
}


void insertSimputExttypeBuffer(void** buffer,
			       const char* const filename,
			       const int type,
			       int* const status)
{
  struct SimputExttypeBuffer** eb=(struct SimputExttypeBuffer**)buffer;

  // Check if this is a leaf.
  if (NULL==*eb) {
    // Insert the extension type here.
    *eb=newSimputExttypeBuffer(status);
    CHECK_STATUS_VOID(*status);
    (*eb)->fileref=(char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_VOID((*eb)->fileref, *status, 
		    "memory allocation for file reference failed");
    strcpy((*eb)->fileref, filename);
    (*eb)->type=type;
  } else if (strcmp(filename, (*eb)->fileref) < 0){
    insertSimputExttypeBuffer((void**)&((*eb)->left), filename, type, status);
    CHECK_STATUS_VOID(*status);
  } else {
    insertSimputExttypeBuffer((void**)&((*eb)->right), filename, type, status);
    CHECK_STATUS_VOID(*status);
  }
}


SimputMIdpSpec* newSimputMIdpSpec(int* const status)
{
  SimputMIdpSpec* spec=
    (SimputMIdpSpec*)malloc(sizeof(SimputMIdpSpec));
  CHECK_NULL_RET(spec, *status, 
		 "memory allocation for SimputMIdpSpec failed", spec);

  // Initialize elements.
  spec->nentries   =0;
  spec->energy     =NULL;
  spec->fluxdensity=NULL;
  spec->name       =NULL;
  spec->fileref    =NULL;

  return(spec);  
}


void freeSimputMIdpSpec(SimputMIdpSpec** const spec)
{
  if (NULL!=*spec) {
    if (NULL!=(*spec)->energy) {
      free((*spec)->energy);
    }
    if (NULL!=(*spec)->fluxdensity) {
      free((*spec)->fluxdensity);
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


void getSimputMIdpSpecVal(const SimputMIdpSpec* const spec,
			  const long row,
			  float* const energy, 
			  float* const fluxdensity,
			  int* const status)
{
  if (row>=spec->nentries) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("row number exceeds number of bins in spectrum");
    return;
  }

  *energy     =spec->energy[row];
  *fluxdensity=spec->fluxdensity[row];
}


struct SimputMIdpSpecBuffer* newSimputMIdpSpecBuffer(int* const status)
{
  struct SimputMIdpSpecBuffer *specbuff=
    (struct SimputMIdpSpecBuffer*)malloc(sizeof(struct SimputMIdpSpecBuffer));

  CHECK_NULL_RET(specbuff, *status, 
		 "memory allocation for SimputMIdpSpecBuffer failed", specbuff);
    
  specbuff->spectrum=NULL;
  specbuff->left    =NULL;
  specbuff->right   =NULL;

  return(specbuff);
}


void freeSimputMIdpSpecBuffer(struct SimputMIdpSpecBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->spectrum) {
      freeSimputMIdpSpec(&((*sb)->spectrum));
    }
    if (NULL!=(*sb)->left) {
      freeSimputMIdpSpecBuffer(&((*sb)->left));
    }
    if (NULL!=(*sb)->right) {
      freeSimputMIdpSpecBuffer(&((*sb)->right));
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputMIdpSpec* searchSimputMIdpSpecBuffer(void* buffer, 
					   const char* const filename)
{
  struct SimputMIdpSpecBuffer* sb=(struct SimputMIdpSpecBuffer*)buffer;

  // Check if this is a leaf.
  if (NULL==sb) {
    return(NULL);
  }
  
  // Compare with the current spectrum.
  int cmp=strcmp(filename, sb->spectrum->fileref);
  if (0==cmp) {
    return(sb->spectrum);
  } else if (cmp<0) {
    return(searchSimputMIdpSpecBuffer(sb->left, filename));
  } else {
    return(searchSimputMIdpSpecBuffer(sb->right, filename));
  }
}


void insertSimputMIdpSpecBuffer(void** buffer,
				SimputMIdpSpec* const spec,
				int* const status)
{
  struct SimputMIdpSpecBuffer** sb=(struct SimputMIdpSpecBuffer**)buffer;

  // Check if this is a leaf.
  if (NULL==*sb) {
    *sb=newSimputMIdpSpecBuffer(status);
    CHECK_STATUS_VOID(*status);
    (*sb)->spectrum=spec;
  } else if (strcmp(spec->fileref, (*sb)->spectrum->fileref) < 0){
    insertSimputMIdpSpecBuffer((void**)&((*sb)->left), spec, status);
    CHECK_STATUS_VOID(*status);
  } else {
    insertSimputMIdpSpecBuffer((void**)&((*sb)->right), spec, status);
    CHECK_STATUS_VOID(*status);
  }
}


static long MIdpSpecPartition(SimputMIdpSpec** const spectra, 
			      const long left, const long right, 
			      const long pivotIndex)
{
  // Move pivot to end.
  SimputMIdpSpec* buffer;
  buffer=spectra[pivotIndex];
  spectra[pivotIndex]=spectra[right];
  spectra[right]=buffer;

  long storeIndex=left;
  long ii;
  for (ii=left; ii<right; ii++) { // left ≤ i < right  
    if (strcmp(spectra[ii]->fileref, spectra[right]->fileref)<0) {
      if (ii>storeIndex) {
	buffer=spectra[storeIndex];
	spectra[storeIndex]=spectra[ii];
	spectra[ii]=buffer;
      }
      storeIndex++;
    }
  }

  // Move pivot to its final place
  buffer=spectra[storeIndex];
  spectra[storeIndex]=spectra[right];
  spectra[right]=buffer;

  return(storeIndex);
}


static void quicksortMIdpSpec(SimputMIdpSpec** const spectra, 
			      const long left, 
			      const long right)
{
  if (right>left) {
    // Select a pivot index (e.g. pivotIndex := left+(right-left)/2)
    int pivotIndex=left+(right-left)/2;
    int pivotNewIndex=MIdpSpecPartition(spectra, left, right, pivotIndex);
    quicksortMIdpSpec(spectra, left, pivotNewIndex-1);
    quicksortMIdpSpec(spectra, pivotNewIndex+1, right);
  }
}


void buildSimputMIdpSpecBuffer(void** buffer,
			       SimputMIdpSpec** const spectra,
			       const long nspectra,
			       const int sorted,
			       int* const status)
{
  struct SimputMIdpSpecBuffer** sb=(struct SimputMIdpSpecBuffer**)buffer;

  // The binary tree must be empty.
  assert(*sb==NULL);

  // Sort the spectra according to their fileref's.
  if (0==sorted) {
    quicksortMIdpSpec(spectra, 0, nspectra-1);
  }

  // Determine the median.
  long median=(long)nspectra/2;

  // Insert the median element at the current node.
  *sb=newSimputMIdpSpecBuffer(status);
  CHECK_STATUS_VOID(*status);
  (*sb)->spectrum=spectra[median];

  if (median>0) {
    buildSimputMIdpSpecBuffer((void**)&((*sb)->left), spectra, median, 1, status);
    CHECK_STATUS_VOID(*status);
  }

  if (median<nspectra-1) {
    buildSimputMIdpSpecBuffer((void**)&((*sb)->right), &(spectra[median+1]), 
			      nspectra-1-median, 1, status);
    CHECK_STATUS_VOID(*status);
  }
}


SimputSpec* newSimputSpec(int* const status)
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


struct SimputSpecBuffer* newSimputSpecBuffer(int* const status)
{
  struct SimputSpecBuffer *specbuff=
    (struct SimputSpecBuffer*)malloc(sizeof(struct SimputSpecBuffer));

  CHECK_NULL_RET(specbuff, *status, 
		 "memory allocation for SimputSpecBuffer failed", specbuff);
    
  specbuff->spectrum=NULL;
  specbuff->left    =NULL;
  specbuff->right   =NULL;

  return(specbuff);
}


void freeSimputSpecBuffer(struct SimputSpecBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->spectrum) {
      freeSimputSpec(&((*sb)->spectrum));
    }
    if (NULL!=(*sb)->left) {
      freeSimputSpecBuffer(&((*sb)->left));
    }
    if (NULL!=(*sb)->right) {
      freeSimputSpecBuffer(&((*sb)->right));
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputSpec* searchSimputSpecBuffer(void* buffer, 
				   const char* const filename)
{
  struct SimputSpecBuffer* sb=(struct SimputSpecBuffer*)buffer;

  // Check if this is a leaf.
  if (NULL==sb) {
    return(NULL);
  }
  
  // Compare with the current spectrum.
  int cmp=strcmp(filename, sb->spectrum->fileref);
  if (0==cmp) {
    return(sb->spectrum);
  } else if (cmp<0) {
    return(searchSimputSpecBuffer(sb->left, filename));
  } else {
    return(searchSimputSpecBuffer(sb->right, filename));
  }
}


void insertSimputSpecBuffer(void** buffer,
			    SimputSpec* const spec,
			    int* const status)
{
  struct SimputSpecBuffer** sb=(struct SimputSpecBuffer**)buffer;

  // Check if this is a leaf.
  if (NULL==*sb) {
    *sb=newSimputSpecBuffer(status);
    CHECK_STATUS_VOID(*status);
    (*sb)->spectrum=spec;
  } else if (strcmp(spec->fileref, (*sb)->spectrum->fileref) < 0){
    insertSimputSpecBuffer((void**)&((*sb)->left), spec, status);
    CHECK_STATUS_VOID(*status);
  } else {
    insertSimputSpecBuffer((void**)&((*sb)->right), spec, status);
    CHECK_STATUS_VOID(*status);
  }
}


SimputLC* newSimputLC(int* const status)
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
  lc->dperiod =0.;
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


struct SimputLCBuffer* newSimputLCBuffer(int* const status)
{
  struct SimputLCBuffer *lcbuff= 
    (struct SimputLCBuffer*)malloc(sizeof(struct SimputLCBuffer));

  CHECK_NULL_RET(lcbuff, *status, 
		 "memory allocation for SimputLCBuffer failed", lcbuff);
    
  lcbuff->nlcs=0;
  lcbuff->clc =0;
  lcbuff->lcs =NULL;

  return(lcbuff);
}


void freeSimputLCBuffer(struct SimputLCBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->lcs) {
      long ii;
      for (ii=0; ii<(*sb)->nlcs; ii++) {
	freeSimputLC(&((*sb)->lcs[ii]));
      }
      free((*sb)->lcs);
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputPSD* newSimputPSD(int* const status)
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


struct SimputPSDBuffer* newSimputPSDBuffer(int* const status)
{
  struct SimputPSDBuffer *psdbuff = 
    (struct SimputPSDBuffer*)malloc(sizeof(struct SimputPSDBuffer));

  CHECK_NULL_RET(psdbuff, *status, 
		 "memory allocation for SimputPSDBuffer failed", psdbuff);
    
  psdbuff->npsds=0;
  psdbuff->psds =NULL;

  return(psdbuff);
}


void freeSimputPSDBuffer(struct SimputPSDBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->psds) {
      long ii;
      for (ii=0; ii<(*sb)->npsds; ii++) {
	freeSimputPSD(&((*sb)->psds[ii]));
      }
      free((*sb)->psds);
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputImg* newSimputImg(int* const status)
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


struct SimputImgBuffer* newSimputImgBuffer(int* const status)
{
  struct SimputImgBuffer *imgbuff = 
    (struct SimputImgBuffer*)malloc(sizeof(struct SimputImgBuffer));
  CHECK_NULL_RET(imgbuff, *status, 
		 "memory allocation for SimputImgBuffer failed", imgbuff);
    
  imgbuff->nimgs=0;
  imgbuff->imgs =NULL;

  return(imgbuff);
}


void freeSimputImgBuffer(struct SimputImgBuffer** sb)
{
  if (NULL!=*sb) {
    if (NULL!=(*sb)->imgs) {
      long ii;
      for (ii=0; ii<(*sb)->nimgs; ii++) {
	freeSimputImg(&((*sb)->imgs[ii]));
      }
      free((*sb)->imgs);
    }
    free(*sb);
    *sb=NULL;
  }
}


SimputPhList* newSimputPhList(int* const status)
{
  SimputPhList* phl=(SimputPhList*)malloc(sizeof(SimputPhList));
  CHECK_NULL_RET(phl, *status, 
		 "memory allocation for SimputPhList failed", phl);

  // Initialize elements.
  phl->fptr    =NULL;
  phl->nphs    =0;
  phl->nrphs   =0;
  phl->cra     =0;
  phl->cdec    =0;
  phl->cenergy =0;
  phl->ctime   =0;
  phl->fra     =0.;
  phl->fdec    =0.;
  phl->fenergy =0.;
  phl->ftime   =0.;
  phl->refarea =0.;
  phl->accrate =0.;
  phl->mjdref  =0.;
  phl->timezero=0.;
  phl->tstart  =0.;
  phl->tstop   =0.;
  phl->currrow =0;
  phl->fileref=NULL;

  return(phl);
}


void freeSimputPhList(SimputPhList** const phl, int* const status)
{
  if (NULL!=*phl) {
    if (NULL!=(*phl)->fileref) {
      free((*phl)->fileref);
    }
    if (NULL!=(*phl)->fptr) {
      fits_close_file((*phl)->fptr, status);
    }
    free(*phl);
    *phl=NULL;
  }
}


struct SimputPhListBuffer* newSimputPhListBuffer(int* const status)
{
  struct SimputPhListBuffer *phlbuff = 
    (struct SimputPhListBuffer*)malloc(sizeof(struct SimputPhListBuffer));
  CHECK_NULL_RET(phlbuff, *status, 
		 "memory allocation for SimputPhListBuffer failed", phlbuff);
    
  phlbuff->nphls=0;
  phlbuff->phls =NULL;

  return(phlbuff);
}


void freeSimputPhListBuffer(struct SimputPhListBuffer** pb,
			    int* const status)
{
  if (NULL!=*pb) {
    if (NULL!=(*pb)->phls) {
      long ii;
      for (ii=0; ii<(*pb)->nphls; ii++) {
	freeSimputPhList(&((*pb)->phls[ii]), status);
      }
      free((*pb)->phls);
    }
    free(*pb);
    *pb=NULL;
  }
}


