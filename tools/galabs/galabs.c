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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include <stdio.h>
#include <string.h>
#include <pil.h>
#include <math.h>
#include <assert.h>
#include <simput.h>
#include "labnh.h"
#include "spectree.h"

#include <time.h>

#define MAXSTRING (2048)

int searchlev=0;

typedef struct{
  char cat[MAXSTRING];
  char tempfile[MAXSTRING];
  char ocat[MAXSTRING];
  char ospec[MAXSTRING];
  int catspec;    //=0: both in one file, =1: seperate files
  float NHmin;
  float NHmax;
  float NHstep;
  int verbosity;
  int sampling;
}par;

static float *nhgrid=NULL;	// nh-grid for output
static int nhints;		// number of nh-intervals
static float nhstep;		// stepwidth in nh-grid

static SimputSrc **src=NULL;	// simput source block
const long long blockl=100000;	// length of source block
static long long donesrc=0, srcnmbr=0;	//number of sources in source block/total

static SimputMIdpSpec **spc=NULL;	// simut spec block
const long long sblockl=100000;		// length of spec block
static long long donespec=0;		// number of spectra in spec block

static int specfilenentries=0;		// number of files with spectra
static char **specfiles=NULL;		// array of input files with spectra
struct speclist **specs=NULL;		// list of absorbed spectra
struct speclist **oldspecs=NULL;	// list with input spectra

static int NHsteps;		// number of internal nh steps
static int Esteps;		// number of internal energy steps
static float* NHarray=NULL;	// grid for nh (internal)
static float* Earray=NULL;	// grid for energy (internal)
static float** absorptionrel=NULL;	// matrix with flux relations between unabsorbed and absorbed fluxes

void xsphab_(const float* energyArray,     //the energy ranges on which to calculate the model
	     const int *Nenergy,           //the number of energy ranges
	     const float* parameterValues, //the H column density in 10^22 cm^-2
	     const int *spectrumNumber,    //the dataset
	     float* flux,                  //fractional transmission
	     float* fluxError);            //

void FNINIT();
void FPSOLR(const char* table, int* ierr);

static long long newload=0;		// counts number of loaded spectra
static long long cacheload=0;		// counts number of buffered spectra calls
static long long newabsorbed=0;		// counts number of absorbed spectra
static long long cacheabsorbed=0;	// counts number of absorbed buffered spectra calls



int SimputGetSpecFilename(char *specref, char **filename){
//This function determines the name of the file, in which
//the spectrum pointed by specref is included.
//If the spectrum is in the same file as the source and no
//filename is given according to the Extended Filename Syntax,
//the function returns 1. 0 is returned if a filename is found.
//An error will cause -1 as return value.

	if(specref==NULL){
		return -1;
	}

	*filename=(char*)malloc((1+(int)strlen(specref))*sizeof(char));
	strcpy(*filename, specref);

	char *search=strchr(*filename, '[');
	*search='\0';

	if(search==*filename){
        	return 1;
	}

	else{
        return 0;
	}

}

int SimputGetSpecFileExt(char *specref, char **filename){
//This function determines the name of the file, in which
//the spectrum pointed by specref is included.
//If the spectrum is in the same file as the source and no
//filename is given according to the Extended Filename Syntax,
//the function returns 1. 0 is returned if a filename is found.
//An error will cause -1 as return value.

	if(specref==NULL){
		return -1;
	}

	*filename=(char*)malloc((1+(int)strlen(specref))*sizeof(char));
	strcpy(*filename, specref);

	char *search=strchr(*filename, ']');
	*(search+1)='\0';

	if(search==*filename){
        	return 1;
	}

	else{
        return 0;
	}

}

int SimputGetSpecName(char* specref, char** name){
//This function determines the name of the spectrum
//pointed by specref.
//If no name is found, 1 is returned. An error results
//in -1. If a name was found, 0 is the return value.

	if(specref==NULL){
		return -1;
	}

	char* refcpy=(char*)malloc((1+strlen(specref))*sizeof(char));
	strcpy(refcpy, specref);

	char* search=strstr(refcpy, "NAME==");

	if(search==NULL){
        	return 1;
	}

	search=search+6;
	char quot=*search;
	search++;

	char *nameend=strchr(search, quot);
	*nameend='\0';

	*name=(char*)malloc((1+(int)strlen(search))*sizeof(char));
	strcpy(*name, search);

	return 0;

}

int find_grid(float *array, int length, float number){
//this function returns the index, where "number" is found in "array".
//When "number" is outside the array, the index of the closer border is returned.

  int low, high, mid, diff=100;
  low=0;
  high=length-1;
  mid=((high-low)/2)+low;

  if(array[high]<=number){
    return high;
    }
  else if(array[low]>=number){
    return low;
    }

  while(diff!=1){

    diff=high-low;

    if(array[mid]>number){
      high=mid;
      }

    else{
      low=mid;
      }
    mid=(high+low)/2;
  }

  return low;
}


void absorptionInit(float NHmax, float NHmin, float dNH, float Emin, float Emax, float dE){
// Precalculates absorption ratios on nh-E-grid

  int ii, jj;
  float *Elo=NULL;
  float *Ehi=NULL;
  float *absorption=NULL;
  float *absorptionratio=NULL;
  float *earr=NULL;
  float snh, fsnh;
  NHsteps=(int)((NHmax-NHmin)/dNH)+1;
  Esteps=(int)((Emax-Emin)/dE);
  nhstep=dNH;

  Earray=(float*)malloc(Esteps*sizeof(float));
  assert(Earray);
  Elo=(float*)malloc(Esteps*sizeof(float));
  assert(Elo);
  Ehi=(float*)malloc(Esteps*sizeof(float));
  assert(Ehi);
  for(ii=0; ii<Esteps; ii++){
    Elo[ii]=Emin+ii*dE;
    Ehi[ii]=Emin+(ii+1)*dE;
    Earray[ii]=(Elo[ii]+Ehi[ii])/2.;
    }
  absorption=(float*)malloc(Esteps*sizeof(float));
  assert(absorption);
  absorptionratio=(float*)malloc(Esteps*sizeof(float));
  assert(absorptionratio);
  earr=(float*)malloc((Esteps+1)*sizeof(float));
  assert(earr);
  earr[0]=Emin;
  for(ii=0; ii<Esteps; ii++){
    earr[ii+1]=Ehi[ii];
    }

  NHarray=(float*)malloc(NHsteps*sizeof(float));
  assert(NHarray);
  absorptionrel=(float**)malloc(NHsteps*sizeof(float*));
  assert(absorptionrel);
  for(ii=0; ii<NHsteps; ii++){
    absorptionrel[ii]=(float*)malloc(Esteps*sizeof(float));
    assert(absorptionrel[ii]);
    }

  //calculate absorption
  for(ii=0; ii<NHsteps; ii++){
    snh=pow(10., (NHmin+(float)ii*dNH));
    fsnh=snh/1.e22;
    printf("NH=%.5e\n", snh);
    xsphab_(earr, &Esteps, &fsnh, &Esteps, absorption, absorptionratio);
    NHarray[ii]=(NHmin+(float)ii*dNH);
    for(jj=0; jj<Esteps; jj++){
      absorptionrel[ii][jj]=absorption[jj];
      }
    }
}

void Specabsorb(int nhind, struct speclist *ispec, SimputMIdpSpec** ospec){
// generates a new absorbed spectrum and puts it into the tree

  //int nhind;
  int eind, ii;
  float abs;
  assert(*ospec);
  (*ospec)->nentries=ispec->nentries;
  (*ospec)->energy=(float*)malloc(ispec->nentries*sizeof(float));
  if((*ospec)->energy==NULL){
    puts("Memory allocation for ospec energy array failed.");
    exit(1);
    }
  (*ospec)->pflux=(float*)malloc(ispec->nentries*sizeof(float));
  if((*ospec)->pflux==NULL){
    puts("Memory allocation for ospec pflux array failed.");
    exit(1);
    }

  //nhind=find_grid(NHarray, NHsteps, (float)lognh);
  eind=find_grid(Earray, Esteps, ispec->energy[0]);

  for(ii=0; ii<ispec->nentries; ii++){
    while(Earray[eind]>ispec->energy[ii]){
      eind++;
      }
    abs=absorptionrel[nhind][eind]+(ispec->energy[ii]-Earray[eind])/(Earray[eind+1]-Earray[eind])*(absorptionrel[nhind][eind+1]-absorptionrel[nhind][eind]);
    (*ospec)->energy[ii]=ispec->energy[ii];
    (*ospec)->pflux[ii]=ispec->pflux[ii]*abs;

    }
}

void getSpecfile(char* reference, char **filename){
// returns the filename which is included in a simput spectrum reference

  int length=0;
  int found=0;
  int ii;

  while(reference[length]!='\0' && found!=1){
    if(reference[length]=='[' && reference[length+1]=='N' && reference[length+2]=='A' && reference[length+3]=='M' && reference[length+4]=='E'){
      found=1;
      }
    else
      length++;
    }

  *filename=(char*)malloc((length+1)*sizeof(char));
  if(*filename==NULL){
    puts("ERROR: unable to allocate memory for filename.");
    exit(1);
    }

  for(ii=0; ii<length; ii++){
    (*filename)[ii]=reference[ii];
    }
  (*filename)[length]='\0';

}

int findSpecfile(char* filename){
//returns the index of the file in the buffer

  if(specfilenentries==0)
    return -1;

  int fileindex=0;
  while(fileindex<specfilenentries && strcmp(specfiles[fileindex], filename)!=0){
    fileindex++;
    }

  if(fileindex==specfilenentries){
    fileindex=-1;
    printf("findSpecfile: specfiles[0]=%s\n", specfiles[0]);
  }

  return fileindex;

}

void insertNewSpecfile(char* filename){
// inserts a new input file into the buffer

  specfilenentries++;

  specfiles=(char**)realloc(specfiles, specfilenentries*sizeof(char*));
  if(specfiles==NULL){
    puts("ERROR: wasn't able to (re)allocate memory for specfiles-array");
    exit(1);
    }

  specfiles[specfilenentries-1]=(char*)malloc(MAXSTRING*sizeof(char));
  if(specfiles[specfilenentries-1]==NULL){
    puts("ERROR: wasn't able to allocate memory for new specfilename");
    exit(1);
    }
  strcpy(specfiles[specfilenentries-1], filename);

  specs=(struct speclist**)realloc(specs, specfilenentries*sizeof(struct speclist*));
  if(specs==NULL){
    puts("ERROR: wasn't able to (re)allocate memory for specs-array");
    exit(1);
    }
  specs[specfilenentries-1]=return_nil();

  oldspecs=(struct speclist**)realloc(oldspecs, specfilenentries*sizeof(struct speclist*));
  if(oldspecs==NULL){
    puts("ERROR: wasn't able to (re)allocate memory for oldspecs-array");
    exit(1);
    }
  oldspecs[specfilenentries-1]=return_nil();

}

void getParameters(int argc, char **argv, par *pars){
//reads input parameters via PIL
  int parameter_read;

  if(PIL_OK!=(parameter_read=PILInit(argc, argv))){
    printf("ERROR: PILInit failed. %s\n", PIL_err_handler(parameter_read));
    exit(1);
    }
  puts("PIL successfully initialized.\n");

  if(PILGetString("input_catalog", pars->cat)<0){
    printf("ERROR: wasn't able to read input catalog name.\n");
    exit(1);
    }
  printf("input catalog name = \"%s\" .\n", pars->cat);

  if(PILGetInt("twofiles", &pars->catspec)<0){
    printf("ERROR: wasn't able to read twofiles value.\n");
    exit(1);
    }
  printf("twofiles = \"%d\" .\n", pars->catspec);

  if(PILGetString("output_catalog", pars->ocat)<0){
    printf("ERROR: wasn't able to read filename for output catalog.\n");
    exit(1);
    }
  printf("output catalog name = \"%s\" .\n", pars->ocat);
  if(strcmp(pars->cat, pars->ocat)==0){
    puts("invalid name for output catalog.\ninput and output must not have the same name.");
    }

  if(pars->catspec==1){
    if(PILGetString("output_spec", pars->ospec)<0){
      printf("ERROR: wasn't able to read filename for output spectra.\n");
      exit(1);
      }
      printf("filename for output spectra = \"%s\" .\n", pars->ospec);
    }
  else{
    strcpy(pars->ospec, pars->ocat);
    if(pars->ospec==NULL){
      printf("ERROR: wasn't able to determine pars->ospec\n");
      exit(1);
      }
    }

  if(pars->catspec==0){
    if(PILGetString("temp_filename", pars->tempfile)<0){
      printf("ERROR: wasn't able to read filename for temporary data.\n");
      exit(1);
      }
    }
  else{
    strcpy(pars->tempfile, pars->ospec);
    }
  printf("filename for temporary data = \"%s\" .\n", pars->tempfile);


  if(PILGetReal4("logNH_min", &pars->NHmin)<0){
    printf("ERROR: wasn't able to read minimum log NH value.\n");
    exit(1);
    }
  printf("logNH_min = \"%f\" .\n", pars->NHmin);

  if(PILGetReal4("logNH_max", &pars->NHmax)<0){
    printf("ERROR: wasn't able to read maximum log NH value.\n");
    exit(1);
    }
  printf("logNH_max = \"%f\" .\n", pars->NHmax);

  if(PILGetReal4("dlogNH", &pars->NHstep)<0){
    printf("ERROR: wasn't able to read log NH grid width.\n");
    exit(1);
    }
  printf("dlogNH = \"%f\" .\n", pars->NHstep);

  if(PILGetInt("sampling", &pars->sampling)<0){
    printf("ERROR: wasn't able to read sampling value.\n");
    exit(1);
    }
  printf("sampling = \"%d\" .\n", pars->sampling);

  if(PILGetInt("verbosity", &pars->verbosity)<0){
    printf("ERROR: wasn't able to read verbosity value.\n");
    exit(1);
    }
  printf("verbosity = \"%d\" .\n", pars->verbosity);

}

void NHgridInit(float NHmin, float NHmax, float NHstep){
//initializes nhgrid

  int ii;

  nhints=(int)ceil((NHmax-NHmin)/NHstep);
  nhgrid=(float*)malloc(nhints*sizeof(float));

  for(ii=0; ii<nhints; ii++){
    nhgrid[ii]=NHmin+(float)ii*NHstep;
    //printf("nhgrid[%d]=%.3e\n", ii, nhgrid[ii]);
    }
}

double integrateBSpec(long nentries, float *energy, float *pflux, float a, float b){
  //This function integrates the the spectrum from a to b
  static int isav=0;
  if ( isav>=nentries || !(energy[isav]<=a && energy[isav+1]>a)) {
  	isav=find_grid(energy, nentries, a);
   	}
  int ii=isav;
  double integral=0.;
  while(energy[ii+1]<=b && ii<nentries-2){
    integral+=(energy[ii+1]-energy[ii])*(pflux[ii]+pflux[ii+1])*(energy[ii]+energy[ii+1]);
    ii++;
   }
   integral=integral*1.6022e-9/4.; // adjust erg->eV
   return integral;
}

double integrateSpec(SimputMIdpSpec* spec, float a, float b){
   return integrateBSpec(spec->nentries, spec->energy, spec->pflux, a, b);
}

void doSrc(SimputCtlg* icat, SimputCtlg** ocat, int srcind, char *specfilename, char *tempfilename, int *status){
// this function absorbes one source at a time

  long long jj;
  double brightness_old, brightness_new;
  char *specfile=NULL;
  int exist;

  //allocate memory for source block
  if(src==NULL){
    puts("Allocate memory for source block...");
    src=(SimputSrc**)malloc(blockl*sizeof(SimputSrc*));
    if(src==NULL){
      puts("Memory allocation for Source Block failed");
      exit(1);
      }
    }

  //allocate memory for spec block
  if(spc==NULL){
    puts("Allocate memory for spec block...");
    spc=(SimputMIdpSpec**)malloc(sblockl*sizeof(SimputMIdpSpec*));
    if(spc==NULL){
      puts("Memory allocation for Spec Block failed");
      exit(1);
      }
    }

  //load source, spectrum and construct new name
  src[donesrc]=loadSimputSrc(icat, srcind+1, status);
  fits_report_error(stderr, *status);
  //evaluates NH value on source position
  double ra=360./2./M_PI*(float)src[donesrc]->ra;
  double dec=360./2./M_PI*(float)src[donesrc]->dec;
  double snh=nh_equ(ra, dec);
  double lognh=log10(snh);

  //NH-grid
  int gridind=(int)((lognh-nhgrid[0])/nhstep);
  lognh=nhgrid[gridind];

  //get name of spectrum file
  SimputGetSpecFileExt(src[donesrc]->spectrum, &specfile);
  //check if file is already known
  int fileknown=findSpecfile(specfile);
  //if file is not yet known, append to array
  if(fileknown==-1){
    puts("new file detected.");
    puts(specfile);
    insertNewSpecfile(specfile);
    //load all spectra of file
    loadCacheAllSimputMIdpSpec(icat, specfile, status);
    fits_report_error(stderr, *status);
    fileknown=specfilenentries-1;
    }

  char* inspecname=NULL;
  SimputGetSpecName(src[donesrc]->spectrum, &inspecname);

  struct speclist* findoldspec=check_if_exists(oldspecs[fileknown], inspecname, &exist);
  searchlev=0;

  SimputMIdpSpec* ispec =NULL;
  char specname[2*MAXSTRING], specref[2*MAXSTRING], tempspec[2*MAXSTRING];

  if(exist==1){
    cacheload++;
    if(findoldspec->emin==src[donesrc]->e_min && findoldspec->emax==src[donesrc]->e_max){
      brightness_old=findoldspec->integral;
      }
    else{
      brightness_old=integrateBSpec(findoldspec->nentries, findoldspec->energy, findoldspec->pflux, src[donesrc]->e_min, src[donesrc]->e_max);
      }
    }
  else{
    newload++;
    ispec=getSimputSrcMIdpSpec(icat, src[donesrc], 0., 0., status);
    brightness_old=integrateSpec(ispec, src[donesrc]->e_min, src[donesrc]->e_max);
    oldspecs[fileknown]=insert_spec(oldspecs[fileknown], ispec->name, src[donesrc]->spectrum, src[donesrc]->e_min, src[donesrc]->e_max, brightness_old, ispec->nentries, ispec->energy, ispec->pflux);
    findoldspec=check_if_exists(oldspecs[fileknown], inspecname, &exist);
    if(exist!=1){
      puts("Insertation of spectrum failed. Exit.");
      exit(1);
      }
    }
  searchlev=0;

  sprintf(specname, "f%d%sgal%2.3f", fileknown, inspecname, lognh);
  sprintf(specref,"%s[SPEC][NAME=='%s']", specfilename, specname);
  sprintf(tempspec,"%s[SPEC][NAME=='%s']", tempfilename, specname);
  struct speclist* findspec=check_if_exists(specs[fileknown], specname, &exist);
  searchlev=0;

  //if spectrum already exists, take its integral or integrate again
  if(exist==1){
    cacheabsorbed++;
    if(findspec->emin==src[donesrc]->e_min && findspec->emax==src[donesrc]->e_max){
      brightness_new=findspec->integral;
      }
    else{
      brightness_new=integrateBSpec(findspec->nentries, findspec->energy, findspec->pflux, src[donesrc]->e_min, src[donesrc]->e_max);
      }
    }
  //if spectrum not exists, make it
  else{
    spc[donespec]=newSimputMIdpSpec(status);
    if(spc[donespec]==NULL){
      puts("malloc failed. exit.");
      exit(1);
    }
    spc[donespec]->name=(char*)malloc(2*MAXSTRING*sizeof(char));
    if(spc[donespec]->name==NULL){
      puts("Memory allocation for spectrum name failed.");
      exit(1);
    }
    strcpy(spc[donespec]->name, specname);
    //absorb
    newabsorbed++;
    Specabsorb(gridind, findoldspec, &(spc[donespec]));

    assert(spc[donespec]->pflux);
    brightness_new=integrateSpec(spc[donespec], src[donesrc]->e_min, src[donesrc]->e_max);
    donespec++;

    specs[fileknown]=insert_spec(specs[fileknown], spc[donespec-1]->name, specref, src[donesrc]->e_min, src[donesrc]->e_max, brightness_new, spc[donespec-1]->nentries, spc[donespec-1]->energy, spc[donespec-1]->pflux);
    findspec=check_if_exists(specs[fileknown], spc[donespec-1]->name, &exist);
    if(exist!=1){
      puts("Insertation of spectrum failed. Exit.");
      exit(1);
    }
  }
  searchlev=0;
  src[donesrc]->eflux=src[donesrc]->eflux*brightness_new/brightness_old;
  free(src[donesrc]->spectrum);
  src[donesrc]->spectrum=(char*)malloc(2*MAXSTRING*sizeof(char));
  if(src[donesrc]->spectrum==NULL){
    puts("Memory allocation for src spec reference failed.");
    exit(1);
    }
  strcpy(src[donesrc]->spectrum, specref);
  srcnmbr++;
  donesrc++;

  fits_report_error(stderr, *status);
// if spec block is full, save it and free the individuals
  if(donespec==sblockl){
    saveSimputMIdpSpecBlock(spc, donespec, tempfilename, "SPEC", 1, status);
    fits_report_error(stderr, *status);
    puts("SpecBlock appended successfully.");
    for(jj=0; jj<sblockl; jj++){
      if(spc[jj]==NULL){
        printf("spc[%lld] was not valid. Appending failed.", jj);
        exit(1);
      }
      freeSimputMIdpSpec(&(spc[jj]));
      spc[jj]=NULL;
    }
    donespec=0;
  }
  fits_report_error(stderr, *status);
  if(*status!=0){
    puts("exit");
    exit(1);
    }
// if source block is full and catalog is not yet finished, free the source block
  if(donesrc==blockl && srcnmbr!=icat->nentries-1){
    donesrc=0;
    puts("begin append...");
    appendSimputSrcBlock(*ocat, src, blockl, status);
    fits_report_error(stderr, *status);
    puts("appended successfully.");
    for(jj=0; jj<blockl; jj++){
      if(src[jj]==NULL){
        printf("src[%lld] was not valid. Appending failed.", jj);
        exit(1);
        }
      freeSimputSrc(&src[jj]);
      }
    }
// if the catalog is finished, save the remaining spurces and free the simput blocks
  if(srcnmbr==icat->nentries){
    if(donesrc!=0){
      appendSimputSrcBlock(*ocat, src, donesrc, status);
      fits_report_error(stderr, *status);
      }
    for(jj=0; jj<blockl; jj++) {
      freeSimputSrc(&src[jj]);
      }
    free(src);
    puts("All sources appended.");
    if(donespec>0){
      saveSimputMIdpSpecBlock(spc, donespec, tempfilename, "SPEC", 1, status);
      fits_report_error(stderr, *status);
      puts("SpecBlock appended successfully.");
      for(jj=0; jj<donespec; jj++){
        if(spc[jj]==NULL){
          printf("spc[%lld] was not valid. Appending failed.", jj);
          exit(1);
        }
        freeSimputMIdpSpec(&spc[jj]);
      }
    }
    free(spc);
    puts("All spectra appended");
  }
}

int main (int argc,char **argv) {

  //define variables
  int status=0;
  int ii;
  long nentries;
  par pars;
  SimputCtlg* rawcatalog=NULL;

  puts("Read Parameters...");
  getParameters(argc, argv, &pars);

  puts("Set nh_verbosity and sampling...");
  nh_verbosity(pars.verbosity);
  nh_sampling(pars.sampling);

  puts("Initialize NH grid");
  NHgridInit(pars.NHmin, pars.NHmax, pars.NHstep);

  puts("Do FNINIT...");
  // initialize xspec routines
  int ierr=0;
  FNINIT();

  puts("Set abundances...");
  // Set abundance set to Wilms abundances
  FPSOLR("wilm",&ierr);

  //initialize absorption arrays
  absorptionInit(pars.NHmax, pars.NHmin, pars.NHstep, 0.001, 1100., 0.01);

  puts("Open input Catalog...");
  //load source catalog
  rawcatalog = openSimputCtlg(pars.cat, READONLY, 0, 0, 0, 0, &status);
  fits_report_error(stderr, status);
  nentries=rawcatalog->nentries;

  puts("Open output Catalog...");
  //open output file
  SimputCtlg *ocat;
  ocat=openSimputCtlg(pars.ocat, READWRITE, 0, 0, 0, 0, &status);
  fits_report_error(stderr, status);

  nil_init();

  puts("Main loop...");
  //main loop: gets NH-values for source positions,
  //loads spectra, calculates absorption, absorbes
  for(ii=0; ii<nentries; ii++){
	searchlev=0;
	doSrc(rawcatalog, &ocat, ii, pars.ospec, pars.tempfile, &status);
    if(ii%(nentries/100)==0){
      printf("%d/%ld (%d %%) sources absorbed.\n", ii, nentries, (int)(ii*100/nentries));
    }
  }
  printf("Spectra rewritten.\n");

  //close simput files
  freeSimputCtlg(&ocat, &status);
  freeSimputCtlg(&rawcatalog, &status);

  fits_report_error(stderr, status);
  puts("copy spectra into file...");

  //copy spectra into output file
  if(pars.catspec==0){
    char tempfileext[MAXSTRING*3];
    sprintf(tempfileext, "%s[SPEC]", pars.tempfile);

    fitsfile *ifptr;
    if(fits_open_file(&ifptr, tempfileext, READONLY, &status)){
      fits_report_error(stderr, status);
      printf("ERROR: Tempfile \"%s\" does not exist.\n", tempfileext);
      exit(1);
      }
    fitsfile *ofptr;

    if(fits_open_file(&ofptr, pars.ocat, READWRITE, &status)){
      fits_report_error(stderr, status);
      printf("ERROR: Output file \"%s\" can't be opened.\n", pars.ocat);
      exit(1);
      }

    if(fits_copy_hdu(ifptr, ofptr, 0, &status)){
      fits_report_error(stderr, status);
      }

    fits_close_file(ifptr, &status);
    fits_close_file(ofptr, &status);

    //delete temporary file

    if(remove(pars.tempfile)!=0)printf("ERROR deleting temporary file\"%s\"\n", pars.tempfile);

  }

  fitsfile *ofptr;
  char oext[3*MAXSTRING];
  sprintf(oext,"%s[SRC_CAT]", pars.ocat);
  puts("test status...");
  fits_report_error(stderr, status);
  puts("status tested.");

  if(fits_open_file(&ofptr, oext, READWRITE, &status)){
      fits_report_error(stderr, status);
      printf("ERROR: File \"%s\" does not exist.\n", oext);
      exit(1);
      }

  fits_write_comment(ofptr, "Galactic absorption with galabs. N_H value from Leiden/Argentine/Bonn Galactic HI Survey", &status);
  fits_update_key(ofptr, TFLOAT, "GalNHmin", &pars.NHmin, "minimum of GalNH-grid", &status);
  fits_update_key(ofptr, TFLOAT, "GalNHmax", &pars.NHmax, "maximum of GalNH-grid", &status);
  fits_update_key(ofptr, TFLOAT, "GalNHstep", &pars.NHstep, "stepsize of GalNH-grid", &status);
  fits_close_file(ofptr, &status);

  printf("Number of input spectrum loads:   %lld\n", newload);
  printf("Number of call of cached spectra: %lld\n", cacheload);
  printf("Number of absorbed Spectra:       %lld\n", newabsorbed);
  printf("Number of call of cached\n              absorbed spectra:   %lld\n\n", cacheabsorbed);

}
