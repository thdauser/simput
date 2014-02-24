/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "simputlc.h"


int simputlc_main() 
{
  // Program parameters.
  struct Parameters par;

  // Input ASCII file containing the light curve.
  FILE* asciilc=NULL;

  // Output SimputLC.
  SimputLC* simputlc=NULL;

  // SimputCtlg the light curve should be attached to.
  SimputCtlg* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputlc");
  set_toolversion("0.02");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputlc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Open the ASCII file with the light curve.
    asciilc=fopen(par.LCFile,"r");
    CHECK_NULL_BREAK(asciilc, status, "could not open input light curve file");

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Determine the number of rows.
    long nlines=0;
    char c=0;
    while(!feof(asciilc)) {
      c=fgetc(asciilc);
      if ('\n'==c) {
	nlines++;
      }
    }
    // Check if the last line has been empty.
    if('\n'==c) {
      nlines--;
    }

    // Allocate memory.
    simputlc=newSimputLC(&status);
    CHECK_STATUS_BREAK(status);
    simputlc->nentries=nlines;
    simputlc->time=(double*)malloc(nlines*sizeof(double));
    CHECK_NULL_BREAK(simputlc->time, status, "memory allocation failed");
    simputlc->flux=(float*)malloc(nlines*sizeof(float));
    CHECK_NULL_BREAK(simputlc->flux, status, "memory allocation failed");

    // Reset the file pointer, read the data and store them in
    // the SimputLC data structure.
    rewind(asciilc);
    long ii;
    for (ii=0; ii<nlines; ii++) {
      if (fscanf(asciilc, "%lf %f\n",
		 &(simputlc->time[ii]), 
		 &(simputlc->flux[ii]))<2) {
	SIMPUT_ERROR("failed reading light curve from ASCII file");
	status=EXIT_FAILURE;
	break;
      }
    }
    CHECK_STATUS_BREAK(status);

    // Store the light curve in the SIMPUT file.
    saveSimputLC(simputlc, par.Simput, par.Extname, par.Extver, &status);
    CHECK_STATUS_BREAK(status);


    // Open the SimputCtlg.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Set the timing reference in the source catalog.
    if (strlen(par.Extname)==0) {
      SIMPUT_ERROR("no EXTNAME specified");
      status=EXIT_FAILURE;
      break;
    }
    if (strlen(par.Extname)>24) {
      SIMPUT_ERROR("EXTNAME too long");
      status=EXIT_FAILURE;
      break;
    }
    if ((par.Extver<=0) || (par.Extver>9999)) {
      char msg[SIMPUT_MAXSTR];
      sprintf(msg, "value for EXTVER outside of allowed limit (%d)", par.Extver);
      SIMPUT_ERROR(msg);
      status=EXIT_FAILURE;
      break;
    }
    char* timeref=(char*)malloc(32*sizeof(char));
    CHECK_NULL_BREAK(timeref, status, "memory allocation failed");
    sprintf(timeref, "[%s,%d]", par.Extname, par.Extver);
    fits_write_col(cat->fptr, TSTRING, cat->ctiming, 1, 1, 1, 
    		   &timeref, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Close open files.
  if (NULL!=asciilc) {
    fclose(asciilc);
    asciilc=NULL;
  }

  // Release memory.
  freeSimputLC(&simputlc);
  freeSimputCtlg(&cat, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int simputlc_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.
  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the SIMPUT catalog failed");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Extname", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the EXTNAME of the generated HDU failed");
    return(status);
  }
  strcpy(par->Extname, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("Extver", &par->Extver);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the EXTVER of the generated HDU failed");
    return(status);
  }

  status=ape_trad_query_file_name("LCFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the light curve file failed");
    return(status);
  }
  strcpy(par->LCFile, sbuffer);
  free(sbuffer);

  return(status);
}


