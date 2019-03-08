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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "simputimg.h"


int simputimg_main()
{
  // Program parameters.
  struct Parameters par;

  // Output SimputImg.
  SimputImg* simputimg=NULL;

  // SimputCtlg the image should be attached to.
  SimputCtlg* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputimg");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputimg_getpar(&par);
    CHECK_STATUS_BREAK(status);

    if ((0==strcmp(par.ImageFile, "none"))||
	(0==strcmp(par.ImageFile, "NONE"))) {
      strcpy(par.ImageFile, "");
    }
    if (strlen(par.ImageFile)==0) {
      SIMPUT_ERROR("no input image specified");
      status=EXIT_FAILURE;
      break;
    }
    // END of checking the input type for the image.

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Load the image from the input file.
    simputimg=loadSimputImg(par.ImageFile, &status);
    CHECK_STATUS_BREAK(status);

    // Store the image in the SIMPUT file.
    saveSimputImg(simputimg, par.Simput, par.Extname, par.Extver, &status);
    CHECK_STATUS_BREAK(status);

    // Open the SimputCtlg.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Set the image reference in the source catalog.
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
    char* imgref=(char*)malloc(32*sizeof(char));
    CHECK_NULL_BREAK(imgref, status, "memory allocation failed");
    sprintf(imgref, "[%s,%d]", par.Extname, par.Extver);
    fits_write_col(cat->fptr, TSTRING, cat->cimage, 1, 1, 1, &imgref, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Release memory.
  freeSimputImg(&simputimg);
  freeSimputCtlg(&cat, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int simputimg_getpar(struct Parameters* const par)
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

  status=ape_trad_query_string("ImageFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the image file failed");
    return(status);
  }
  strcpy(par->ImageFile, sbuffer);
  free(sbuffer);

  return(status);
}
