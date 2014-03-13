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

#include "simputfile.h"


int simputfile_main()
{
  // Program parameters.
  struct Parameters par;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputfile");
  set_toolversion("0.19");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputfile_getpar(&par);
    CHECK_STATUS_BREAK(status);

    if ((0==strcmp(par.LCFile, "none"))||
	(0==strcmp(par.LCFile, "NONE"))) {
      strcpy(par.LCFile, "");
    }
    if ((0==strcmp(par.PSDFile, "none"))||
	(0==strcmp(par.PSDFile, "NONE"))) {
      strcpy(par.PSDFile, "");
    }
    if ((0==strcmp(par.ImageFile, "none"))||
	(0==strcmp(par.ImageFile, "NONE"))) {
      strcpy(par.ImageFile, "");
    }

    // Check the input type for the power spectrum: individual 
    // components or an ASCII file. Only one of these two option 
    // may be used. In case multiple of them exist, throw an error 
    // message and abort.

    int ntoptions=0;
    if (strlen(par.LCFile)>0) {
      ntoptions++;
    }
    if (strlen(par.PSDFile)>0) {
      ntoptions++;
    }
    if ((par.LFQ!=0) || (par.HBOQ!=0) || 
	(par.Q1Q!=0) || (par.Q2Q!=0) || (par.Q3Q!=0)) {
      ntoptions++;
    }
    if (ntoptions>1) {
      SIMPUT_ERROR("specification of multiple timing models not possible");
      status=EXIT_FAILURE;
      break;
    }
    // END of checking the input type for the timing extension.

    // ---- END of Initialization ----


    // ---- Main Part ----

    char command[SIMPUT_MAXSTR];
    char sclobber[SIMPUT_MAXSTR], shistory[SIMPUT_MAXSTR];
    if (0==par.clobber) {
      strcpy(sclobber, "no");
    } else {
      strcpy(sclobber, "yes");
    }
    if (0==par.history) {
      strcpy(shistory, "no");
    } else {
      strcpy(shistory, "yes");
    }

    // Call 'simputsrc' to produce a SIMPUT catalog with a 
    // single point-like source.
    sprintf(command, 
	    "simputsrc Simput=%s Src_ID=%d Src_Name=%s RA=%f Dec=%f "
	    "Emin=%f Emax=%f Flux=%e chatter=%d clobber=%s history=%s",
	    par.Simput, par.Src_ID, par.Src_Name, par.RA, par.Dec,
	    par.Emin, par.Emax, par.srcFlux,
	    par.chatter, sclobber, shistory);
    status=system(command);
    CHECK_STATUS_BREAK(status);

    // Call 'simputspec' to produce a spectrum and assign it
    // to the source.
    sprintf(command, 
	    "simputspec Simput=%s "
	    "plPhoIndex=%f plFlux=%e "
	    "bbkT=%f bbFlux=%e "
	    "flSigma=%e flFlux=%e "
	    "rflSpin=%f rflFlux=%e "
	    "NH=%e "
	    "Emin=%f Emax=%f "
	    "ISISFile=%s XSPECFile=%s PHAFile=%s "
	    "chatter=%d history=%s",
	    par.Simput,
	    par.plPhoIndex, par.plFlux,
	    par.bbkT, par.bbFlux,
	    par.flSigma, par.flFlux,
	    par.rflSpin, par.rflFlux,
	    par.NH,
	    par.Emin, par.Emax, 
	    par.ISISFile, par.XSPECFile, par.PHAFile,
	    par.chatter, shistory);
    status=system(command);
    CHECK_STATUS_BREAK(status);

    // Call 'simputlc' to produce a light curve and assign it 
    // to the source.
    if (strlen(par.LCFile)>0) {
      sprintf(command, 
	      "simputlc Simput=%s LCFile=%s MJDREF=%.10lf chatter=%d history=%s",
	      par.Simput, par.LCFile, par.MJDREF, par.chatter, shistory);
      status=system(command);
      CHECK_STATUS_BREAK(status);
    }

    // Call 'simputpsd' to produce a power spectrum and assign
    // it to the source.
    if ((strlen(par.PSDFile)>0) || (par.LFQ!=0) || (par.HBOQ!=0) || 
	(par.Q1Q!=0) || (par.Q2Q!=0) || (par.Q3Q!=0)) {
      sprintf(command, 
	      "simputpsd Simput=%s "
	      "PSDnpt=%ld PSDfmin=%e PSDfmax=%e "
	      "LFQ=%e LFrms=%e "
	      "HBOf=%e HBOQ=%e HBOrms=%e "
	      "Q1f=%e Q1Q=%e Q1rms=%e "
	      "Q2f=%e Q2Q=%e Q2rms=%e "
	      "Q3f=%e Q3Q=%e Q3rms=%e "
	      "PSDFile=%s chatter=%d history=%s",
	      par.Simput, 
	      par.PSDnpt, par.PSDfmin, par.PSDfmax,
	      par.LFQ, par.LFrms,
	      par.HBOf, par.HBOQ, par.HBOrms,
	      par.Q1f, par.Q1Q, par.Q1rms,
	      par.Q2f, par.Q2Q, par.Q2rms,
	      par.Q3f, par.Q3Q, par.Q3rms,
	      par.PSDFile, par.chatter, shistory);
      status=system(command);
      CHECK_STATUS_BREAK(status);
    }

    // Call 'simputimg' to assign an image to the source.
    if (strlen(par.ImageFile)>0) {
      sprintf(command, 
	      "simputimg Simput=%s "
	      "ImageFile=%s chatter=%d history=%s",
	      par.Simput, 
	      par.ImageFile, par.chatter, shistory);
      status=system(command);
      CHECK_STATUS_BREAK(status);
    }

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int simputfile_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.
  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the output SIMPUT catalog file failed");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("Src_ID", &par->Src_ID);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the source ID failed");
    return(status);
  }

  status=ape_trad_query_string("Src_Name", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the source name failed");
    return(status);
  }
  strcpy(par->Src_Name, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the right ascension failed");
    return(status);
  }

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the declination failed");
    return(status);
  }

  status=ape_trad_query_float("srcFlux", &par->srcFlux);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the srcFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("plPhoIndex", &par->plPhoIndex);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the plPhoIndex parameter failed");
    return(status);
  }

  status=ape_trad_query_float("plFlux", &par->plFlux);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the plFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("bbkT", &par->bbkT);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the bbkT parameter failed");
    return(status);
  }

  status=ape_trad_query_float("bbFlux", &par->bbFlux);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the bbFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("flSigma", &par->flSigma);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the flSigma parameter failed");
    return(status);
  }

  status=ape_trad_query_float("flFlux", &par->flFlux);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the flFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("rflSpin", &par->rflSpin);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the rflSpin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("rflFlux", &par->rflFlux);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the rflFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("NH", &par->NH);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the N_H parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Emin", &par->Emin);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Emin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Emax", &par->Emax);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Emax parameter failed");
    return(status);
  }

  status=ape_trad_query_string("ISISFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the ISIS spectral parameter file failed");
    return(status);
  }
  strcpy(par->ISISFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XSPECFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the XSPEC spectral model file failed");
    return(status);
  }
  strcpy(par->XSPECFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("PHAFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the PHA file failed");
    return(status);
  }
  strcpy(par->PHAFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("LCFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the light curve file failed");
    return(status);
  }
  strcpy(par->LCFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter MJDREF failed");
    return(status);
  }

  status=ape_trad_query_long("PSDnpt", &par->PSDnpt);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter PSDnpt failed");
    return(status);
  }

  status=ape_trad_query_float("PSDfmin", &par->PSDfmin);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter PSDfmin failed");
    return(status);
  }

  status=ape_trad_query_float("PSDfmax", &par->PSDfmax);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter PSDfmax failed");
    return(status);
  }

  status=ape_trad_query_float("LFQ", &par->LFQ);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter LFQ failed");
    return(status);
  }

  status=ape_trad_query_float("LFrms", &par->LFrms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter LFrms failed");
    return(status);
  }

  status=ape_trad_query_float("HBOf", &par->HBOf);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter HBOf failed");
    return(status);
  }

  status=ape_trad_query_float("HBOQ", &par->HBOQ);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter HBOQ failed");
    return(status);
  }

  status=ape_trad_query_float("HBOrms", &par->HBOrms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter HBOrms failed");
    return(status);
  }

  status=ape_trad_query_float("Q1f", &par->Q1f);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q1f failed");
    return(status);
  }

  status=ape_trad_query_float("Q1Q", &par->Q1Q);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q1Q failed");
    return(status);
  }

  status=ape_trad_query_float("Q1rms", &par->Q1rms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q1rms failed");
    return(status);
  }

  status=ape_trad_query_float("Q2f", &par->Q2f);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q2f failed");
    return(status);
  }

  status=ape_trad_query_float("Q2Q", &par->Q2Q);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q2Q failed");
    return(status);
  }

  status=ape_trad_query_float("Q2rms", &par->Q2rms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q2rms failed");
    return(status);
  }

  status=ape_trad_query_float("Q3f", &par->Q3f);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q3f failed");
    return(status);
  }

  status=ape_trad_query_float("Q3Q", &par->Q3Q);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q3Q failed");
    return(status);
  }

  status=ape_trad_query_float("Q3rms", &par->Q3rms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the parameter Q3rms failed");
    return(status);
  }

  status=ape_trad_query_string("PSDFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the PSD file failed");
    return(status);
  }
  strcpy(par->PSDFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("ImageFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the image file failed");
    return(status);
  }
  strcpy(par->ImageFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("chatter", &par->chatter);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the chatter parameter failed");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the clobber parameter failed");
    return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the history parameter failed");
    return(status);
  }

  return(status);
}

