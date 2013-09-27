#include "simputfile.h"


int simputfile_main() 
{
  // Program parameters.
  struct Parameters par;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputfile");
  set_toolversion("0.13");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputfile_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Check the input type for the power spectrum: individual 
    // components or and ASCII file. Only one of these two option 
    // may be used. In case multiple of them exist, throw an error 
    // message and abort.

    if ((0==strcmp(par.LCFile, "none"))||
	(0==strcmp(par.LCFile, "NONE"))) {
      strcpy(par.LCFile, "");
    }
    if ((0==strcmp(par.PSDFile, "none"))||
	(0==strcmp(par.PSDFile, "NONE"))) {
      strcpy(par.PSDFile, "");
    }

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
	    "simputsrc Simput=%s Src_Name=%s RA=%f Dec=%f "
	    "Emin=%f Emax=%f Flux=%e chatter=%d clobber=%s history=%s",
	    par.Simput, par.Src_Name, par.RA, par.Dec,
	    par.Emin, par.Emax, 
	    par.plFlux+par.bbFlux+par.flFlux+par.rflFlux,
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
	    "Emin=%f Emax=%f ISISFile=%s XSPECFile=%s "
	    "chatter=%d history=%s",
	    par.Simput,
	    par.plPhoIndex, par.plFlux,
	    par.bbkT, par.bbFlux,
	    par.flSigma, par.flFlux,
	    par.rflSpin, par.rflFlux,
	    par.NH,
	    par.Emin, par.Emax, par.ISISFile, par.XSPECFile,
	    par.chatter, shistory);
    status=system(command);
    CHECK_STATUS_BREAK(status);

    // Call 'simputlc' to produce a light curve and assign it 
    // to the source.
    if (strlen(par.LCFile)>0) {
      sprintf(command, 
	      "simputlc Simput=%s LCFile=%s chatter=%d history=%s",
	      par.Simput, par.LCFile, par.chatter, shistory);
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

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
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
    SIMPUT_ERROR("reading the name of the XSPEC spectrum file failed");
    return(status);
  }
  strcpy(par->XSPECFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("LCFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the light curve file failed");
    return(status);
  }
  strcpy(par->LCFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_long("PSDnpt", &par->PSDnpt);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the PSDnpt parameter failed");
    return(status);
  }

  status=ape_trad_query_float("PSDfmin", &par->PSDfmin);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the PSDfmin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("PSDfmax", &par->PSDfmax);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the PSDfmax parameter failed");
    return(status);
  }

  status=ape_trad_query_float("LFQ", &par->LFQ);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the LFQ parameter failed");
    return(status);
  }

  status=ape_trad_query_float("LFrms", &par->LFrms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the LFrms parameter failed");
    return(status);
  }

  status=ape_trad_query_float("HBOf", &par->HBOf);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the HBOf parameter failed");
    return(status);
  }

  status=ape_trad_query_float("HBOQ", &par->HBOQ);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the HBOQ parameter failed");
    return(status);
  }

  status=ape_trad_query_float("HBOrms", &par->HBOrms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the HBOrms parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q1f", &par->Q1f);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q1f parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q1Q", &par->Q1Q);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q1Q parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q1rms", &par->Q1rms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q1rms parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q2f", &par->Q2f);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q2f parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q2Q", &par->Q2Q);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q2Q parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q2rms", &par->Q2rms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q2rms parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q3f", &par->Q3f);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q3f parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q3Q", &par->Q3Q);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q3Q parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Q3rms", &par->Q3rms);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Q3rms parameter failed");
    return(status);
  }

  status=ape_trad_query_string("PSDFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the PSD file failed");
    return(status);
  }
  strcpy(par->PSDFile, sbuffer);
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

