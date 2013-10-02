#include "simputsrc.h"


int simputsrc_main()
{
  // Program parameters.
  struct Parameters par;

  // SIMPUT data structures.
  SimputSrc* src=NULL;
  SimputCtlg* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputsrc");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputsrc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Check if the SIMPUT file already exists and remove the old 
    // one if necessary.
    int exists;
    fits_file_exists(par.Simput, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
        // Delete the file.
        remove(par.Simput);
      } else {
        // Throw an error.
        char msg[SIMPUT_MAXSTR];
        sprintf(msg, "file '%s' already exists", par.Simput);
        SIMPUT_ERROR(msg);
        status=EXIT_FAILURE;
        break;
      }
    }

    // Create a new SIMPUT catalog.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Insert a point-like source.
    if ((0==strcmp(par.Src_Name, "none"))||
	(0==strcmp(par.Src_Name, "NONE"))) {
      strcpy(par.Src_Name, "");
    }

    // Get a new source entry.
    src=newSimputSrcV(1, par.Src_Name, par.RA*M_PI/180., par.Dec*M_PI/180.,
		      0., 1., par.Emin, par.Emax, par.Flux,
		      "NULL", "NULL", "NULL", &status);
    CHECK_STATUS_BREAK(status);
    appendSimputSrc(cat, src, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Release memory.
  freeSimputSrc(&src);
  freeSimputCtlg(&cat, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int simputsrc_getpar(struct Parameters* const par)
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

  status=ape_trad_query_float("Flux", &par->Flux);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Flux parameter failed");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the clobber parameter failed");
    return(status);
  }

  return(status);
}

