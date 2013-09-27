#include "simputspec.h"

// TODO Re-structure the code and use the tmpfile command.

int simputspec_main() 
{
  // Program parameters.
  struct Parameters par;

  // Temporary files for ISIS interaction.
  FILE* isiscmdfile=NULL;
  char isiscmdfilename[SIMPUT_MAXSTR]="";
  fitsfile* isisfitsfile=NULL;

  // XSPEC iplot file containing a spectrum.
  FILE* xspecfile=NULL;

  // Flag, whether the spectrum should be constructed from 
  // different components.
  int use_components=0;

  // Output SimputMIdpSpec.
  SimputMIdpSpec* simputspec=NULL;
  SimputMIdpSpec* simputspecbuffer=NULL;

  // SIMPUT catalog the spectrum should be attached to.
  SimputCtlg* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputspec");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputspec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Check the input type for the spectrum.
    // Check the specification of an ISIS parameter file, an
    // XSPEC file, and the individual spectral components.
    // Only one of these 3 option may be used. In case multiple of
    // them exist, throw an error message and abort.

    if ((0==strcmp(par.ISISFile, "none"))||
	(0==strcmp(par.ISISFile, "NONE"))) {
      strcpy(par.ISISFile, "");
    }

    if ((0==strcmp(par.XSPECFile, "none"))||
	(0==strcmp(par.XSPECFile, "NONE"))) {
      strcpy(par.XSPECFile, "");
    }

    int noptions=0;
    if (strlen(par.ISISFile)>0) {
      noptions++;
    }
    if (strlen(par.XSPECFile)>0) {
      noptions++;
    }
    if ((par.plFlux>0.) || (par.bbFlux>0.) || 
        (par.flFlux>0.) || (par.rflFlux>0.)) {
      use_components=1;
      noptions++;
    }
    if (0==noptions) {
      SIMPUT_ERROR("no spectral model specified");
      status=EXIT_FAILURE;
      break;
    }
    if (noptions>1) {
      SIMPUT_ERROR("specification of multiple spectral models");
      status=EXIT_FAILURE;
      break;
    }
    // END of checking the input type for the spectrum.

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Create the spectrum.

    // If individual components or an ISIS .par file are given,
    // we have to run ISIS in order to produce a spectrum.
    if ((strlen(par.ISISFile)>0) || (use_components>0)) {

      // Open the ISIS command file.
      sprintf(isiscmdfilename, "%s.isis", par.Simput);
      isiscmdfile=fopen(isiscmdfilename,"w");
      CHECK_NULL_BREAK(isiscmdfile, status, "opening temporary file failed");

      // Write the header.
      fprintf(isiscmdfile, "require(\"isisscripts\");\n");
      fprintf(isiscmdfile, "()=xspec_abund(\"wilm\");\n");
      fprintf(isiscmdfile, "use_localmodel(\"relline\");\n");

      // Define the energy grid.
      fprintf(isiscmdfile, "variable lo=[0.05:100.0:0.01];\n");
      fprintf(isiscmdfile, "variable hi=make_hi_grid(lo);\n");
      fprintf(isiscmdfile, "variable flux;\n");
      fprintf(isiscmdfile, "variable spec;\n");

      // Distinguish whether the individual spectral components or
      // an ISIS spectral parameter file should be used.
      if (strlen(par.ISISFile)==0) {

        // Loop over the different components of the spectral model.
        int ii;
        for (ii=0; ii<4; ii++) {

          // Define the spectral model and set the parameters.
          switch(ii) {
          case 0:
            fprintf(isiscmdfile, "fit_fun(\"phabs(1)*powerlaw(1)\");\n");
            fprintf(isiscmdfile, "set_par(\"powerlaw(1).PhoIndex\", %e);\n",
		    par.plPhoIndex);
            break;
          case 1:
            fprintf(isiscmdfile, "fit_fun(\"phabs(1)*bbody(1)\");\n");
            fprintf(isiscmdfile, "set_par(\"bbody(1).kT\", %e);\n", par.bbkT);
            break;
          case 2:
            fprintf(isiscmdfile, "fit_fun(\"phabs(1)*egauss(1)\");\n");
            fprintf(isiscmdfile, "set_par(\"egauss(1).center\", 6.4);\n");
            fprintf(isiscmdfile, "set_par(\"egauss(1).sigma\", %e);\n", par.flSigma);
            break;
          case 3:
            fprintf(isiscmdfile, "fit_fun(\"phabs(1)*relline(1)\");\n");
            fprintf(isiscmdfile, "set_par(\"relline(1).lineE\", 6.4);\n");
            fprintf(isiscmdfile, "set_par(\"relline(1).a\", %f);\n", par.rflSpin);
            break;
          default:
            status=EXIT_FAILURE;
            break;
          }
          CHECK_STATUS_BREAK(status);

          // Absorption is the same for all spectral components.
          fprintf(isiscmdfile, "set_par(\"phabs(1).nH\", %e);\n", par.NH);

          // Evaluate the spectral model and store the data in a temporary
          // FITS file.
          fprintf(isiscmdfile, "flux=eval_fun_keV(lo, hi)/(hi-lo);\n");
          fprintf(isiscmdfile, "spec=struct{ENERGY=0.5*(lo+hi), FLUX=flux};\n");
          fprintf(isiscmdfile,
		  "fits_write_binary_table(\"%s.spec%d\",\"SPECTRUM\", spec);\n",
		  par.Simput, ii);
        }
        CHECK_STATUS_BREAK(status);
        // END of loop over the different spectral components.

      } else { 
        // An ISIS parameter file with an explizit spectral
        // model is given.
        fprintf(isiscmdfile, "load_par(\"%s\");\n", par.ISISFile);
        fprintf(isiscmdfile, "flux=eval_fun_keV(lo, hi)/(hi-lo);\n");
        fprintf(isiscmdfile, "spec=struct{ENERGY=0.5*(lo+hi), FLUX=flux};\n");
        fprintf(isiscmdfile,
		"fits_write_binary_table(\"%s.spec0\",\"SPECTRUM\", spec);\n",
		par.Simput);
      }
      // END of using an explicit spectral model given in an ISIS 
      // parameter file.

      fprintf(isiscmdfile, "exit;\n");

      // End of writing the ISIS command file.
      fclose(isiscmdfile);
      isiscmdfile=NULL;

      // Construct the shell command to run ISIS.
      char command[SIMPUT_MAXSTR];
      strcpy(command, "isis ");
      strcat(command, isiscmdfilename);

      // Run ISIS.
      status=system(command);
      CHECK_STATUS_BREAK(status);

    } // END of running ISIS.

    // Add the spectra and store the total spectrum.
    simputspec=newSimputMIdpSpec(&status);
    CHECK_STATUS_BREAK(status);

    // Read the spectrum from the file(s).
    if ((strlen(par.ISISFile)>0) || (use_components>0)) {
      
      // Allocate memory for the buffer.
      simputspecbuffer=newSimputMIdpSpec(&status);
      CHECK_STATUS_BREAK(status);

      // Loop over the different components of the spectral model.
      long nrows=0;
      int ii;
      for (ii=0; ii<4; ii++) {

        // Determine the file name.
        char filename[SIMPUT_MAXSTR];
        sprintf(filename, "%s.spec%d", par.Simput, ii);
        fits_open_table(&isisfitsfile, filename, READONLY, &status);
        CHECK_STATUS_BREAK(status);

        // Load the data from the file.
        int anynull;
        if (0==ii) {
          // Determine the number of rows.
          fits_get_num_rows(isisfitsfile, &nrows, &status);
          CHECK_STATUS_BREAK(status);

          // Allocate memory.
          simputspec->nentries=nrows;
          simputspec->energy=(float*)malloc(nrows*sizeof(float));
          CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");
          simputspec->pflux=(float*)malloc(nrows*sizeof(float));
          CHECK_NULL_BREAK(simputspec->pflux, status, "memory allocation failed");

	  simputspecbuffer->nentries=nrows;
          simputspecbuffer->energy=(float*)malloc(nrows*sizeof(float));
          CHECK_NULL_BREAK(simputspecbuffer->energy, status, "memory allocation failed");
          simputspecbuffer->pflux=(float*)malloc(nrows*sizeof(float));
          CHECK_NULL_BREAK(simputspecbuffer->pflux, status, "memory allocation failed");

          // Read the energy column.
          fits_read_col(isisfitsfile, TFLOAT, 1, 1, 1, nrows, 0, simputspec->energy,
			&anynull, &status);
          CHECK_STATUS_BREAK(status);
	  long jj;
	  for (jj=0; jj<nrows; jj++) {
	    simputspecbuffer->energy[jj]=simputspec->energy[jj];
	  }

        } else {
          // Check whether the number of entries is
          // consistent with the previous files.
          fits_get_num_rows(isisfitsfile, &nrows, &status);
          CHECK_STATUS_BREAK(status);

          if (nrows!=simputspec->nentries) {
            SIMPUT_ERROR("inconsistent sizes of spectra");
            status=EXIT_FAILURE;
            break;
          }
        }

        // Read the flux column.
        fits_read_col(isisfitsfile, TFLOAT, 2, 1, 1, nrows, 0, 
		      simputspecbuffer->pflux, &anynull, &status);
        CHECK_STATUS_BREAK(status);

        fits_close_file(isisfitsfile, &status);
        isisfitsfile=NULL;
        CHECK_STATUS_BREAK(status);

        // If the spectrum is given via individual components, they
        // have to be normalized according to their respective fluxes.
        if (strlen(par.ISISFile)==0) {
          // Determine the required flux in the reference band.
          float shouldflux=0.;
          switch(ii) {
          case 0:
            shouldflux=par.plFlux;
            break;
          case 1:
            shouldflux=par.bbFlux;
            break;
          case 2:
            shouldflux=par.flFlux;
            break;
          case 3:
            shouldflux=par.rflFlux;
            break;
          default:
            status=EXIT_FAILURE;
            break;
          }
          CHECK_STATUS_BREAK(status);

          float factor;
          if (shouldflux==0.) {
            factor=0.;
          } else {
            // Determine the factor between the actual flux in the reference
            // band and the required flux.
            factor=shouldflux/
	      getSimputMIdpSpecBandFlux(simputspecbuffer, par.Emin, par.Emax);
          }

          // Add the normalized component to the total spectrum.
          if (factor>0.) {
            long jj;
            for (jj=0; jj<nrows; jj++) {
              simputspec->pflux[jj]+=simputspecbuffer->pflux[jj]*factor;
            }
          }

        } else {
          // The spectral model is given in an ISIS parameter file.
          // Therefore we do not have to normalize it, but can directly
          // add it to the SIMPUT spectrum.
          long jj;
          for (jj=0; jj<nrows; jj++) {
            simputspec->pflux[jj]=simputspecbuffer->pflux[jj];
          }

          // Since there are no further components, we can skip
          // the further processing of the loop.
          break;
        }
      }
      CHECK_STATUS_BREAK(status);
      // END of loop over the different spectral components.

    } else {
      // The XPSEC spectrum is contained in an ASCII file and has 
      // to be loaded from there.

      // Open the file.
      xspecfile=fopen(par.XSPECFile,"r");
      CHECK_NULL_BREAK(xspecfile, status, "could not open XSPEC file");

      // Determine the number of rows.
      long nlines=0;
      char c=0;
      while(!feof(xspecfile)) {
        c=fgetc(xspecfile);
        if ('\n'==c) {
          nlines++;
        }
      }
      // Check if the last line has been empty.
      if('\n'==c) {
        nlines--;
      }

      // The first 3 lines do not contain data.
      nlines-=3;

      // Allocate memory.
      simputspec->nentries=nlines;
      simputspec->energy=(float*)malloc(nlines*sizeof(float));
      CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");
      simputspec->pflux=(float*)malloc(nlines*sizeof(float));
      CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");

      // Reset the file pointer, read the data, and store them in
      // the SimputMIdpSpec data structure.
      rewind(xspecfile);
      // Read the first three lines.
      char sbuffer1[SIMPUT_MAXSTR], sbuffer2[SIMPUT_MAXSTR];
      int ibuffer;
      if (fscanf(xspecfile, "%s %s %d\n", sbuffer1, sbuffer2, &ibuffer)<3) {
	SIMPUT_ERROR("failed reading data from ASCII file");
	status=EXIT_FAILURE;
	break;
      }
      if (fscanf(xspecfile, "%s\n", sbuffer1)<1) {
	SIMPUT_ERROR("failed reading data from ASCII file");
	status=EXIT_FAILURE;
	break;
      }
      if (fscanf(xspecfile, "%s\n", sbuffer1)<1) {
	SIMPUT_ERROR("failed reading data from ASCII file");
	status=EXIT_FAILURE;
	break;
      }
      // Read the actual data.
      long ii;
      for (ii=0; ii<nlines; ii++) {
	float fbuffer;
        if (fscanf(xspecfile, "%f %f %f\n",
		   &(simputspec->energy[ii]), 
		   &fbuffer, 
		   &(simputspec->pflux[ii]))<3) {
	  SIMPUT_ERROR("failed reading data from ASCII file");
	  status=EXIT_FAILURE;
	  break;
	}
      }
      CHECK_STATUS_BREAK(status);

      // Close the file.
      fclose(xspecfile);
      xspecfile=NULL;

    } // END of loading the spectrum from an XSPEC file.

    long jj;
    for (jj=0; jj<simputspec->nentries; jj++) {
      // Check if the flux has a physically reasonable value.
      if ((simputspec->pflux[jj]<0.)||(simputspec->pflux[jj]>1.e12)) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "flux (%e photons/cm**2/keV) out of boundaries", 
		simputspec->pflux[jj]);
        SIMPUT_ERROR(msg);
        status=EXIT_FAILURE;
	break;
      }
    }
    CHECK_STATUS_BREAK(status);
    saveSimputMIdpSpec(simputspec, par.Simput, "SPECTRUM", 1, &status);
    CHECK_STATUS_BREAK(status);


    // Open the SimputCtlg.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Set the spectrum reference in the source catalog.
    char* specref=(char*)malloc(32*sizeof(char));
    CHECK_NULL_BREAK(specref, status, "memory allocation failed");
    strcpy(specref, "[SPECTRUM,1]");
    fits_write_col(cat->fptr, TSTRING, cat->cspectrum, 1, 1, 1,
		   &specref, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Close open files.
  if (NULL!=xspecfile) {
    fclose(xspecfile);
    xspecfile=NULL;
  }

  // Close the temporary files.
  if (NULL!=isiscmdfile) {
    fclose(isiscmdfile);
    isiscmdfile=NULL;
  }
  if (NULL!=isisfitsfile) {
    fits_close_file(isisfitsfile, &status);
    isisfitsfile=NULL;
  }
  // Remove the temporary files.
  if (strlen(isiscmdfilename)>0) {
    remove(isiscmdfilename);
  }
  if (use_components>0) {
    int ii;
    for (ii=0; ii<4; ii++) {
      char filename[SIMPUT_MAXSTR];
      sprintf(filename, "%s.spec%d", par.Simput, ii);
      remove(filename);
    }      
  }
  if (strlen(par.ISISFile)>0) {
    char filename[SIMPUT_MAXSTR];
    sprintf(filename, "%s.spec0", par.Simput);
    remove(filename);
  }

  // Release memory.
  freeSimputMIdpSpec(&simputspecbuffer);
  freeSimputMIdpSpec(&simputspec);
  freeSimputCtlg(&cat, &status);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int simputspec_getpar(struct Parameters* const par)
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

  return(status);
}

