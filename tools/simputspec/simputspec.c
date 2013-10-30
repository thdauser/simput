#include "simputspec.h"

// TODO Re-structure the code and use the tmpfile command.

int simputspec_main() 
{
  // Program parameters.
  struct Parameters par;

  fitsfile* fitsfile=NULL;

  // Temporary files for ISIS interaction.
  FILE* isiscmdfile=NULL;
  char isiscmdfilename[SIMPUT_MAXSTR]="";

  // XSPEC iplot file containing a spectrum.
  FILE* xspecfile=NULL;

  // Instrument response.
  struct ARF* arf=NULL;
  struct RMF* rmf=NULL;

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
  set_toolversion("0.07");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputspec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Check the input type for the spectrum.
    // Check the specification of an ISIS parameter file, an
    // XSPEC file, a PHA file, and the individual spectral components.
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
    if ((0==strcmp(par.PHAFile, "none"))||
	(0==strcmp(par.PHAFile, "NONE"))) {
      strcpy(par.PHAFile, "");
    }

    int noptions=0;
    if ((par.plFlux>0.) || (par.bbFlux>0.) || 
        (par.flFlux>0.) || (par.rflFlux>0.)) {
      use_components=1;
      noptions++;
    }
    if (strlen(par.ISISFile)>0) {
      noptions++;
    }
    if (strlen(par.XSPECFile)>0) {
      noptions++;
    }
    if (strlen(par.PHAFile)>0) {
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
        fits_open_table(&fitsfile, filename, READONLY, &status);
        CHECK_STATUS_BREAK(status);

        // Load the data from the file.
        int anynull;
        if (0==ii) {
          // Determine the number of rows.
          fits_get_num_rows(fitsfile, &nrows, &status);
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
          fits_read_col(fitsfile, TFLOAT, 1, 1, 1, nrows, 0, simputspec->energy,
			&anynull, &status);
          CHECK_STATUS_BREAK(status);
	  long jj;
	  for (jj=0; jj<nrows; jj++) {
	    simputspecbuffer->energy[jj]=simputspec->energy[jj];
	  }

        } else {
          // Check whether the number of entries is
          // consistent with the previous files.
          fits_get_num_rows(fitsfile, &nrows, &status);
          CHECK_STATUS_BREAK(status);

          if (nrows!=simputspec->nentries) {
            SIMPUT_ERROR("inconsistent sizes of spectra");
            status=EXIT_FAILURE;
            break;
          }
        }

        // Read the flux column.
        fits_read_col(fitsfile, TFLOAT, 2, 1, 1, nrows, 0, 
		      simputspecbuffer->pflux, &anynull, &status);
        CHECK_STATUS_BREAK(status);

        fits_close_file(fitsfile, &status);
        fitsfile=NULL;
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

    } else if (strlen(par.XSPECFile)>0) {
      // The spectrum is contained in an ASCII file produced by XSPEC,
      // and has to be loaded from there.

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

    } else {
      // The spectrum has to be obtained from a PHA file.

      // Load the spectrum from the PHA file.
      long nrows;
      fits_open_table(&fitsfile, par.PHAFile, READONLY, &status);
      CHECK_STATUS_BREAK(status);
      fits_get_num_rows(fitsfile, &nrows, &status);
      CHECK_STATUS_BREAK(status);

      // Allocate memory.
      simputspec->nentries=nrows;
      simputspec->energy=(float*)malloc(nrows*sizeof(float));
      CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");
      simputspec->pflux=(float*)malloc(nrows*sizeof(float));
      CHECK_NULL_BREAK(simputspec->pflux, status, "memory allocation failed");

      // Need to distinguish whether the file contains counts or rate.
      char comment[SIMPUT_MAXSTR];
      char hduclas3[SIMPUT_MAXSTR];
      fits_read_key(fitsfile, TSTRING, "HDUCLAS3", hduclas3, comment, &status);
      if (EXIT_SUCCESS!=status) {
	SIMPUT_ERROR("could not find keyword 'HDUCLAS3' in PHA file");
	break;
      }

      if ((0==strcmp(hduclas3, "COUNT"))||(0==strcmp(hduclas3, "count"))) {
	float exposure;
	fits_read_key(fitsfile, TFLOAT, "EXPOSURE", &exposure, comment, &status);
	if (EXIT_SUCCESS!=status) {
	  SIMPUT_ERROR("could not find keyword 'EXPOSURE' in PHA file");
	  break;
	}

	int ccount;
	fits_get_colnum(fitsfile, CASEINSEN, "COUNTS", &ccount, &status);
	if (EXIT_SUCCESS!=status) {
	  SIMPUT_ERROR("could not find column 'COUNTS' in PHA file");
	  break;
	}

	int anynull=0;
	fits_read_col(fitsfile, TFLOAT, ccount, 1, 1, nrows, 0, simputspec->pflux,
		      &anynull, &status);

	// Divide by exposure time.
	long ii;
	for (ii=0; ii<nrows; ii++) {
	  simputspec->pflux[ii]*=1./exposure;
	}

      } else if ((0==strcmp(hduclas3, "RATE"))||(0==strcmp(hduclas3, "rate"))) {
	int crate;
	fits_get_colnum(fitsfile, CASEINSEN, "RATE", &crate, &status);
	if (EXIT_SUCCESS!=status) {
	  SIMPUT_ERROR("could not find column 'RATE' in PHA file");
	  break;
	}

	int anynull=0;
	fits_read_col(fitsfile, TFLOAT, crate, 1, 1, nrows, 0, simputspec->pflux,
		      &anynull, &status);

      } else {
	SIMPUT_ERROR("invalid value for keyword 'HDUCLAS3'");
	status=EXIT_FAILURE;
	break;
      }


      // Load the ARF and the RMF.
      char ancrfile[SIMPUT_MAXSTR];
      fits_read_key(fitsfile, TSTRING, "ANCRFILE", ancrfile, comment, &status);
      if (EXIT_SUCCESS!=status) {
	SIMPUT_ERROR("could not find keyword 'ANCRFILE' in event file");
	break;
      }
      char respfile[SIMPUT_MAXSTR];
      fits_read_key(fitsfile, TSTRING, "RESPFILE", respfile, comment, &status);
      if (EXIT_SUCCESS!=status) {
	SIMPUT_ERROR("could not find keyword 'RESPFILE' in event file");
	break;
      }
      arf=loadARF(ancrfile, &status);
      CHECK_STATUS_BREAK(status);
      rmf=loadRMF(respfile, &status);
      CHECK_STATUS_BREAK(status);
      loadEbounds(rmf, respfile, &status);
      CHECK_STATUS_BREAK(status);

      // Check that RMF and ARF have the same number of channels.
      if (rmf->NumberEnergyBins!=arf->NumberEnergyBins) {
	SIMPUT_ERROR("ARF and RMF must contain the same number of energy bins");
	status=EXIT_FAILURE;
	break;
      }

      // Deconvolve the data according to the method presented by Nowak (2005).
      long ii;
      for (ii=0; ii<simputspec->nentries; ii++) {
	// Store the energy.
	float lo, hi;
	getEBOUNDSEnergyLoHi(ii, rmf, &lo, &hi, &status);
	CHECK_STATUS_BREAK(status);
	simputspec->energy[ii]=0.5*(lo+hi);

	// Calculate the integral \int R(h,E)A(E)dE.
	float area=0.;
	long kk;
	for (kk=0; kk<arf->NumberEnergyBins; kk++) {
	  area+=ReturnRMFElement(rmf, ii, kk)*arf->EffArea[kk];
	}
	
	// Divide by the area and the width of the energy bin.
	simputspec->pflux[ii]*=1./area/(hi-lo);	
      }
      CHECK_STATUS_BREAK(status);

    }

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
    saveSimputMIdpSpec(simputspec, par.Simput, par.Extname, par.Extver, &status);
    CHECK_STATUS_BREAK(status);


    // Open the SimputCtlg.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Set the spectrum reference in the source catalog.
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
    char* specref=(char*)malloc(32*sizeof(char));
    CHECK_NULL_BREAK(specref, status, "memory allocation failed");
    sprintf(specref, "[%s,%d]", par.Extname, par.Extver);
    fits_write_col(cat->fptr, TSTRING, cat->cspectrum, 1, 1, 1,
		   &specref, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the source flux in the catalog is set. If not (value=0.0),
    // set the flux according to the spectrum.
    int anynul=0;
    float srcflux=0.0;
    fits_read_col(cat->fptr, TFLOAT, cat->cflux, 1, 1, 1,
		  &srcflux, &srcflux, &anynul, &status);
    CHECK_STATUS_BREAK(status);

    if (0.0==srcflux) {
      // Determine the reference band.
      float srce_min=0.0, srce_max=0.0;
      fits_read_col(cat->fptr, TFLOAT, cat->ce_min, 1, 1, 1,
		    &srce_min, &srce_min, &anynul, &status);
      fits_read_col(cat->fptr, TFLOAT, cat->ce_max, 1, 1, 1,
		    &srce_max, &srce_max, &anynul, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the flux in the reference band.
      srcflux=getSimputMIdpSpecBandFlux(simputspec, srce_min, srce_max);

      // Store the flux in the source catalog.
      fits_write_col(cat->fptr, TFLOAT, cat->cflux, 1, 1, 1,
		     &srcflux, &status);
      CHECK_STATUS_BREAK(status);
    }

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
  if (NULL!=fitsfile) {
    fits_close_file(fitsfile, &status);
    fitsfile=NULL;
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
  freeRMF(rmf);
  freeARF(arf);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int simputspec_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

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

  status=ape_trad_query_string("PHAFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the PHA file failed");
    return(status);
  }
  strcpy(par->PHAFile, sbuffer);
  free(sbuffer);

  return(status);
}

