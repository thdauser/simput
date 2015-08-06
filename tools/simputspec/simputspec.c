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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "simputspec.h"


int simputspec_main() 
{
  // Program parameters.
  struct Parameters par;

  fitsfile* fptr=NULL;

  // Temporary files for ISIS and Xspec interaction.
  FILE* cmdfile=NULL;
  char cmdfilename[L_tmpnam]="";

  // XSPEC .qdp file containing the spectrum.
  FILE* xspecfile=NULL;

  // ASCII file containing the spectrum.
  FILE* asciifile=NULL;

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
  set_toolversion("0.11");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputspec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Check if the specified energy ranges are reasonable.
    if (par.Elow>par.Emin) {
      SIMPUT_ERROR("parameter 'Emin' must be higher than 'Elow'");
      status=EXIT_FAILURE;
      break;
    }
    if (par.Eup<par.Emax) {
      SIMPUT_ERROR("parameter 'Emax' may not exceed 'Eup'");
      status=EXIT_FAILURE;
      break;
    }
    if (par.Estep>par.Eup-par.Elow) {
      SIMPUT_ERROR("parameter 'Estep' may not exceed difference "
		   "between 'Eup' and 'Elow'");
      status=EXIT_FAILURE;
      break;
    }

    // Check the input type for the spectrum.
    // Check the specification of an ISIS parameter file, an
    // XSPEC file, a PHA file, and the individual spectral components.
    // Only one of these 3 option may be used. In case multiple of
    // them exist, throw an error message and abort.

    if ((0==strcmp(par.ISISFile, "none"))||
	(0==strcmp(par.ISISFile, "NONE"))) {
      strcpy(par.ISISFile, "");
    }
    if ((0==strcmp(par.ISISPrep, "none"))||
	(0==strcmp(par.ISISPrep, "NONE"))) {
      strcpy(par.ISISPrep, "");
    }
    if ((0==strcmp(par.XSPECFile, "none"))||
	(0==strcmp(par.XSPECFile, "NONE"))) {
      strcpy(par.XSPECFile, "");
    }
    if ((0==strcmp(par.PHAFile, "none"))||
	(0==strcmp(par.PHAFile, "NONE"))) {
      strcpy(par.PHAFile, "");
    }
    if ((0==strcmp(par.ASCIIFile, "none"))||
	(0==strcmp(par.ASCIIFile, "NONE"))) {
      strcpy(par.ASCIIFile, "");
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
    if (strlen(par.ASCIIFile)>0) {
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

    	write_isisSpec_fits_file(par.Simput, par.ISISFile, par.ISISPrep,
    			par.ISISPostCmd, par.Elow, par.Eup, par.Estep,
				par.plPhoIndex, par.bbkT, par.flSigma, par.rflSpin, par.NH,
				&status);

    } // END of running ISIS.


    // If an Xspec .xcm file is given, we have to run Xspec in order 
    // to produce the spectrum.
    if (strlen(par.XSPECFile)>0) {

    	write_xspecSpec_file(par.Simput, par.XSPECFile, par.XSPECPostCmd,
    			par.Elow, par.Eup, par.Estep, &status);

    } // END of running Xspec.

    // Read the spectrum from the temporary file and store
    // it in the SIMPUT file.
    simputspec=newSimputMIdpSpec(&status);
    CHECK_STATUS_BREAK(status);

    if ((strlen(par.ISISFile)>0) || (use_components>0)) {
      
    	read_isisSpec_fits_file(par.Simput, simputspec,
    			par.ISISFile, par.Emin, par.Emax,
				par.plFlux, par.bbFlux, par.flFlux, par.rflFlux,
				&status);

    } else if (strlen(par.XSPECFile)>0) {
      // The spectrum is contained in a .qdp file produced by XSPEC/PLT,
      // and has to be loaded from there.

    	read_xspecSpec_file(par.Simput, simputspec, &status);

    } else if (strlen(par.ASCIIFile)>0) {
        // The spectrum is contained in a ascii file,
        // and has to be loaded from there.

        // Open the file.
        asciifile=fopen(par.ASCIIFile, "r");
        CHECK_NULL_BREAK(asciifile, status, "could not open ascii file");

        // Determine the number of rows.
        long nlines=0;
        char c=0;
        while(!feof(asciifile)) {
          c=fgetc(asciifile);
          if ('\n'==c) {
            nlines++;
          }
        }
        // Check if the last line has been empty.
        if('\n'==c) {
          nlines--;
        }

        // Allocate memory.
        simputspec->nentries=nlines;
        simputspec->energy=(float*)malloc(nlines*sizeof(float));
        CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");
        simputspec->fluxdensity=(float*)malloc(nlines*sizeof(float));
        CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");

        // Reset the file pointer, read the data, and store them in
        // the SimputMIdpSpec data structure.
        rewind(asciifile);
        // Read the actual data.
        long ii;
        for (ii=0; ii<nlines; ii++) {
  	char linebuffer[SIMPUT_MAXSTR];
  	if(fgets(linebuffer, SIMPUT_MAXSTR, asciifile)!=NULL){
  	  if(sscanf(linebuffer, "%f %f",
  		 &(simputspec->energy[ii]),
  		 &(simputspec->fluxdensity[ii]))!=2) {
  	  SIMPUT_ERROR("failed reading data from ASCII file");
  	  status=EXIT_FAILURE;
  	  break;
  	  }
  	}
        }
        CHECK_STATUS_BREAK(status);

        // Close the file.
        fclose(asciifile);
        asciifile=NULL;

    } else {
      // The spectrum has to be obtained from a PHA file.

      // Load the spectrum from the PHA file.
      long nrows;
      fits_open_table(&fptr, par.PHAFile, READONLY, &status);
      CHECK_STATUS_BREAK(status);
      fits_get_num_rows(fptr, &nrows, &status);
      CHECK_STATUS_BREAK(status);

      // Allocate memory.
      simputspec->nentries=nrows;
      simputspec->energy=(float*)malloc(nrows*sizeof(float));
      CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");
      simputspec->fluxdensity=(float*)malloc(nrows*sizeof(float));
      CHECK_NULL_BREAK(simputspec->fluxdensity, status, "memory allocation failed");

      // Need to distinguish whether the file contains counts or rate.
      char comment[SIMPUT_MAXSTR];
      char hduclas3[SIMPUT_MAXSTR];
      fits_read_key(fptr, TSTRING, "HDUCLAS3", hduclas3, comment, &status);
      if (EXIT_SUCCESS!=status) {
	SIMPUT_ERROR("could not find keyword 'HDUCLAS3' in PHA file");
	break;
      }

      if ((0==strcmp(hduclas3, "COUNT"))||(0==strcmp(hduclas3, "count"))) {
	float exposure;
	fits_read_key(fptr, TFLOAT, "EXPOSURE", &exposure, comment, &status);
	if (EXIT_SUCCESS!=status) {
	  SIMPUT_ERROR("could not find keyword 'EXPOSURE' in PHA file");
	  break;
	}

	int ccount;
	fits_get_colnum(fptr, CASEINSEN, "COUNTS", &ccount, &status);
	if (EXIT_SUCCESS!=status) {
	  SIMPUT_ERROR("could not find column 'COUNTS' in PHA file");
	  break;
	}

	int anynull=0;
	fits_read_col(fptr, TFLOAT, ccount, 1, 1, nrows, 0,
		      simputspec->fluxdensity, &anynull, &status);

	// Divide by exposure time.
	long ii;
	for (ii=0; ii<nrows; ii++) {
	  simputspec->fluxdensity[ii]*=1./exposure;
	}

      } else if ((0==strcmp(hduclas3, "RATE"))||(0==strcmp(hduclas3, "rate"))) {
	int crate;
	fits_get_colnum(fptr, CASEINSEN, "RATE", &crate, &status);
	if (EXIT_SUCCESS!=status) {
	  SIMPUT_ERROR("could not find column 'RATE' in PHA file");
	  break;
	}

	int anynull=0;
	fits_read_col(fptr, TFLOAT, crate, 1, 1, nrows, 0, 
		      simputspec->fluxdensity, &anynull, &status);

      } else {
	SIMPUT_ERROR("invalid value for keyword 'HDUCLAS3'");
	status=EXIT_FAILURE;
	break;
      }


      // Load the ARF and the RMF.
      char ancrfile[SIMPUT_MAXSTR];
      fits_read_key(fptr, TSTRING, "ANCRFILE", ancrfile, comment, &status);
      if (EXIT_SUCCESS!=status) {
	SIMPUT_ERROR("could not find keyword 'ANCRFILE' in event file");
	break;
      }
      char respfile[SIMPUT_MAXSTR];
      fits_read_key(fptr, TSTRING, "RESPFILE", respfile, comment, &status);
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
	simputspec->fluxdensity[ii]*=1./area/(hi-lo);	
      }
      CHECK_STATUS_BREAK(status);

    }

    long jj;
    for (jj=0; jj<simputspec->nentries; jj++) {
      // Check if the flux has a physically reasonable value.
      if ((simputspec->fluxdensity[jj]<0.)||(simputspec->fluxdensity[jj]>1.e12)) {
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "flux (%e photons/cm**2/keV) out of limits", 
		simputspec->fluxdensity[jj]);
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
  if (NULL!=cmdfile) {
    fclose(cmdfile);
    cmdfile=NULL;
  }
  if (NULL!=fptr) {
    fits_close_file(fptr, &status);
    fptr=NULL;
  }
  // Remove the temporary files.
  if (strlen(cmdfilename)>0) {
    remove(cmdfilename);
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
  if (strlen(par.XSPECFile)>0) {
    char filename[SIMPUT_MAXSTR];
    sprintf(filename, "%s.qdp", par.Simput);
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

  status=ape_trad_query_float("Elow", &par->Elow);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Elow parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Eup", &par->Eup);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Eup parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Estep", &par->Estep);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the Estep parameter failed");
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

  status=ape_trad_query_string("ISISPrep", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the ISIS prep file failed");
    return(status);
  }
  strcpy(par->ISISPrep, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("ISISPostCmd", &sbuffer);
    if (EXIT_SUCCESS!=status) {
      SIMPUT_ERROR("reading the name of the ISIS post cmd file failed");
      return(status);
    }
    strcpy(par->ISISPostCmd, sbuffer);
    free(sbuffer);



  status=ape_trad_query_string("XSPECFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the XSPEC spectrum file failed");
    return(status);
  }
  strcpy(par->XSPECFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XSPECPostCmd", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the XSPEC spectrum file failed");
    return(status);
  }
  strcpy(par->XSPECPostCmd, sbuffer);
  free(sbuffer);


  status=ape_trad_query_string("ASCIIFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the name of the ASCII spectrum file failed");
    return(status);
  }
  strcpy(par->ASCIIFile, sbuffer);
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

