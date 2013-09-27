#include "simputpsd.h"


int simputpsd_main() 
{
  // Program parameters.
  struct Parameters par;

  // Input ASCII file containing the PSD.
  FILE* asciipsd=NULL;

  // Output SimputPSD.
  SimputPSD* simputpsd=NULL;

  // SimputCtlg the PSD should be attached to.
  SimputCtlg* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputpsd");
  set_toolversion("0.02");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputpsd_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Check the input type for the power spectrum: individual 
    // components or and ASCII file. Only one of these two option 
    // may be used. In case multiple of them exist, throw an error 
    // message and abort.

    if ((0==strcmp(par.PSDFile, "none"))||
	(0==strcmp(par.PSDFile, "NONE"))) {
      strcpy(par.PSDFile, "");
    }

    int noptions=0;
    if (strlen(par.PSDFile)>0) {
      noptions++;
    }
    if ((par.LFQ!=0) || (par.HBOQ!=0) || 
	(par.Q1Q!=0) || (par.Q2Q!=0) || (par.Q3Q!=0)) {
      noptions++;
    }
    if (0==noptions) {
      SIMPUT_ERROR("no PSD model specified");
      status=EXIT_FAILURE;
      break;
    }
    if (noptions>1) {
      SIMPUT_ERROR("specification of multiple PSD models");
      status=EXIT_FAILURE;
      break;
    }
    // END of checking the input type for the power spectrum.

    // ---- END of Initialization ----


    // ---- Main Part ----

    simputpsd=newSimputPSD(&status);
    CHECK_STATUS_BREAK(status);

    if (strlen(par.PSDFile)>0) {
    
      // Open the ASCII file with the PSD.
      asciipsd=fopen(par.PSDFile,"r");
      CHECK_NULL_BREAK(asciipsd, status, "could not open input PSD file");

      // Determine the number of rows.
      long nlines=0;
      char c=0;
      while(!feof(asciipsd)) {
	c=fgetc(asciipsd);
	if ('\n'==c) {
	  nlines++;
	}
      }
      // Check if the last line has been empty.
      if('\n'==c) {
	nlines--;
      }

      // Allocate memory.
      simputpsd->nentries=nlines;
      simputpsd->frequency=(float*)malloc(nlines*sizeof(float));
      CHECK_NULL_BREAK(simputpsd->frequency, status, "memory allocation failed");
      simputpsd->power=(float*)malloc(nlines*sizeof(float));
      CHECK_NULL_BREAK(simputpsd->power, status, "memory allocation failed");

      // Reset the file pointer, read the data and store them in
      // the SimputPSD data structure.
      rewind(asciipsd);
      long ii;
      for (ii=0; ii<nlines; ii++) {
	if (fscanf(asciipsd, "%f %f\n",
		   &(simputpsd->frequency[ii]), 
		   &(simputpsd->power[ii]))<2) {
	  SIMPUT_ERROR("failed reading PSD from ASCII file");
	  status=EXIT_FAILURE;
	  break;
	}
      }
      CHECK_STATUS_BREAK(status);

    } else { // Assemble the PSD from individual components.

      // Allocate memory.
      simputpsd->nentries=par.PSDnpt;
      simputpsd->frequency=(float*)malloc(par.PSDnpt*sizeof(float));
      CHECK_NULL_BREAK(simputpsd->frequency, status, "memory allocation failed");
      simputpsd->power=(float*)malloc(par.PSDnpt*sizeof(float));
      CHECK_NULL_BREAK(simputpsd->power, status, "memory allocation failed");

      // Generate log-scaled frequency grid.
      long ii;
      for(ii=0; ii<par.PSDnpt; ii++) {
        simputpsd->frequency[ii]=
	  exp(log(par.PSDfmin)+ii*(log(par.PSDfmax/par.PSDfmin)/par.PSDnpt));
      }

      // Calculate Lorentzians using Formula (5.1) in Pottschmidt, K.: 
      // Accretion Disk Weather of Black Hole X-Ray Binaries (2002), p. 95.
      float* Lzero=NULL;
      float* LHBO =NULL;
      float* LQ1  =NULL;
      float* LQ2  =NULL;
      float* LQ3  =NULL;

      // Zero order Lorentzian.
      if(par.LFQ!=0) {
        float zNorm=par.LFrms/sqrt(0.5-(atan(par.LFQ*(-1))/M_PI));
        Lzero=(float*)malloc(par.PSDnpt*sizeof(float));
	CHECK_NULL_BREAK(Lzero, status, "memory allocation failed");
        for(ii=0; ii<par.PSDnpt; ii++) {
          Lzero[ii]=
	    (1.0/M_PI)*
	    ((pow(zNorm, 2)*par.LFQ*1e-5)/
	     (pow(1e-5, 2)+(pow(par.LFQ, 2)*pow((simputpsd->frequency[ii]-1e-5), 2))));
          simputpsd->power[ii]+=Lzero[ii];
        }
      }

      // HBO Lorentzian.
      if(par.HBOf!=0) {
        float HBONorm=par.HBOrms/sqrt(0.5-(atan(par.HBOQ*(-1))/M_PI));
        LHBO=(float*)malloc(par.PSDnpt*sizeof(float));
	CHECK_NULL_BREAK(LHBO, status, "memory allocation failed");
        for(ii=0; ii<par.PSDnpt; ii++) {
          LHBO[ii]=
	    (1.0/M_PI)*
	    ((pow(HBONorm, 2)*par.HBOQ*par.HBOf)/
	     (pow(par.HBOf, 2)+
	      (pow(par.HBOQ, 2)*pow((simputpsd->frequency[ii]-par.HBOf), 2))));
          simputpsd->power[ii]+=LHBO[ii];
        }
      }

      // QPO1 Lorentzian.
      if(par.Q1f!=0) {
        float Q1Norm=par.Q1rms/sqrt(0.5-(atan(par.Q1Q*(-1))/M_PI));
        LQ1=(float*)malloc(par.PSDnpt*sizeof(float));
	CHECK_NULL_BREAK(LQ1, status, "memory allocation failed");
        for(ii=0; ii<par.PSDnpt; ii++) {
          LQ1[ii]=(1.0/M_PI)*
	    ((pow(Q1Norm, 2)*par.Q1Q*par.Q1f)/
	     (pow(par.Q1f, 2)+
	      (pow(par.Q1Q, 2)*pow((simputpsd->frequency[ii]-par.Q1f), 2))));
          simputpsd->power[ii]+=LQ1[ii];
        }
      }

      // QPO2 Lorentzian.
      if(par.Q2f!=0) {
        float Q2Norm=par.Q2rms/sqrt(0.5-(atan(par.Q2Q*(-1))/M_PI));
        LQ2=(float*)malloc(par.PSDnpt*sizeof(float));
	CHECK_NULL_BREAK(LQ2, status, "memory allocation failed");
        for(ii=0; ii<par.PSDnpt; ii++) {
          LQ2[ii]=(1.0/M_PI)*
	    ((pow(Q2Norm, 2)*par.Q2Q*par.Q2f)/
	     (pow(par.Q2f, 2)+
	      (pow(par.Q2Q, 2)*pow((simputpsd->frequency[ii]-par.Q2f), 2))));
          simputpsd->power[ii]+=LQ2[ii];
        }
      }

      // QPO3 Lorentzian.
      if(par.Q3f!=0) {
        float Q3Norm=par.Q3rms/sqrt(0.5-(atan(par.Q3Q*(-1))/M_PI));
        LQ3=(float*)malloc(par.PSDnpt*sizeof(float));
	CHECK_NULL_BREAK(LQ3, status, "memory allocation failed");
        for(ii=0; ii<par.PSDnpt; ii++) {
          LQ3[ii]=(1.0/M_PI)*
	    ((pow(Q3Norm, 2)*par.Q3Q*par.Q3f)/
	     (pow(par.Q3f, 2)+
	      (pow(par.Q3Q, 2)*pow((simputpsd->frequency[ii]-par.Q3f), 2))));
          simputpsd->power[ii]+=LQ3[ii];
        }
      }

      // Release memory.
      free(Lzero); 
      free(LHBO); 
      free(LQ1); 
      free(LQ2); 
      free(LQ3);
    }

    // Store the PSD in the SIMPUT file.
    saveSimputPSD(simputpsd, par.Simput, "TIMING", 1, &status);
    CHECK_STATUS_BREAK(status);


    // Open the SimputCtlg.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Set the timing reference in the source catalog.
    char* timeref=(char*)malloc(32*sizeof(char));
    CHECK_NULL_BREAK(timeref, status, "memory allocation failed");
    strcpy(timeref, "[TIMING,1]");
    fits_write_col(cat->fptr, TSTRING, cat->ctiming, 1, 1, 1, 
    		   &timeref, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Close open files.
  if (NULL!=asciipsd) {
    fclose(asciipsd);
    asciipsd=NULL;
  }

  // Release memory.
  freeSimputPSD(&simputpsd);
  freeSimputCtlg(&cat, &status);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int simputpsd_getpar(struct Parameters* const par)
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

  return(status);
}


