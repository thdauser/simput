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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "simputverify.h"


int simputverify_main()
{
  // Program parameters.
  struct Parameters par;

  // Simput data structures (used as buffers).
  SimputCtlg* cat =NULL;
  SimputMIdpSpec* spec=NULL;
  SimputImg* img=NULL;
  SimputLC* lc=NULL;
  SimputPSD* psd=NULL;
  SimputPhList* phl=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputverify");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputverify_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Open the catalog.
    cat=openSimputCtlg(par.Simput, READONLY, 0, 0, 0, 0, &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all entries in the source catalog.
    long ii;
    for (ii=0; ii<cat->nentries; ii++) {

    	// Flag, whether the source refers to a spectrum.
    	int has_spectrum=0;

    	// Load the source from the catalog.
    	SimputSrc* src=loadSimputSrc(cat, ii+1, &status);
    	CHECK_STATUS_BREAK(status);

    	// Try and load the spectrum.
    	char filename[SIMPUT_MAXSTR];
    	if ((strlen(src->spectrum)>0) &&
    			(0!=strcmp(src->spectrum, "NULL"))) {
    		if ('['==src->spectrum[0]) {
    			strcpy(filename, par.Simput);
    			strcat(filename, src->spectrum);
    		} else {
    			strcpy(filename, src->spectrum);
    		}


    		// Determine the type of the spectrum extension.
    		int exttype=getSimputExtType(cat, filename, &status);
    		CHECK_STATUS_BREAK(status);
    		switch (exttype) {
    		case EXTTYPE_MIDPSPEC:
    			spec=loadSimputMIdpSpec(filename, &status);
    			CHECK_STATUS_BREAK(status);
    			freeSimputMIdpSpec(&spec);
    			break;

    		case EXTTYPE_PHLIST:
    			phl=openSimputPhList(filename, READONLY, &status);
    			CHECK_STATUS_BREAK(status);
    			freeSimputPhList(&phl, &status);
    			CHECK_STATUS_BREAK(status);
    			break;

    		default:
    			status=EXIT_FAILURE;
    			char msg[SIMPUT_MAXSTR];
    			sprintf(msg, "spectrum reference '%s' in catalog line '%ld' "
    					"refers to unknown extension type", src->spectrum, ii+1);
    			SIMPUT_ERROR(msg);
    		}
    		CHECK_STATUS_BREAK(status);

    		has_spectrum=1;
    	}

    	// Try and load the image.
    	if ((strlen(src->image)>0) &&
    			(0!=strcmp(src->image, "NULL"))) {
    		if ('['==src->image[0]) {
    			strcpy(filename, par.Simput);
    			strcat(filename, src->image);
    		} else {
    			strcpy(filename, src->image);
    		}

    		// Determine the type of the spectrum extension.
    		int exttype=getSimputExtType(cat, filename, &status);
    		CHECK_STATUS_BREAK(status);
    		switch (exttype) {
    		case EXTTYPE_IMAGE:
    			img=loadSimputImg(filename, &status);
    			CHECK_STATUS_BREAK(status);
    			freeSimputImg(&img);
    			break;

    		case EXTTYPE_PHLIST:
    			phl=openSimputPhList(filename, READONLY, &status);
    			CHECK_STATUS_BREAK(status);
    			freeSimputPhList(&phl, &status);
    			CHECK_STATUS_BREAK(status);
    			break;

    		default:
    			status=EXIT_FAILURE;
    			char msg[SIMPUT_MAXSTR];
    			sprintf(msg, "image reference '%s' in catalog line '%ld' "
    					"refers to unknown extension type", src->image, ii+1);
    			SIMPUT_ERROR(msg);
    		}
    		CHECK_STATUS_BREAK(status);
    	}

    	// Try and load the timing extension.
    	if ((strlen(src->timing)>0) &&
    			(0!=strcmp(src->timing, "NULL"))) {
    		if ('['==src->timing[0]) {
    			strcpy(filename, par.Simput);
    			strcat(filename, src->timing);
    		} else {
    			strcpy(filename, src->timing);
    		}

    		// Determine the type of the timing extension.
    		int exttype=getSimputExtType(cat, filename, &status);
    		CHECK_STATUS_BREAK(status);
    		switch (exttype) {
    		case EXTTYPE_LC:
    			lc=loadSimputLC(filename, &status);
    			CHECK_STATUS_BREAK(status);

    			// Check whether the light curve refers to further SPECTRUM
    			// and IMAGE extensions.
    			if (NULL!=lc->spectrum) {
    				long jj;
    				for (jj=0; jj<lc->nentries; jj++) {
    					if ((strlen(lc->spectrum[jj])>0) &&
    							(0!=strcmp(lc->spectrum[jj], "NULL"))) {
    						char filename2[SIMPUT_MAXSTR];
    						if ('['==lc->spectrum[jj][0]) {
    							strcpy(filename2, lc->fileref);
    							char* firstbracket=strchr(filename2, '[');
    							if (NULL!=firstbracket) {
    								firstbracket[0]='\0';
    							}
    							strcat(filename2, lc->spectrum[jj]);
    						} else {
    							strcpy(filename2, lc->spectrum[jj]);
    						}

    						// Determine extension type and open the extension.
    						int exttype=getSimputExtType(cat, filename2, &status);
    						CHECK_STATUS_BREAK(status);
    						switch (exttype) {
    						case EXTTYPE_MIDPSPEC:
    							spec=loadSimputMIdpSpec(filename2, &status);
    							CHECK_STATUS_BREAK(status);
    							freeSimputMIdpSpec(&spec);
    							break;

    						case EXTTYPE_PHLIST:
    							phl=openSimputPhList(filename2, READONLY, &status);
    							CHECK_STATUS_BREAK(status);
    							freeSimputPhList(&phl, &status);
    							CHECK_STATUS_BREAK(status);
    							break;

    						default:
    							status=EXIT_FAILURE;
    							char msg[SIMPUT_MAXSTR];
    							sprintf(msg, "spectrum reference '%s' in line '%ld' of "
    									"light curve '%s' refers to unknown extension type",
										lc->spectrum[jj], jj+1, filename);
    							SIMPUT_ERROR(msg);
    						}
    						CHECK_STATUS_BREAK(status);

    						has_spectrum=1;
    					}
    				}
    				CHECK_STATUS_BREAK(status);
    			}
    			if (NULL!=lc->image) {
    				long jj;
    				for (jj=0; jj<lc->nentries; jj++) {
    					if ((strlen(lc->image[jj])>0) &&
    							(0!=strcmp(lc->image[jj], "NULL"))) {
    						char filename2[SIMPUT_MAXSTR];
    						if ('['==lc->image[jj][0]) {
    							strcpy(filename2, lc->fileref);
    							char* firstbracket=strchr(filename2, '[');
    							if (NULL!=firstbracket) {
		    firstbracket[-1]='\0';
		  }
		  strcat(filename2, lc->image[jj]);
		} else {
		  strcpy(filename2, lc->image[jj]);
		}

		// Determine extension type and open the extension.
		int exttype=getSimputExtType(cat, filename2, &status);
		CHECK_STATUS_BREAK(status);
		switch (exttype) {
		case EXTTYPE_IMAGE:
		  img=loadSimputImg(filename2, &status);
		  CHECK_STATUS_BREAK(status);
		  freeSimputImg(&img);
		  break;

		case EXTTYPE_PHLIST:
		  phl=openSimputPhList(filename2, READONLY, &status);
		  CHECK_STATUS_BREAK(status);
		  freeSimputPhList(&phl, &status);
		  CHECK_STATUS_BREAK(status);
		  break;

		default:
		  status=EXIT_FAILURE;
		  char msg[SIMPUT_MAXSTR];
		  sprintf(msg, "image reference '%s' in line '%ld' of "
			  "light curve '%s' refers to unknown extension type",
			  lc->image[jj], jj+1, filename);
		  SIMPUT_ERROR(msg);
		}
		CHECK_STATUS_BREAK(status);
	      }
	    }
	    CHECK_STATUS_BREAK(status);
	  }

	  freeSimputLC(&lc);
	  break;

	case EXTTYPE_PSD:
	  psd=loadSimputPSD(filename, &status);
	  CHECK_STATUS_BREAK(status);
	  freeSimputPSD(&psd);
	  break;

	case EXTTYPE_PHLIST:
	  phl=openSimputPhList(filename, READONLY, &status);
	  CHECK_STATUS_BREAK(status);
	  freeSimputPhList(&phl, &status);
	  CHECK_STATUS_BREAK(status);
	  break;

	default:
	  status=EXIT_FAILURE;
	  char msg[SIMPUT_MAXSTR];
	  sprintf(msg, "timing reference '%s' in catalog line '%ld' "
		  "refers to unknown extension type", src->timing, ii+1);
	  SIMPUT_ERROR(msg);
	}
	CHECK_STATUS_BREAK(status);
      }

      // Check if there has been at least one SPECTRUM reference, either
      // in the source catalog or in the light curve.
      if (0==has_spectrum) {
	status=EXIT_FAILURE;
	char msg[SIMPUT_MAXSTR];
	sprintf(msg, "source in catalog line '%ld' has no spectrum", ii+1);
	SIMPUT_ERROR(msg);
	break;
      }
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all entries in the source catalog.

  } while(0); // END of ERROR HANDLING Loop.

  // --- Clean up ---
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeSimputCtlg(&cat, &status);
  freeSimputMIdpSpec(&spec);
  freeSimputImg(&img);
  freeSimputLC(&lc);
  freeSimputPSD(&psd);
  freeSimputPhList(&phl, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(0, "verification was successful!\n\n");
    return(EXIT_SUCCESS);
  } else {
    headas_chat(0, "verification failed!\n\n");
    return(EXIT_FAILURE);
  }
}


int simputverify_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("failed reading the name of the SIMPUT file");
    return(status);
  }
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  return(status);
}
