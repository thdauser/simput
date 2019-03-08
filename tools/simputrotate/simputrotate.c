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


   Copyright 2015 Thorsten Brand, FAU
   Copyright 2016-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "simputrotate.h"

int simputrotate_getpar(Parameters *par){

  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("Infile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("failed reading the name of the input file.");
    return(status);
  }
  strcpy(par->incat, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Outfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("failed reading the name of the output file.");
    return(status);
  }
  strcpy(par->outcat, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("C1_RA", &par->c1_ra);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the C1_RA parameter failed");
    return(status);
  }

  status=ape_trad_query_float("C1_Dec", &par->c1_dec);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the C1_Dec parameter failed");
    return(status);
  }

  status=ape_trad_query_float("C2_RA", &par->c2_ra);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the C2_RA parameter failed");
    return(status);
  }

  status=ape_trad_query_float("C2_Dec", &par->c2_dec);
  if (EXIT_SUCCESS!=status) {
    SIMPUT_ERROR("reading the C2_Dec parameter failed");
    return(status);
  }

  return(status);

}

int simputrotate_main()
{
  // Program parameters.
  Parameters par;

  fitsfile *ifptr=NULL;
  fitsfile *ofptr=NULL;
  float *ra=NULL;
  float *dec=NULL;
  float *rra=NULL;
  float *rdec=NULL;

  // Register HEATOOL
  set_toolname("simputrotate");
  set_toolversion("0.00");

  int status=EXIT_SUCCESS;

  do { // Beginning of ERROR HANDLING Loop.

    // Read the parameters using PIL.
    status=simputrotate_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Copy catalog as a whole
    if(fits_open_file(&ifptr, par.incat, READONLY, &status)){
      puts("***ERROR: Unable to open input file.");
      break;
    }

    if(fits_create_file(&ofptr, par.outcat, &status)){
      puts("***ERROR: Unable to open output file. Already exists?");
      break;
    }

    if(fits_copy_file(ifptr, ofptr, 1, 1, 1, &status)){
      puts("***ERROR: Unable to copy HDUs from the original catalog.");
      break;
    }

    // Move to the SRC_CAT-extension
    if(fits_movnam_hdu(ofptr, BINARY_TBL, "SRC_CAT", 0, &status)){
      puts("***ERROR: Unable to move to the SRC_CAT extension in the output file.");
      break;
    }

    // Read the RA, Dec-columns
    int cra, cdec;
    long nc;
    if(fits_get_colnum(ofptr, CASEINSEN, "RA", &cra, &status)){
      puts("***ERROR: Unable to find the RA-column.");
      break;
    }
    if(fits_get_colnum(ofptr, CASEINSEN, "DEC", &cdec, &status)){
      puts("***ERROR: Unable to find the DEC-column.");
      break;
    }
    if(fits_get_num_rows(ofptr, &nc, &status)){
      puts("***ERROR: Unable to read number of rows.");
      break;
    }

    ra=(float*)malloc(nc*sizeof(float));
    if(ra==NULL){
      puts("***ERROR: Unable to allocate memory.");
      break;
    }
    dec=(float*)malloc(nc*sizeof(float));
    if(dec==NULL){
      puts("***ERROR: Unable to allocate memory.");
      break;
    }
    rra=(float*)malloc(nc*sizeof(float));
    if(rra==NULL){
      puts("***ERROR: Unable to allocate memory.");
      break;
    }
    rdec=(float*)malloc(nc*sizeof(float));
    if(rdec==NULL){
      puts("***ERROR: Unable to allocate memory.");
      break;
    }

    float nulval=0.;
    int anynul=0;

    if(fits_read_col(ofptr, TFLOAT, cra, 1, 1, nc, &nulval, ra, &anynul, &status)){
      puts("***ERROR: Unable to read RA column.");
      break;
    }
    if(fits_read_col(ofptr, TFLOAT, cdec, 1, 1, nc, &nulval, dec, &anynul, &status)){
      puts("***ERROR: Unable to read DEC column.");
      break;
    }

    // Rotate the coordinates
    rotate_coord_system(par.c1_ra, par.c1_dec,
			par.c2_ra, par.c2_dec,
			ra, dec,
			rra,rdec,
			nc);

    // Overwrite the old coordinates
    if(fits_write_col(ofptr, TFLOAT, cra, 1, 1, nc, rra, &status)){
      puts("***ERROR: Unable to write RA column.");
      break;
    }
    if(fits_write_col(ofptr, TFLOAT, cdec, 1, 1, nc, rdec, &status)){
      puts("***ERROR: Unable to write DEC column.");
      break;
    }

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  headas_chat(3, "\ncleaning up ...\n");

  if(ofptr!=NULL){
    fits_close_file(ofptr, &status);
    ofptr=NULL;
  }
  if(ifptr!=NULL){
    fits_close_file(ifptr, &status);
    ifptr=NULL;
  }

  if(ra!=NULL){
    free(ra);
    ra=NULL;
  }
  if(dec!=NULL){
    free(dec);
    dec=NULL;
  }
  if(rra!=NULL){
    free(rra);
    rra=NULL;
  }
  if(rdec!=NULL){
    free(rdec);
    rdec=NULL;
  }

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }

}
