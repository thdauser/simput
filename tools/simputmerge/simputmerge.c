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


   Copyright 2015 Thomas Dauser, FAU
   Copyright 2016-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "simputmerge.h"
#include "parinput.h"

// TODO: - add PSD to possible extensions?
// TODO: - Replcament for global counters (Do not use)
int GLOBAL_SPEC_COUNTER = 0;
int GLOBAL_IMG_COUNTER = 0;
int GLOBAL_LC_COUNTER = 0;
int GLOBAL_COUNTER = 0;

static char** get_infile_names(int* num_cat, struct Parameters par, int* status){

  // separates comma separated file list into individual filenames
	int infiles_capacity = 100; // inital capacity to include a maximum of 100 infilenames
  char** infilenames = malloc(infiles_capacity*sizeof(char*));
	CHECK_NULL_RET(infilenames, *status, "memory allocation for infilenames failed", NULL);
	int num_infiles_read = 0; // number of infiles read
  char* p = strtok(par.Infiles, ",");

  while (p != NULL)
  {
		if(num_infiles_read == infiles_capacity-1){
			infiles_capacity += 100;
			infilenames = realloc(infilenames, infiles_capacity*sizeof(char*));
			CHECK_NULL_RET(infilenames, *status, "memory reallocation for infilenames failed", NULL);
		}
		infilenames[num_infiles_read] = malloc((strlen(p)+1)*sizeof(char));
		CHECK_NULL_RET(infilenames[num_infiles_read], *status, "memory allocation failed", NULL);
		strcpy(infilenames[num_infiles_read], p);
    p = strtok(NULL, ",");
		num_infiles_read++;
  }

	*num_cat = num_infiles_read;
	return infilenames;
}

static void show_progress(SimputCtlg* outcat, int num_incat) {
	// Running from 0 to 100
	headas_chat(3, "\r%.0lf%%", outcat->nentries*100./num_incat);
	fflush(NULL);

	return;
}

static char* get_unused_ref(int type){

	char ref[SIMPUT_MAXSTR];

	switch (type) {
	case SIMPUT_SPEC_TYPE:
		sprintf(ref, "[SPECTRUM,1][NAME=='spec_%010i']", GLOBAL_SPEC_COUNTER);
		break;
	case SIMPUT_IMG_TYPE:
		sprintf(ref, "[IMG_%010i,1]", GLOBAL_IMG_COUNTER);
		break;
	case (SIMPUT_LC_TYPE):
		sprintf(ref, "[TIM_%010i,1]", GLOBAL_LC_COUNTER);
		break;
	case (SIMPUT_PSD_TYPE):
		sprintf(ref, "[TIM_%010i,1]", GLOBAL_LC_COUNTER);
		break;
	}

	return strdup(ref);
}

int comp_elem(void* d1, void* d2){
	tdata *cd1 = (tdata*)d1;
	tdata *cd2 = (tdata*)d2;
	return strcmp(cd1->orig_ref, cd2->orig_ref);
}

void free_elem(void* d){
	tdata* cd = (tdata*)d;
	free(cd->orig_ref);
	free(cd->ref);
	free(cd);
	return;
}

static void append_src_to_buffer(simput_data* src, simput_data*** buf, int* status){
	GLOBAL_COUNTER++;
	*buf = (simput_data**)realloc(*buf, GLOBAL_COUNTER * sizeof(simput_data*));
	CHECK_NULL_VOID(*buf, *status, "memory (re)allocation failed")

	(*buf)[GLOBAL_COUNTER-1] = src;
	return;
}

static tdata* get_empty_elem(char* ref, int* status){
	tdata *src = NULL;
	src = (tdata*)malloc(sizeof(tdata));
	CHECK_NULL_RET(src, *status, "memory allocation failed", NULL);
	src->orig_ref = strdup(ref);
	return src;
}

static char* add_single_data_to_buffer(char* filename, char* ref,
		tree* mtree, simput_data*** buf, int type, int* status){

	char new_ref[SIMPUT_MAXSTR];

	// first make sure we have the absolute path:
	if ('[' == ref[0]) {
		strcpy(new_ref, filename);
		strcat(new_ref, ref);
	} else {
		strcpy(new_ref, ref);
	}

	// this information is not used currently, as the function is not working yet
	// uniqueSimputident *ident = get_simput_ident(new_ref,type,status);

	tdata* new_elem = get_empty_elem(new_ref, status);

	node* ptr = NULL;
	ptr = find_elmt(mtree, new_elem);

	if(ptr == NULL){
		// load source

		simput_data* src = (simput_data*)malloc(sizeof(simput_data));
		switch (type) {

		case SIMPUT_SPEC_TYPE:
			GLOBAL_SPEC_COUNTER++;
			src->data = (SimputMIdpSpec*)malloc(sizeof(SimputMIdpSpec));
			CHECK_NULL_RET(src->data, *status, "memory allocation failed", NULL);
			SimputMIdpSpec* buf = loadSimputMIdpSpec(new_ref, status);
			buf->name = realloc(buf->name, sizeof(char)*100);
			CHECK_NULL_RET(buf->name, *status, "memory allocation failed", NULL);

			sprintf(buf->name, "spec_%010i", GLOBAL_SPEC_COUNTER);
			CHECK_STATUS_RET(*status, NULL);
			src->data = buf;
			break;

		case SIMPUT_IMG_TYPE:
			GLOBAL_IMG_COUNTER++;
			src->data = (SimputImg*)malloc(sizeof(SimputImg));
			CHECK_NULL_RET(src->data, *status, "memory allocation failed", NULL);
			src->data = loadSimputImg(new_ref, status);
			CHECK_STATUS_RET(*status, NULL);
			break;

		case SIMPUT_LC_TYPE:
			GLOBAL_LC_COUNTER++;
			src->data = (SimputLC*)malloc(sizeof(SimputLC));
			CHECK_NULL_RET(src->data, *status, "memory allocation failed", NULL);
			src->data = loadSimputLC(new_ref, status);
			CHECK_STATUS_RET(*status, NULL);
			break;

		case SIMPUT_PSD_TYPE:
			GLOBAL_LC_COUNTER++;   // use the same counter??? (extension is named the same!
			src->data = (SimputPSD*)malloc(sizeof(SimputPSD));
			CHECK_NULL_RET(src->data, *status, "memory allocation failed", NULL);
			src->data = loadSimputPSD(new_ref, status);
			CHECK_STATUS_RET(*status, NULL);
			break;
		}
		src->type = type;
		// src->ident = ident;

		// need to find an unused reference for this data
		new_elem->ref = get_unused_ref(type);
		src->ref = strdup(new_elem->ref);

		ptr = add_elmt(mtree, new_elem);
		CHECK_NULL_RET(ptr, *status, "failed to add element", NULL);

		append_src_to_buffer(src, buf, status);

	} else {
		//free(new_elem);
		//tdata *
		new_elem = (tdata*)ptr->data;
	}
	return ((tdata*)new_elem)->ref;
}

static int is_fileref_given(char* str){
	if ((strcmp(str, "") == 0) || (strcmp(str, "NULL") == 0) ){
		return 0;
	} else {
		return 1;
	}
}

simput_refs* add_data_to_buffer(SimputSrc* insrc, SimputCtlg* incat, tree* mtree,
		simput_data*** data_buffer, int* status){

	simput_refs* refs = (simput_refs*)malloc(sizeof(simput_refs));
	CHECK_NULL_RET(refs, *status, "memory allocation for SimputSrc failed", NULL);

        // get the full filename (with path)
	char fname[SIMPUT_MAXSTR];
        strcpy(fname, incat->filepath);
        strcat(fname, incat->filename);

	// need to do this for all possible extensions (lightcurves as well?)

	// SPECTRUM
	if (is_fileref_given(insrc->spectrum) == 1){
		refs->spectrum = add_single_data_to_buffer(fname, insrc->spectrum, mtree, data_buffer, SIMPUT_SPEC_TYPE, status);
		CHECK_STATUS_RET(*status, NULL);
	} else {
		refs->spectrum = strdup(insrc->spectrum);
	}

	// IMAGE
	if (is_fileref_given(insrc->image) == 1){
		refs->image = add_single_data_to_buffer(fname, insrc->image, mtree, data_buffer, SIMPUT_IMG_TYPE, status);
		CHECK_STATUS_RET(*status, NULL);
	} else {
		refs->image = strdup(insrc->image);
	}

	// LIGHTCURVE or PSD
	if (is_fileref_given(insrc->timing) == 1){

		fitsfile* fptr = NULL;

		char new_ref[SIMPUT_MAXSTR];
		strcpy(new_ref, fname);
		strcat(new_ref, insrc->timing);

		fits_open_table(&fptr, new_ref, READONLY, status);
	    if (EXIT_SUCCESS != *status) {
	      char msg[SIMPUT_MAXSTR];
	      sprintf(msg, "could not open FITS table in file '%s'", insrc->timing);
	      SIMPUT_ERROR(msg);
	      return NULL;
	    }
	    int opt_status = EXIT_SUCCESS;
	    int ncol;
	    fits_get_colnum(fptr, CASEINSEN, "FLUX", &ncol, &opt_status);

	    if (NULL != fptr) fits_close_file(fptr, status);

		if (opt_status == EXIT_SUCCESS) {
			refs->timing = add_single_data_to_buffer(fname, insrc->timing, mtree, data_buffer, SIMPUT_LC_TYPE, status);
		} else { // assume oif not LC it's PSD
			refs->timing = add_single_data_to_buffer(fname, insrc->timing, mtree, data_buffer, SIMPUT_PSD_TYPE, status);
		}
		CHECK_STATUS_RET(*status, NULL);
	} else {
		refs->timing = strdup(insrc->timing);
	}

	return refs;
}

static char* extract_extname_link(char* str, int* status){

	char* pch;
	pch = strchr(str, ']');

	CHECK_NULL_RET(pch, *status, "getting extension substring failed", NULL);

	int n = pch-str+1;
	char* substr;

	substr = (char*)malloc((n+1)*sizeof(char));
	CHECK_NULL_RET(substr, *status, "memory allocation failed", NULL);

	strncpy(substr, str, n);
	substr[n] = '\0';

	return substr;
}
static char* extract_extname(char *str_init, int* status){

	char* str;
	str = extract_extname_link(str_init, status);

	char*pch1,*pch2;
	pch2 = strchr(str, ',');
	pch1 = strchr(str, '[');

	CHECK_NULL_RET(pch1, *status, "getting extension substring failed", NULL);
	CHECK_NULL_RET(pch2, *status, "getting extension substring failed", NULL);

	int n = pch2-pch1-1;
	char* substr;

	substr= (char*)malloc((n+1)*sizeof(char));
	CHECK_NULL_RET(substr, *status, "memory allocation failed", NULL);

	strncpy(substr, pch1+1, n);
	substr[n] = '\0';

	return substr;
}

static void write_single_merge_data(simput_data* buf, char* filename, int* status){

	int extver = 1;

	switch (buf->type){
	case SIMPUT_SPEC_TYPE:
		saveSimputMIdpSpec(buf->data, filename,
					extract_extname(buf->ref, status),
					extver, status);
		break;
	case SIMPUT_IMG_TYPE:
		saveSimputImg(buf->data, filename,
					extract_extname(buf->ref, status),
					extver, status);
		break;
	case SIMPUT_LC_TYPE:
		saveSimputLC(buf->data, filename,
				extract_extname(buf->ref, status),
				extver, status);
		break;

	case SIMPUT_PSD_TYPE:
		saveSimputPSD(buf->data, filename,
				extract_extname(buf->ref, status),
				extver, status);
		break;

	}

}

static void write_merge_data(simput_data** buf, char* filename, int* status){

	for (int ii = 0; ii < GLOBAL_COUNTER; ii++){
		write_single_merge_data(buf[ii], filename, status);
	}

}

static void freeSimputData(simput_data** dat){
	if (dat != NULL){
		for (int ii = 0; ii < GLOBAL_COUNTER; ii++){
			if (dat[ii] != NULL){
				switch (dat[ii]->type) {
				case SIMPUT_SPEC_TYPE:
					// SimputMIdpSpec** buf = &(dat[ii]->data);
					freeSimputMIdpSpec( (SimputMIdpSpec**) &(dat[ii]->data) );
					break;
				case SIMPUT_IMG_TYPE:
					freeSimputImg( (SimputImg**) &(dat[ii]->data));
					break;
				case SIMPUT_LC_TYPE:
					freeSimputLC( (SimputLC**) &(dat[ii]->data));
					break;
				case SIMPUT_PSD_TYPE:
					freeSimputPSD( (SimputPSD**) &(dat[ii]->data));
					break;
				}
				free(dat[ii]->ref);
			}
			free(dat[ii]);
		}
	}

}

void merge_simput_files(int incat_buffer_length, SimputCtlg** incat_buffer, SimputCtlg* outcat, char* outfile, char FetchExtensions,
												long* src_id, int* status){

		simput_data** data_buffer = NULL;
		tree* mtree = get_tree(&comp_elem, &free_elem);
		do{
			// Loop over all source catalogs.
			for (int ii = 0; ii < incat_buffer_length; ii++) {
				for (long jj = 0; jj < incat_buffer[ii]->nentries; jj++) {

					SimputSrc* insrc = getSimputSrc(incat_buffer[ii], jj+1, status);
					CHECK_STATUS_BREAK(*status);

					// Check whether the extensions should remain in their current
					// place or if they should by copied to the new output file.
					if (0==FetchExtensions) {
						printf("Not implemented yet!\n");
						break;
					}

					simput_refs* refs=NULL;
					refs = add_data_to_buffer(insrc, incat_buffer[ii], mtree, &data_buffer, status);
					CHECK_NULL_BREAK(refs, *status, "adding data to buffer failed");

					// Copy the entry from the input to the output catalog.
					SimputSrc* outsrc=newSimputSrcV(*src_id,
							insrc->src_name,
							insrc->ra,
							insrc->dec,
							insrc->imgrota,
							insrc->imgscal,
							insrc->e_min,
							insrc->e_max,
							insrc->eflux,
							refs->spectrum, refs->image, refs->timing,
							status);
					CHECK_STATUS_BREAK(*status);

					appendSimputSrc(outcat, outsrc, status);
					CHECK_STATUS_BREAK(*status);

					// needs to be a unique identifier
					(*src_id)++;
				}
				CHECK_STATUS_BREAK(*status);
				// END of loop over all entries in the source catalog.
			}
			CHECK_STATUS_BREAK(*status);
			// END of loop over both source catalogs.

			// Copy the used extensions to the new output file.
			if (0 != FetchExtensions) {
				write_merge_data(data_buffer, outfile, status); //CHECK
			}
		}while(0);

		// Reset data buffer
		freeSimputData(data_buffer);

		//free mtree
		free_tree(mtree);

		// Reset global counter
		GLOBAL_COUNTER = 0;
}

void freeSimputFilenames(char* infile1, char* infile2, char* infiles, char** infilenames, int num_incat, char* outfile){

	freeSimputFile(&infile1);
	freeSimputFile(&infile2);
	freeSimputFile(&infiles);
	for (int ii = 0; ii < num_incat; ii++) {
		freeSimputFile(&infilenames[ii]);
	}
	freeSimputFile(&outfile);
}

int simputmerge_main()
{
	// Program parameters.
	struct Parameters par;

	// Filenames of the input catalogs.
	int num_incat = 0; // Total number of input catalogs (updated after user enters infile names)
	char** infilenames = NULL;

  // Number of input catalogs loaded and merged into output
	// Hard coded to 1 for now. Possible point for future optimization
	const int incat_buffer_length = 1;

	// SIMPUT source catalogs.
	SimputCtlg* outcat  = NULL;
	SimputCtlg* incat_buffer[incat_buffer_length];

	// Array of already used IDs.
	long src_id = 1;

	// Error status.
	int status = EXIT_SUCCESS;

	// Register HEATOOL
	set_toolname("simputmerge");
	set_toolversion("0.03");

	do { // Beginning of ERROR HANDLING Loop.

		// ---- Initialization ----

		// Read the parameters using PIL.
		status = simputmerge_getpar(&par);
		CHECK_STATUS_BREAK(status);

		infilenames = get_infile_names(&num_incat, par, &status);

		// Get an empty output catalog. (TODO: Do proper checK??)
		remove(par.Outfile);
		outcat = openSimputCtlg(par.Outfile, READWRITE, 0, 0, 0, 0, &status);
		CHECK_STATUS_BREAK(status);

		headas_chat(3, "\nmerging simput files ...\n");

		// Loop over all source catalogs. Always load and merge up
		// to #incat_buffer_length catalogs at once.
		for (int ii = 0; ii < num_incat; ii += incat_buffer_length) {
			int n_incat_loaded = 0;
			for (int jj = 0; jj < incat_buffer_length; jj++){
				// Check if catalog exists (the last batch might contain less than
				// incat_buffer_length entries).
				if (ii + jj < num_incat) {
					incat_buffer[jj] = openSimputCtlg(infilenames[ii + jj], READONLY, 0, 0, 0, 0, &status);
					CHECK_STATUS_BREAK(status);
					n_incat_loaded++;
				}
			}

			merge_simput_files(n_incat_loaded, incat_buffer, outcat, par.Outfile, par.FetchExtensions, &src_id, &status);
			show_progress(outcat, num_incat);

			for (int jj = 0; jj < n_incat_loaded; jj++){
				freeSimputCtlg(&incat_buffer[jj], &status);
			}
		}

	} while(0); // END of ERROR HANDLING Loop.

	// --- Clean up ---
	headas_chat(3, "\n\ncleaning up ...\n");

	for(int ii = 0; ii < incat_buffer_length; ii++){
		freeSimputCtlg(&incat_buffer[ii], &status);
	}
	freeSimputCtlg(&outcat, &status);
	freeSimputFilenames(par.Infile1, par.Infile2, par.Infiles, infilenames, num_incat, par.Outfile);

	if (EXIT_SUCCESS==status) {
		headas_chat(3, "finished successfully!\n\n");
		return(EXIT_SUCCESS);
	} else {
		return(EXIT_FAILURE);
	}
}

int simputmerge_getpar(struct Parameters* const par)
{

	// Error status.
	int status=EXIT_SUCCESS;

	// Read all parameters via the ape_trad_ routines.
	query_simput_parameter_file_name("Infile1", &(par->Infile1), &status);
	query_simput_parameter_file_name("Infile2", &(par->Infile2), &status);

	if (par->Infile1 && par->Infile2) {
		// combine the two infiles names into a comma separated file name
		size_t combined_size = strlen(par->Infile1) + strlen(par->Infile2) + 2; //+2 for comma and terminating NULL
		par->Infiles = malloc(combined_size*sizeof(char));
		CHECK_NULL_RET(par->Infiles, status, "memory allocation failed", EXIT_FAILURE);
		snprintf(par->Infiles, combined_size, "%s,%s", par->Infile1, par->Infile2);
	} else {
		query_simput_parameter_file_name("Infiles", &(par->Infiles), &status);
		}

  query_simput_parameter_file_name("Outfile", &(par->Outfile), &status);

	status = ape_trad_query_bool("FetchExtensions", &par->FetchExtensions);
	if (EXIT_SUCCESS != status) {
		SIMPUT_ERROR("failed reading the FetchExtensions parameter");
		return(status);
	}

	status = ape_trad_query_bool("clobber", &par->clobber);
	if (EXIT_SUCCESS != status) {
		SIMPUT_ERROR("failed reading the clobber parameter");
		return(status);
	}

	return(status);
}
