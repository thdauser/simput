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
*/

#include "simputmerge.h"
#include "simput_tree.h"

// TODO: - add PSD to possible extensions?

int GLOBAL_SPEC_COUNTER=0;
int GLOBAL_IMG_COUNTER=0;
int GLOBAL_LC_COUNTER=0;
int GLOBAL_COUNTER=0;

static void get_infile_names(char ** infilenames, int num_cat, struct Parameters par, int *status){
	// Open the input catalogs. TODO: Need general routine
	infilenames[0]= (char*)malloc((strlen(par.Infile1)+1)*sizeof(char));
	CHECK_NULL_VOID(infilenames[0], *status, "memory allocation failed");
	infilenames[1]= (char*)malloc((strlen(par.Infile2)+1)*sizeof(char));
	CHECK_NULL_VOID(infilenames[1], *status, "memory allocation failed");

	if(num_cat != 2){
		*status=EXIT_FAILURE;
		printf(" num_cat!=2 : not implemented yet \n");
	}

	strcpy(infilenames[0], par.Infile1);
	strcpy(infilenames[1], par.Infile2);
	return ;
}
static void show_progress(SimputCtlg* outcat, SimputCtlg** incat){
	// Output of progress.
	// TODO: need this for arbitrary number num_cat
	if (0==outcat->nentries % 1000) {
		headas_chat(1, "\r%ld/%ld (%.1lf%%) entries",
				outcat->nentries, incat[0]->nentries+incat[1]->nentries,
				outcat->nentries*100./(incat[0]->nentries+incat[1]->nentries));
		fflush(NULL);
	}
	return;
}
static char* get_unused_ref(int type){

	char ref[SIMPUT_MAXSTR];

	switch (type) {
	case SIMPUT_SPEC_TYPE:
		sprintf(ref,"[SPECTRUM,1][NAME=='spec_%010i']",GLOBAL_SPEC_COUNTER);
		break;
	case SIMPUT_IMG_TYPE:
		sprintf(ref,"[IMG_%010i,1]",GLOBAL_IMG_COUNTER);
		break;
	case (SIMPUT_LC_TYPE):
		sprintf(ref,"[TIM_%010i,1]",GLOBAL_LC_COUNTER);
		break;
	case (SIMPUT_PSD_TYPE):
		sprintf(ref,"[TIM_%010i,1]",GLOBAL_LC_COUNTER);
		break;
	}

	return strdup(ref);
}

/**
static int compare_unique_ident(uniqueSimputident *id1, uniqueSimputident *id2){
	if (id1->io_pos != id2->io_pos){
		return -1;
	}

	if (strcmp(id1->filename,id2->filename) != 0){
		return -1;
	} else {
		return 0;
	}
} **/


int comp_elem(void* d1,void* d2){
	tdata *cd1 = (tdata*)d1;
	tdata *cd2 = (tdata*)d2;
	return strcmp(cd1->orig_ref,cd2->orig_ref);
}

void free_elem(void *d){
	tdata *cd = (tdata*)d;
	free(cd->orig_ref);
	free(cd->ref);
	free(cd);
	return;
}

/**
static simput_data* find_src_in_buffer(char *ref, simput_data** buf, int *status){


	for (int ii=0; ii< GLOBAL_COUNTER; ii++){
		if (strcmp(buf[ii]->orig_ref,ref) == 0){
			return buf[ii];
		}
	}
	CHECK_STATUS_RET(*status,NULL);
	// if we do not find the src in the buffer, we return NULL
	return NULL;
}
**/
static void append_src_to_buffer(simput_data* src,simput_data*** buf,int *status){
	GLOBAL_COUNTER++;
	*buf = (simput_data**) realloc(*buf,GLOBAL_COUNTER * sizeof(simput_data*));
	CHECK_NULL_VOID(*buf, *status, "memory (re)allocation failed")

	(*buf)[GLOBAL_COUNTER-1] = src;
	return;
}

static tdata* get_empty_elem(char* ref,int *status){
	tdata *src = NULL;
	src = (tdata*) malloc (sizeof(tdata));
	CHECK_NULL_RET(src,*status,"memory allocation failed",NULL);
	src->orig_ref = strdup(ref);
	return src;
}

static char* add_single_data_to_buffer(char* filename, char *ref,
		tree* mtree, simput_data*** buf,int type, int *status){

	char new_ref[SIMPUT_MAXSTR];

	// first make sure we have the absolute path:
	if ('['==ref[0]) {
		strcpy(new_ref, filename);
		strcat(new_ref, ref);
	} else {
		strcpy(new_ref, ref);
	}

	// this information is not used currently, as the function is not working yet
	// uniqueSimputident *ident = get_simput_ident(new_ref,type,status);

	tdata *new_elem = get_empty_elem(new_ref,status);

	node *ptr=NULL;
	ptr = find_elmt(mtree, new_elem);

	if(	 ptr == NULL){
		// load source

		simput_data* src = (simput_data*) malloc (sizeof(simput_data));
		switch (type) {

		case SIMPUT_SPEC_TYPE:
			GLOBAL_SPEC_COUNTER++;
			src->data = (SimputMIdpSpec*) malloc(sizeof(SimputMIdpSpec));
			CHECK_NULL_RET(src->data,*status,"memory allocation failed",NULL);
			SimputMIdpSpec* buf = loadSimputMIdpSpec(new_ref,status);
			sprintf(buf->name,"spec_%010i",GLOBAL_SPEC_COUNTER);
			CHECK_STATUS_RET(*status,NULL);
			src->data = buf;
			break;

		case SIMPUT_IMG_TYPE:
			GLOBAL_IMG_COUNTER++;
			src->data = (SimputImg*) malloc(sizeof(SimputImg));
			CHECK_NULL_RET(src->data,*status,"memory allocation failed",NULL);
			src->data = loadSimputImg(new_ref, status);
			CHECK_STATUS_RET(*status,NULL);
			break;

		case SIMPUT_LC_TYPE:
			GLOBAL_LC_COUNTER++;
			src->data = (SimputLC*) malloc(sizeof(SimputLC));
			CHECK_NULL_RET(src->data,*status,"memory allocation failed",NULL);
			src->data = loadSimputLC(new_ref, status);
			CHECK_STATUS_RET(*status,NULL);
			break;

		case SIMPUT_PSD_TYPE:
			GLOBAL_LC_COUNTER++;   // use the same counter??? (extension is named the same!
			src->data = (SimputPSD*) malloc(sizeof(SimputPSD));
			CHECK_NULL_RET(src->data,*status,"memory allocation failed",NULL);
			src->data = loadSimputPSD(new_ref, status);
			CHECK_STATUS_RET(*status,NULL);
			break;
		}
		src->type=type;
		// src->ident = ident;

		// need to find an unused reference for this data
		new_elem->ref = get_unused_ref(type);
		src->ref = strdup(new_elem->ref);

		ptr = add_elmt(mtree, new_elem);
		CHECK_NULL_RET(ptr,*status,"failed to add element",NULL);

		append_src_to_buffer(src,buf,status);

	} else {
		//free(new_elem);
		//tdata *
		new_elem = (tdata*)ptr->data;
	}
	return ((tdata*)new_elem)->ref;
}

static int is_fileref_given(char *str){
	if ((strcmp(str,"") == 0) || (strcmp(str,"NULL") == 0) ){
		return 0;
	} else {
		return 1;
	}
}

simput_refs* add_data_to_buffer(SimputSrc* insrc, SimputCtlg* incat, tree* mtree,
		simput_data*** data_buffer, int *status){

	simput_refs* refs=(simput_refs*)malloc(sizeof(simput_refs));
	CHECK_NULL_RET(refs, *status,"memory allocation for SimputSrc failed", NULL);

	// need to do this for all possible extensions (lightcurves as well?)

	// SPECTRUM
	if (is_fileref_given(insrc->spectrum) == 1){
		refs->spectrum = add_single_data_to_buffer(incat->filename,insrc->spectrum,
				mtree,data_buffer,SIMPUT_SPEC_TYPE,status);
		CHECK_STATUS_RET(*status,NULL);
	} else {
		refs->spectrum = strdup(insrc->spectrum);
	}

	// IMAGE
	if (is_fileref_given(insrc->image) == 1){
		refs->image = add_single_data_to_buffer(incat->filename,insrc->image,
				mtree,data_buffer,SIMPUT_IMG_TYPE,status);
		CHECK_STATUS_RET(*status,NULL);
	} else {
		refs->image = strdup(insrc->image);
	}

	// LIGHTCURVE or PSD
	if (is_fileref_given(insrc->timing) == 1){

		fitsfile* fptr=NULL;

		char new_ref[SIMPUT_MAXSTR];
		strcpy(new_ref, incat->filename);
		strcat(new_ref, insrc->timing);

		fits_open_table(&fptr,new_ref, READONLY, status);
	    if (EXIT_SUCCESS!=*status) {
	      char msg[SIMPUT_MAXSTR];
	      sprintf(msg, "could not open FITS table in file '%s'", insrc->timing);
	      SIMPUT_ERROR(msg);
	      return NULL;
	    }
	    int opt_status, ncol;
	    fits_get_colnum(fptr, CASEINSEN, "FLUX", &ncol, &opt_status);
		if (NULL!=fptr) fits_close_file(fptr, status);


		if (opt_status==EXIT_SUCCESS) {
			refs->timing = add_single_data_to_buffer(incat->filename,insrc->timing,
					mtree,data_buffer,SIMPUT_LC_TYPE,status);
		} else { // assume oif not LC it's PSD
			refs->timing = add_single_data_to_buffer(incat->filename,insrc->timing,
					mtree,data_buffer,SIMPUT_PSD_TYPE,status);
		}
		CHECK_STATUS_RET(*status,NULL);
	} else {
		refs->timing = strdup(insrc->timing);
	}

	return refs;
}

static char* extract_extname_link(char *str,int *status){

	char *pch;
	pch = strchr(str,']');

	CHECK_NULL_RET(pch,*status,"getting extension substring failed",NULL);

	int n = pch-str+1;
	char* substr;

	substr= (char*)malloc((n+1)*sizeof(char));
	CHECK_NULL_RET(substr, *status, "memory allocation failed",NULL);

	strncpy(substr,str,n);
	substr[n] = '\0';

	return substr;
}
static char* extract_extname(char *str_init,int *status){

	char *str;
	str = extract_extname_link(str_init,status);

	char *pch1,*pch2;
	pch2 = strchr(str,',');
	pch1 = strchr(str,'[');

	CHECK_NULL_RET(pch1,*status,"getting extension substring failed",NULL);
	CHECK_NULL_RET(pch2,*status,"getting extension substring failed",NULL);

	int n = pch2-pch1-1;
	char* substr;

	substr= (char*)malloc((n+1)*sizeof(char));
	CHECK_NULL_RET(substr, *status, "memory allocation failed",NULL);

	strncpy(substr,pch1+1,n);
	substr[n] = '\0';

	return substr;
}

static void write_single_merge_data(simput_data* buf, char* filename, int *status){

	int extver = 1;

	switch (buf->type){
	case SIMPUT_SPEC_TYPE:
		saveSimputMIdpSpec(buf->data, filename,
					extract_extname(buf->ref,status),
					extver, status);
		break;
	case SIMPUT_IMG_TYPE:
		saveSimputImg(buf->data, filename,
					extract_extname(buf->ref,status),
					extver, status);
		break;
	case SIMPUT_LC_TYPE:
		saveSimputLC(buf->data, filename,
				extract_extname(buf->ref,status),
				extver, status);
		break;

	case SIMPUT_PSD_TYPE:
		saveSimputPSD(buf->data, filename,
				extract_extname(buf->ref,status),
				extver, status);
		break;

	}

}

static void write_merge_data(simput_data** buf, char *filename, int *status){

	for (int ii=0; ii<GLOBAL_COUNTER;ii++){
		write_single_merge_data(buf[ii],filename,status);
	}

}

static void freeSimputData(simput_data** dat){
	if (dat != NULL){
		for (int ii=0; ii<GLOBAL_COUNTER; ii++){
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


int simputmerge_main() 
{
	// Program parameters.
	struct Parameters par;

	// Filenames of the input catalogs.
	const int num_cat = 2; // TODO: only hard coded for now
	char* infilenames[num_cat];


	// SIMPUT source catalogs.
	SimputCtlg* outcat  = NULL;
	SimputCtlg* incat[num_cat];

	// Array of already used IDs.
	long src_id=1;

	// Error status.
	int status=EXIT_SUCCESS;

	// Register HEATOOL
	set_toolname("simputmerge");
	set_toolversion("0.03");

	simput_data** data_buffer=NULL;
	tree* mtree=get_tree(&comp_elem, &free_elem);

	do { // Beginning of ERROR HANDLING Loop.

		// ---- Initialization ----

		// Read the parameters using PIL.
		status=simputmerge_getpar(&par);
		CHECK_STATUS_BREAK(status);

		get_infile_names(infilenames, num_cat, par, &status);

		// Get an empty output catalog. (TODO: Do proper checK??)
		remove(par.Outfile);
		outcat=openSimputCtlg(par.Outfile, READWRITE, 0, 0, 0, 0, &status);
		CHECK_STATUS_BREAK(status);

		// Loop over both source catalogs.
		for (int ii=0; ii<num_cat; ii++) {
			// Loop over all entries in the source catalog.

			incat[ii]=openSimputCtlg(infilenames[ii], READONLY, 0, 0, 0, 0, &status);
			CHECK_STATUS_BREAK(status);

			long jj;
			for (jj=0; jj<incat[ii]->nentries; jj++) {

				SimputSrc* insrc=getSimputSrc(incat[ii], jj+1, &status);
				CHECK_STATUS_BREAK(status);

				// Check whether the extensions should remain in their current
				// place or if they should by copied to the new output file.
				if (0==par.FetchExtensions) {
					printf("Not implemented yet!\n");
					break;
				}

				simput_refs *refs=NULL;
				refs = add_data_to_buffer(insrc,incat[ii], mtree, &data_buffer, &status);
				CHECK_NULL_BREAK(refs,status,"adding data to buffer failed");

				// Copy the entry from the input to the output catalog.
				SimputSrc* outsrc=newSimputSrcV(src_id,
						insrc->src_name,
						insrc->ra,
						insrc->dec,
						insrc->imgrota,
						insrc->imgscal,
						insrc->e_min,
						insrc->e_max,
						insrc->eflux,
						refs->spectrum, refs->image, refs->timing,
						&status);
				CHECK_STATUS_BREAK(status);

				appendSimputSrc(outcat, outsrc, &status);
				CHECK_STATUS_BREAK(status);

				// needs to be a unique identifier
				src_id++;

				show_progress(outcat,incat);
			}
			CHECK_STATUS_BREAK(status);
			// END of loop over all entries in the source catalog.
		}
		CHECK_STATUS_BREAK(status);
		headas_chat(1, "\n");
		// END of loop over both source catalogs.

		// Copy the used extensions to the new output file.
		if (0!=par.FetchExtensions) {
			write_merge_data(data_buffer,par.Outfile,&status);
		}
		// END of copy extensions to the new output file.

	} while(0); // END of ERROR HANDLING Loop.

	// --- Clean up ---
	headas_chat(3, "\ncleaning up ...\n");

	freeSimputData(data_buffer);

	for(int ii=0;ii<num_cat;ii++){
		freeSimputCtlg(&incat[ii], &status);
	}
	freeSimputCtlg(&outcat, &status);


	if (EXIT_SUCCESS==status) {
		headas_chat(3, "finished successfully!\n\n");
		return(EXIT_SUCCESS);
	} else {
		return(EXIT_FAILURE);
	}
}


int simputmerge_getpar(struct Parameters* const par)
{
	// String input buffer.
	char* sbuffer=NULL;

	// Error status.
	int status=EXIT_SUCCESS;

	// Read all parameters via the ape_trad_ routines.

	status=ape_trad_query_file_name("Infile1", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIMPUT_ERROR("failed reading the name of the input file 1");
		return(status);
	}
	strcpy(par->Infile1, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Infile2", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIMPUT_ERROR("failed reading the name of the input file 2");
		return(status);
	}
	strcpy(par->Infile2, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Outfile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIMPUT_ERROR("failed reading the name of the output file");
		return(status);
	}
	strcpy(par->Outfile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_bool("FetchExtensions", &par->FetchExtensions);
	if (EXIT_SUCCESS!=status) {
		SIMPUT_ERROR("failed reading the FetchExtensions parameter");
		return(status);
	}

	status=ape_trad_query_bool("clobber", &par->clobber);
	if (EXIT_SUCCESS!=status) {
		SIMPUT_ERROR("failed reading the clobber parameter");
		return(status);
	}

	return(status);
}

