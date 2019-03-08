/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2015 Philippe Peille, IRAP
   Copyright 2016-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "simputmulticell.h"

/** Function to calibrate the parameter bins */
static void cal_par_arrays(par_info *par, struct param_input *ipar, ParamReader* reader, int num_param, int* status){

	int ii;
	double* min_array=(double*)malloc(reader->npara*sizeof(double));
	double* max_array=(double*)malloc(reader->npara*sizeof(double));
	reader->get_minmax(reader->paraminfo,min_array,max_array,status);

	for (ii=0; ii<num_param; ii++){

		ipar[ii].minPar = min_array[ii];
		ipar[ii].maxPar = max_array[ii];

		if (ipar[ii].minPar >= ipar[ii].maxPar){
			printf("for Param%i parMin (%.4f) >= parMax (%.4f)",ii+1,ipar[ii].minPar,ipar[ii].maxPar);
			*status=EXIT_FAILURE;
			break;
		}
		//        printf("\n------- %.3f %.3f (%i) -------- \n",ipar[ii].minPar,ipar[ii].maxPar,ipar[ii].num_values);

		par[ii].pvals = (double *)malloc(ipar[ii].num_values*(sizeof(double)));
		CHECK_MALLOC_VOID(par[ii].pvals);

		cal_single_par_array(ipar[ii].minPar,ipar[ii].maxPar,ipar[ii].num_values,
				par[ii].pvals,ipar[ii].logScale);

		// set the structure
		par[ii].num_param = num_param;
		par[ii].num_pvals = ipar[ii].num_values;
		par[ii].par_names = ipar[ii].param_names;

		print_par_array(par[ii], ii, ipar[ii].logScale);
	}
	free(min_array);
	free(max_array);

}

static void get_par_id_from_par_array(int id_array[], double* par_array, const struct param_input ipar[],
		const int num_param){

	int kk;
	for ( kk=0; kk < num_param; kk++) {
		// get the parameter ID
		id_array[kk] = get_single_par_id(ipar[kk].minPar,ipar[kk].maxPar,ipar[kk].num_values,
				par_array[kk],ipar[kk].logScale);

	}

}

/** Function to set up the param_input structure (save for all param the number of values and the name, initialize min/max) */
static void set_input_par(struct param_input *ipar, struct Parameters* par, int num_par, int *status) {

    char **pnames = NULL;
    char **pnumval = NULL;
    char **plog = NULL;

    int ii;

    pnames = parse_string2array(par->ParamNames,num_par,status);
    CHECK_NULL_VOID(pnames,*status,"Error when parsing ParamNames");

    for (ii=0; ii < num_par ; ii++ ){
        	ipar[ii].param_names = pnames[ii];
    }

    // we now do allow NULL here and then set the values to the standard value
    pnumval = parse_string2array(par->ParamsNumValues,num_par,status);
    for (ii=0; ii < num_par ; ii++ ){
    	if (pnumval == NULL){
    		ipar[ii].num_values = DEFAULT_NUM_VALUES;
    	} else {
    	    ipar[ii].num_values = (int) (atof(pnumval[ii]));
    	    // make sure that if a string (like "none") is given, we use the default value
    	    if (ipar[ii].num_values == 0 ){
    	        ipar[ii].num_values = DEFAULT_NUM_VALUES;
    	    }
    	}
    }

    plog = parse_string2array(par->ParamsLogScale,num_par,status);

    for (ii=0; ii < num_par ; ii++ ){
    	if (plog == NULL){
    		ipar[ii].logScale = 0;
    	} else {
    		if ( strcmp(plog[ii],"yes") == 0 ){
    			ipar[ii].logScale = 1;
    		} else {
    			ipar[ii].logScale = 0;
    		}
    	}
    }

    // set minial and maximal parameter value
    // TODO: do we want to allow them to be set in the par file?
    for (ii=0; ii < num_par ; ii++ ){
        ipar[0].minPar = 0.0;
    	ipar[0].maxPar = 0.0;
    }
    for (ii=0;ii<num_par;ii++){
    	if (pnumval!=NULL){
    		free(pnumval[ii]);
    	}
    	if (plog!=NULL){
    		free(plog[ii]);
    	}
    }
    free(pnumval);
    free(plog);
    free(pnames);

}


/** Open an image FITS file and prepare to read it as an ImageCellDistrib structure */
static ImageCellDistrib* newImageCellDistrib(char* paramFile,char* paramIndexes,int fluxIndex,long startIndex,int* status){
	// Malloc structure
	ImageCellDistrib* distrib=(ImageCellDistrib*)malloc(sizeof(*distrib));
	CHECK_MALLOC_RET_NULL(distrib);

	// Load image
	distrib->image= loadSimputImg_map(paramFile, status);
	CHECK_STATUS_RET(*status, NULL);

	// Initialize row counter
	distrib->row = startIndex-1;

	// Create the index array
	distrib->num_par = parse_num_param(paramIndexes);
	distrib->index_array = (int*)malloc(distrib->num_par*sizeof(int));
	CHECK_MALLOC_RET_NULL(distrib->index_array);
	char** pindexes = NULL;
	pindexes = parse_string2array(paramIndexes,distrib->num_par,status);
	CHECK_MALLOC_RET_NULL(pindexes);
	for (int i=0;i<distrib->num_par;i++){
		distrib->index_array[i] = atoi(pindexes[i]);
		free(pindexes[i]);
	}
	free(pindexes);
	distrib->flux_index = fluxIndex;
	return(distrib);
}

/** Function to get the parameters of the next cell from an image */
int getNextCellFromImage(void *paraminfo, double* para_array, double* ra, double* dec, double* flux, int* status){
	CHECK_STATUS_RET(*status,0);
	// Cast param info as an ImageCellDistrib
	ImageCellDistrib* distrib = (ImageCellDistrib*)paraminfo;

	// Get the next image row
	if (distrib->row<distrib->image->naxis1){
		for (int i=0;i<distrib->num_par;i++){
			para_array[i] = distrib->image->dist[distrib->row][distrib->index_array[i]];
		}
		*ra = distrib->image->dist[distrib->row][0];
		*dec = distrib->image->dist[distrib->row][1];
		*flux = distrib->image->dist[distrib->row][distrib->flux_index];
		distrib->row++;
		return(1);
	} else {
		return(0);
	}
}

/** Function to get the min and max values for all needed parameters */
void getMinMaxFromImage(void *paraminfo, double* min_array, double* max_array, int* status){
	// Cast param info as an ImageCellDistrib
	ImageCellDistrib* distrib = (ImageCellDistrib*)paraminfo;

	// Get image first line
	long current_row=distrib->row;
	distrib->row=0;
	double* para_array=(double*)malloc(distrib->num_par*sizeof(*para_array));
	CHECK_MALLOC_VOID(para_array);
	double ra=0,dec=0,flux=0;
	getNextCellFromImage(distrib,para_array,&ra,&dec,&flux,status);

	// Read all the data and find min/max for all parameters
	for (int i=0;i<distrib->num_par;i++){
		min_array[i] = para_array[i];
		max_array[i] = para_array[i];
	}
	while (getNextCellFromImage(distrib,para_array,&ra,&dec,&flux,status)){
		for (int i=0;i<distrib->num_par;i++){
			min_array[i] = min(min_array[i],para_array[i]);
			max_array[i] = max(max_array[i],para_array[i]);
		}
	}

	// Reboot image
	distrib->row=current_row;
	free(para_array);
}

/** Free a an ImageCellDistrib structure */
static void freeImageCellDistrib(void** paraminfo, int* status){
	CHECK_STATUS_VOID(*status);
	// Cast param info as an ImageCellDistrib
	ImageCellDistrib* distrib = (ImageCellDistrib*)*paraminfo;

	// Free memory
	if (distrib!=NULL){
		freeSimputImg(&(distrib->image));
		free(distrib->index_array);
		free(distrib);
		*paraminfo=NULL;
	}
}

/** Open an image FITS file and prepare to read it as an ImageCellDistrib structure */
static TableCellDistrib* newTableCellDistrib(char* paramFile,char* paramInputNames,long startIndex,int* status){
	// Malloc structure
	TableCellDistrib* distrib=(TableCellDistrib*)malloc(sizeof(*distrib));
	CHECK_MALLOC_RET_NULL(distrib);

	// Open table
	fits_open_file(&(distrib->fptr), paramFile, READONLY, status);
	CHECK_STATUS_RET(*status, NULL);

	//Move to the binary table
	fits_movnam_hdu(distrib->fptr,ANY_HDU,"SIMPARA",0, status);
	CHECK_STATUS_RET(*status, NULL);

	//Get number of rows
	char comment[FLEN_COMMENT];
	fits_read_key(distrib->fptr, TLONGLONG, "NAXIS2", &(distrib->nrows), comment, status);
	CHECK_STATUS_RET(*status, NULL);

	// Create the index array
	distrib->num_par = parse_num_param(paramInputNames);
	distrib->col_indexes = (int*)malloc(distrib->num_par*sizeof(int));
	CHECK_MALLOC_RET_NULL(distrib->col_indexes);
	char** pinput_names= NULL;
	pinput_names = parse_string2array(paramInputNames,distrib->num_par,status);
	CHECK_MALLOC_RET_NULL(pinput_names);

	//Associate column numbers
	for (int ii=0;ii<distrib->num_par;ii++){
		fits_get_colnum(distrib->fptr, CASEINSEN,pinput_names[ii], &(distrib->col_indexes[ii]), status);
		free(pinput_names[ii]);
	}
	free(pinput_names);
	CHECK_STATUS_RET(*status, NULL);
	fits_get_colnum(distrib->fptr, CASEINSEN,"FLUX", &(distrib->col_flux_index), status);
	fits_get_colnum(distrib->fptr, CASEINSEN,"RA", &(distrib->col_ra_index), status);
	fits_get_colnum(distrib->fptr, CASEINSEN,"DEC", &(distrib->col_dec_index), status);
	CHECK_STATUS_RET(*status, NULL);

	// Initialize row counter
	distrib->row =startIndex;

	return(distrib);
}

/** Function to get the parameters of the next cell from an image */
int getNextCellFromTable(void *paraminfo, double* para_array, double* ra, double* dec, double* flux, int* status){
	// Cast param info as an TableCellDistrib
	TableCellDistrib* distrib = (TableCellDistrib*)paraminfo;

	// Get the next image row
	int anynul=0;
	if (distrib->row<=distrib->nrows){
		for (int ii=0;ii<distrib->num_par;ii++){
			fits_read_col(distrib->fptr, TDOUBLE, distrib->col_indexes[ii],distrib->row,1,1,0,&(para_array[ii]), &anynul,status); // For the moment, it is much simpler to only consider doubles columns
		}
		CHECK_STATUS_RET(*status,0);
		fits_read_col(distrib->fptr, TDOUBLE, distrib->col_ra_index,distrib->row,1,1,0,ra, &anynul,status);
		fits_read_col(distrib->fptr, TDOUBLE, distrib->col_dec_index,distrib->row,1,1,0,dec, &anynul,status);
		fits_read_col(distrib->fptr, TDOUBLE, distrib->col_flux_index,distrib->row,1,1,0,flux, &anynul,status);
		CHECK_STATUS_RET(*status,0);
		distrib->row++;
		return(1);
	} else {
		return(0);
	}
}

/** Function to get the min and max values for all needed parameters */
void getMinMaxFromTable(void *paraminfo, double* min_array, double* max_array, int* status){
	// Cast param info as an ImageCellDistrib
	TableCellDistrib* distrib = (TableCellDistrib*)paraminfo;

	// Allocate column array
	double* para_array=(double*)malloc(distrib->nrows*sizeof(*para_array));
	CHECK_MALLOC_VOID(para_array);

	// Read each parameter column and get min/max from it
	int anynul=0;
	for (int ii=0;ii<distrib->num_par;ii++){
		fits_read_col(distrib->fptr, TDOUBLE, distrib->col_indexes[ii],1,1,distrib->nrows,0,para_array, &anynul,status); // For the moment, it is much simpler to only consider doubles columns
		CHECK_STATUS_VOID(*status);
		min_max(para_array, distrib->nrows, &min_array[ii], &max_array[ii]);
	}

	// Free memory
	free(para_array);
}

/** Free a TableCellDistrib structure */
static void freeTableCellDistrib(void** paraminfo, int* status){
	// Cast param info as an TableCellDistrib
	TableCellDistrib* distrib = (TableCellDistrib*)*paraminfo;

	// Free memory
	if (distrib!=NULL){
		fits_close_file(distrib->fptr, status);
		free(distrib->col_indexes);
		free(distrib);
		*paraminfo=NULL;
	}
}

/** Function to initialize the cell parameters reading */
static ParamReader* initialize_paramdistrib(struct Parameters* par, int* status){

	// Allocate memory for the ParamDistrib structure
	ParamReader* reader=(ParamReader*)malloc(sizeof(*reader));
	CHECK_MALLOC_RET_NULL(reader);

	// If the chosen method is through a parameters image
	if (!strcmp(par->InputType, "IMAGE")) {
		ImageCellDistrib* distrib = newImageCellDistrib(par->ParamFile,par->ParamInputNames,par->FluxIndex,par->StartIndex,status);
		reader->paraminfo=distrib;
		reader->npara = distrib->num_par;
		reader->get_next_cell = &getNextCellFromImage;
		reader->get_minmax = &getMinMaxFromImage;
		reader->free_param_info = &freeImageCellDistrib;
	} else if (!strcmp(par->InputType, "TABLE")){
		TableCellDistrib* distrib = newTableCellDistrib(par->ParamFile,par->ParamInputNames,par->StartIndex,status);
		reader->paraminfo=distrib;
		reader->npara = distrib->num_par;
		reader->get_next_cell = &getNextCellFromTable;
		reader->get_minmax = &getMinMaxFromTable;
		reader->free_param_info = &freeTableCellDistrib;
	} else {
		printf("Error! Unrecognized input format!!\n");
		*status = EXIT_FAILURE;
		return(NULL);
	}

	return(reader);
}

/** Destructor of the ParamReader function */
static void freeParamReader(ParamReader* reader, int* status){
	if (reader!=NULL){
		reader->free_param_info(&(reader->paraminfo),status);
		free(reader);
		reader=NULL;
	}
}

static void addSpecCombiToCat(struct Parameters* par, par_info *data_par, img_list *li, int *status){

	SimputMIdpSpec* spec = newSimputMIdpSpec(status);

	// Set par string Xspec
	printf("Loading Parameter File in XSPEC: %s \n",par->XSPECFile);
	char *XSPECsetPar;
	get_setPar_string_xspec(&XSPECsetPar, data_par, li, status);
	printf("%s\n",XSPECsetPar);
	CHECK_STATUS_VOID(*status);

	// Write xspec file
	write_xspecSpec_file(par->Simput, par->XSPECFile, par->XSPECPrep, XSPECsetPar,
			par->Elow, par->Eup, par->nbins, par->logegrid, status);
	free(XSPECsetPar);
	CHECK_STATUS_VOID(*status);

	// Read result
	read_xspecSpec_file(par->Simput, spec ,status);
	CHECK_STATUS_VOID(*status);

	// Set the name of the spectrum
	spec->name = (char*) malloc (maxStrLenCat*sizeof(char));
	CHECK_NULL_VOID(spec->name,*status,"memory allocation failed");
	strcpy(spec->name,"");
	snprintf(spec->name,maxStrLenCat,"spec");
	for (int ii=0;ii<li->num_param;ii++){
		sprintf(spec->name,"%s_%d",spec->name,li->pval_ar[ii]);
	}
	if ((int)strlen(spec->name) > maxStrLenCat) {
		SIMPUT_ERROR("'NAME' of spectrum contains more than 64 characters");
		*status = EXIT_FAILURE;
	}
	CHECK_STATUS_VOID(*status);

	// Finally store it in the FITS table
	saveSimputMIdpSpec(spec,par->Simput,"SPECTRUM",1,status);

	// delete temp files and free memory
	char filename[SIMPUT_MAXSTR];
	sprintf(filename, "%s.qdp", par->Simput);
	remove(filename);
	freeSimputMIdpSpec(&spec);
}


static void addSpectrumExtFromList(struct Parameters* par,img_list* li,par_info* data_par,int* status){

	addSpecCombiToCat(par,data_par,li,status);

	while (li->next != NULL){
		CHECK_STATUS_BREAK(*status);
		li = li->next;
		addSpecCombiToCat(par,data_par,li,status);
	}
}

static void addSpecToCat(struct Parameters* par,struct param_input *ipar,double* par_array,int num_par,int id,int* status){
	SimputMIdpSpec* spec = newSimputMIdpSpec(status);

	// Set par string Xspec
	printf("Loading Parameter File in XSPEC: %s \n",par->XSPECFile);
	char *XSPECsetPar;
	get_setPar_string_xspec_nohist(&XSPECsetPar, ipar, par_array, num_par,status);
	printf("%s\n",XSPECsetPar);
	CHECK_STATUS_VOID(*status);

	// Write xspec file
	write_xspecSpec_file(par->Simput, par->XSPECFile, par->XSPECPrep, XSPECsetPar,
			par->Elow, par->Eup, par->nbins, par->logegrid, status);
	free(XSPECsetPar);
	CHECK_STATUS_VOID(*status);

	// Read result
	read_xspecSpec_file(par->Simput, spec ,status);
	CHECK_STATUS_VOID(*status);

	// Set the name of the spectrum
	spec->name = (char*) malloc (maxStrLenCat*sizeof(char));
	CHECK_NULL_VOID(spec->name,*status,"memory allocation failed");
	snprintf(spec->name,maxStrLenCat,"spec_%d",id);
	printf("%s\n",spec->name);

	// Finally store it in the FITS table
	saveSimputMIdpSpec(spec,par->Simput,"SPECTRUM",1,status);

	// delete temp files and free memory
	char filename[SIMPUT_MAXSTR];
	sprintf(filename, "%s.qdp", par->Simput);
	remove(filename);
	freeSimputMIdpSpec(&spec);
}


int simputmulticell_main() {
	// Program parameters.
	struct Parameters par;

	// Pointers declarations
	ParamReader* reader=NULL;
	SimputCtlg* cat=NULL;
	double* par_array=NULL;
	param_node* root=NULL;
	img_list *li = NULL;
	SimputSrc *src = NULL;

	// Error status.
	int status=EXIT_SUCCESS;

	// Register HEATOOL
	set_toolname("simputmulticell");
	set_toolversion("0.00");


	do { // Beginning of ERROR HANDLING Loop.

		// Read in tool parameters
		simputmulticell_getpar(&par,&status);
		CHECK_STATUS_BREAK(status);

		// Initialize param reading
		reader = initialize_paramdistrib(&par,&status);
		CHECK_STATUS_BREAK(status);

		// Setup param_input structure
		int num_param = reader->npara; // yes, not really necessary but more readable
		struct param_input ipar[num_param];
		set_input_par(ipar, &par, num_param, &status);
		CHECK_STATUS_BREAK(status);

		// Calibrate parameters arrays (min/max)
		par_info parinf[num_param];
		if (!par.DirectMatch) {
			cal_par_arrays(parinf,ipar,reader,num_param,&status);
			CHECK_STATUS_BREAK(status);
		}

		// Open source catalog
		check_if_output_exists(par.Simput,par.clobber,&status);
		cat=openSimputCtlg(par.Simput, READWRITE,maxStrLenCat,maxStrLenCat,maxStrLenCat,maxStrLenCat, &status);
		CHECK_STATUS_BREAK(status);

		// Populate parameter tree and source catalog
		par_array=(double*)malloc(num_param*sizeof(double));
		double ra=0;
		double dec=0;
		double flux=0.;
		double dummy_double = 0;double* dummy_double_ptr=&dummy_double; double** dummy_img=&dummy_double_ptr;
		int src_id=par.StartIndex;
		char path_fspec[SIMPUT_MAXSTR];
		char spec_name[maxStrLenCat];
		char src_name[SIMPUT_MAXSTR];

		while(reader->get_next_cell(reader->paraminfo,par_array,&ra,&dec,&flux,&status)){
			// get the IDs of the closest grid point
			// matching the given values in the parameter images
			int id_array[num_param];
			if (par.DirectMatch){
				addSpecToCat(&par,ipar,par_array,num_param,src_id,&status);
				snprintf(path_fspec,SIMPUT_MAXSTR,"[SPECTRUM,1][NAME=='spec_%d']",src_id);
			} else {
				get_par_id_from_par_array(id_array,par_array,ipar,num_param);

				// add parameter combination to parameter tree
				param_node *tmp_ptr = get_pointer_to_leave(id_array, parinf, &root, &status);
				tmp_ptr->img = dummy_img; // this allows to use the linked list structure defined for simputmultispec (hacky but prevents code rewriting)

				// add source to catalog
				// Get BASE name fomr the Param Combi (in the list)
				strcpy(spec_name,"");
				for (int ii=0;ii<num_param;ii++){
					sprintf(spec_name,"%s_%d",spec_name,id_array[ii]);
				}
				snprintf(path_fspec,SIMPUT_MAXSTR,"[SPECTRUM,1][NAME=='spec%s']",spec_name);
			}

			if ((int)strlen(spec_name) > maxStrLenCat) {
				SIMPUT_ERROR("'NAME' of spectrum contains more than 64 characters");
				status = EXIT_FAILURE;
			}
			CHECK_STATUS_BREAK(status);

			snprintf(src_name,SIMPUT_MAXSTR,"src_%d",src_id);
			// Add source
			if (!strcmp(par.InputType, "IMAGE")) flux*=1e-14; //Cannot have very small numbers in image due to SIMPUT automatic image summation
			src = newSimputSrcV(src_id,src_name, ra*M_PI/180., dec*M_PI/180.,0., 1., par.Emin, par.Emax,flux,path_fspec, "NULL", "NULL", &status);
			CHECK_STATUS_BREAK(status);
			appendSimputSrc(cat,src,&status);
			CHECK_STATUS_BREAK(status);
			freeSimputSrc(&src);// A shame to do it each time, but the way newSimputSrcV allocates memory for the strings forces me to do so (would not like to write another similar function)

			src_id++;
		}

		if (!par.DirectMatch){
			// Gather the spectrum combinations in a linked list
			li = create_img_list(root,parinf,num_param,&status);
			CHECK_STATUS_BREAK(status);

			// Compute the spectrum extension in the source catalog
			addSpectrumExtFromList(&par,li,parinf,&status);

			// Free memory depending on run
			free_param_tree(&root,parinf);
			free_par_info_array(parinf,num_param);
		}

	} while(0); // END of error handling loop.

	// Free catalog
	freeSimputCtlg(&cat, &status);
	freeParamReader(reader, &status);
	free(par_array);
	if (!par.DirectMatch){
		free_img_list(&li,0);
	}
	freeParStrings(&par);

	if (EXIT_SUCCESS==status) {
		headas_chat(3, "finished successfully!\n\n");
		return(EXIT_SUCCESS);
	}
	return(EXIT_FAILURE);
}

void simputmulticell_getpar(struct Parameters* const par, int* status) {

	query_simput_parameter_file_name("simput",&(par->Simput), status);
	query_simput_parameter_string("ISISFile", &(par->ISISFile), status );
	query_simput_parameter_string("XSPECFile", &(par->XSPECFile), status );
    query_simput_parameter_string("XSPECPrep", &(par->XSPECPrep), status );

	// *** PARAMETER information *** //
	query_simput_parameter_string("ParamFile", &(par->ParamFile), status );
	query_simput_parameter_string("ParamNames", &(par->ParamNames), status );
	query_simput_parameter_string("ParamInputNames", &(par->ParamInputNames), status );
	query_simput_parameter_int("FluxIndex", &(par->FluxIndex), status );
	query_simput_parameter_string("ParamsNumValues", &(par->ParamsNumValues), status );
	query_simput_parameter_string("ParamsLogScale", &(par->ParamsLogScale), status );
	query_simput_parameter_string("InputType", &(par->InputType), status );
	query_simput_parameter_bool("DirectMatch", &(par->DirectMatch), status );
	query_simput_parameter_long("StartIndex", &(par->StartIndex), status );

	query_simput_parameter_int("chatter", &par->chatter, status );
	query_simput_parameter_bool("clobber", &par->clobber, status );
	query_simput_parameter_bool("history", &par->history, status );
	// ***************************** //

	// *** PARAMETER information *** //
	query_simput_parameter_float("Elow", &par->Elow, status );
	query_simput_parameter_float("Eup", &par->Eup, status );
	query_simput_parameter_float("Estep", &par->Estep, status );

	query_simput_parameter_float("Emin", &par->Emin, status );
	query_simput_parameter_float("Emax", &par->Emax, status );

	query_simput_parameter_bool("logEgrid", &par->logegrid, status );
	query_simput_parameter_int("Nbins", &par->nbins, status );


	// need to check how the energy grid for the spectrum should be calculated
	if (par->Estep > 1e-6){
		SIMPUT_WARNING(" ** deprecated use of the Estep parameter ** \n    use Nbins instead to define the energy grid.");
		par->nbins = (par->Eup - par->Elow) / par->Estep;
		par->logegrid = 0;
		printf(" -> given Estep=%.4e converted to nbins=%i on a linear grid \n",par->Estep,par->nbins);
	}

}

void freeParStrings(struct Parameters* par){
	free(par->Simput);
	free(par->ISISFile);
	free(par->XSPECFile);
        free(par->XSPECPrep);
	free(par->ParamFile);
	free(par->ParamNames);
	free(par->ParamInputNames);
	free(par->ParamsNumValues);
	free(par->ParamsLogScale);
	free(par->InputType);
}
