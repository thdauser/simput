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


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "simputmultispec.h"

// FIXME: can we do this globally?
#ifndef min
#define min(a,b)        (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b)        (((a)>(b))?(a):(b))
#endif

static void min_max(double* array, int n, double* minn, double* maxx)  {
	*minn = array[0];
	*maxx = array[0];
   for(int i = 1; i < n; i++ )  {
	   *minn = min(array[i],*minn);
	   *maxx = max(array[i],*maxx);
	}
}
static void min_max_2d(double** array, int n1, int n2, double* minn, double* maxx){
	// initialize
	min_max(array[0], n2, minn, maxx);
	double minn_tmp=0., maxx_tmp=0.;
	for(int i = 1; i < n1; i++ )  {
		min_max(array[i], n2, &minn_tmp, &maxx_tmp);
		*minn = min(*minn,minn_tmp);
		*maxx = max(*maxx,maxx_tmp);
	}
}

static SimputImg* loadSimputImg_map(const char* const parname, int* const status) {
	// ***** load the Simput Image ****
	// (need this function, as the function in simput.h only returns
	// the distribution of the map and not the map itself) Generally:
	// The data in the SimputImg still represents the
	// probability distribution stored in the image data structure.
	// However, we need the actual image,
	// NOT the distribution function. Therefore we have to inverted the
	// summing process.
    SimputImg* img = loadSimputImg(parname, status);
	long ii;
	double sum=0.;
	for (ii=0; ii<img->naxis1; ii++) {
		long jj;
		for (jj=0; jj<img->naxis2; jj++) {
			double buffer = img->dist[ii][jj];
			img->dist[ii][jj] -= sum;
			sum = buffer;
		}
	}
	return img;
}
static void saveSimputImg_map(SimputImg* const img,
		   const char* const filename,
		   char* const extname, int extver,
		   int* const status) {
	// ***** save the Simput Image ****
	// (need this function, as the function in simput.h takes
	// the distribution of the map and not the map itself)

    double sum=0.;
    long ii,jj;
    for(ii=0; ii<img->naxis1; ii++) {
      for(jj=0; jj<img->naxis2; jj++) {
    	  sum += img->dist[ii][jj];
    	  img->dist[ii][jj] = sum;
	  }
    }

    // now save it
    saveSimputImg(img, filename, extname, extver, status);

}

static double get_interp_from_img(const SimputImg* ctsImg, SimputImg* parImg,
		int indX, int indY, int* outside, int* status){
	// interpolate the parameter value for a certain pixel (indX, indY) in the cts-Image

	// 1.0 coord is exactly in the middle of pix[0]
	double xpix = (double) indX + 1.0 ;
	double ypix = (double) indY + 1.0;

	// get sky coords from the brigthness image
	double ra = 0.0, dec = 0.0;
	simput_p2s(ctsImg->wcs,  xpix, ypix, &ra, &dec, status);
	CHECK_STATUS_RET(*status, -1.0);

	// look up to which pixels these refer in the param maps
	// FIXME: need check what happens if images don't coincide
	double xpixPar = 0.0, ypixPar = 0.0;
	simput_s2p(parImg->wcs,  &xpixPar, &ypixPar, ra, dec, status);
	CHECK_STATUS_RET(*status, -1.0);

	// get the indicies in the parameter image for the given (ra,dec)
	int indXpar = lrint(xpixPar-1.0);
	int indYpar = lrint(ypixPar-1.0);

	// set to 0 if outside of the image
	double interpolParVal = 0.;
	*outside = 0;
	if (! ( indX < 0 || indX > ctsImg->naxis1-1 || indY < 0 || indY > ctsImg->naxis2-1  )) {
		double facX = (xpixPar - xpix);
		double facY = (ypixPar - ypix);

		// extrapolate if the pixel is really close to the boundary
		if (indX == ctsImg->naxis1-1) {
			facX += 1.0;
			indXpar -= 1;
		}
		if (indY == ctsImg->naxis2-1) {
			facY += 1.0;
			indYpar -= 1;
		}

		if( indXpar < 0 || indXpar > parImg->naxis1-1 || indYpar < 0 || indYpar > parImg->naxis2-1 ){
			// Throw an error.
			char msg[SIMPUT_MAXSTR];
			sprintf(msg,
					"ParamFile (%s) Index out of bounds: ix=%d not in [0,%ld] or iy=%d not in [0,%ld]! Check Image coverage or WCS coordinates",
					parImg->fileref, indXpar, parImg->naxis1-1, indYpar, parImg->naxis2-1
					);
			SIMPUT_ERROR(msg);
			*status = EXIT_FAILURE;
			return interpolParVal;
		}

		// FIXME: choose a "proper" value here
		// (conservative choice right now, to avoid missing rows if the grids match)
		const double tol = 0.001;
		if (facX <= 1.0+tol && facY <= 1.0+tol){
		interpolParVal = parImg->dist[indXpar][indYpar]*(1.-facX)*(1.-facY) +
				parImg->dist[indXpar+1][indYpar]*facX*(1.-facY) +
				parImg->dist[indXpar][indYpar+1]*(1.-facX)*facY +
				parImg->dist[indXpar+1][indYpar+1]*facX*facY;
		} else {
			*outside = 1;
			interpolParVal = 0.0; // need a proper values here
		}

	}

	return interpolParVal;
}

static void cal_single_par_array(double pmin, double pmax, int npar, double *pvals,
		int logScale){

	// FIXME: what do we do for negative parameter values??
	if(logScale){
		pmin = log(pmin);
		pmax = log(pmax);
	}
	// TODO: why does this not work?
	//pvals = (double *)malloc((npar)*sizeof(double));
	//CHECK_NULL_VOID(pvals, *status,"memory allocation failed");

	for (int ii=0; ii<npar; ii++){
		pvals[ii] = ii/(npar-1.0)*(pmax-pmin) + pmin;
		if (logScale){
			pvals[ii] = exp(pvals[ii]);
		}
	}

}

static int get_single_par_id(double pmin, double pmax, int npar, double pval, int logScale){
	// this should be (much) faster than binary search

	// FIXME: what do we do for negative parameter values??
	if(logScale){
		pmin = log(pmin);
		pmax = log(pmax);
		pval = log(pval);
	}

	double id_val = (pval-pmin) / (pmax - pmin) * (npar-1);

	// we round to the nearest integer, as we want the closest grid point
	return (int) lrint(id_val);
}

static void print_par_array(par_info par, int ii, int logScale){
	int jj;
	printf("\nValues for Parameter %i (%s): \n",ii+1,par.par_names);
	for (jj=0; jj < par.num_pvals; jj++){
	    printf("  %.3e",par.pvals[jj]);
	}
	if (logScale) {
	    printf("   (log scale)");
	}
	printf("\n\n");
}

static void cal_par_arrays(par_info *par, struct param_input *ipar, SimputImg** img, int num_param, int* status){

	int ii;
	for (ii=0; ii<num_param; ii++){

	    min_max_2d(img[ii]->dist,img[ii]->naxis1,img[ii]->naxis2,
	            &ipar[ii].minPar,&ipar[ii].maxPar);

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

}


static char* parse_string(char *str){
	char *pch;
	pch = strtok (str,";");
	return pch;
}
static int parse_num_param(char *paramstr){

	int num_param = 0;

	char buffer[SIMPUT_MAXSTR];
	strcpy(buffer,paramstr); // make sure we do not destroy the string

	char *pch;
	pch = parse_string(buffer);
	while (pch != NULL) {
		pch = parse_string(NULL);
		num_param++;
	}
	return num_param;
}
static char** parse_string2array(char *paramstr, int n, int *status){

	char str[SIMPUT_MAXSTR];
	strcpy(str,paramstr); // make sure we do not destroy the string

	// simply return if the desired number of parameters is not given
	// TODO: Test case if only 1 parameter is given
	if ( n != parse_num_param(str)){
		return NULL;
	}

	// allocate array
	char **strarray = NULL;
	strarray = (char**) malloc( n * sizeof(char*) );
	CHECK_NULL_RET(strarray,*status,"malloc failed",NULL);
	int ii = 0;
	for (ii=0; ii < n; ii++){
		strarray[ii] = (char*) malloc (SIMPUT_MAXSTR * sizeof(char));
		CHECK_NULL_RET(strarray[ii],*status,"malloc failed",NULL);
	}

	char *pch = NULL;
	pch = parse_string(str);
	ii=0;
	strcpy(strarray[ii],pch);

	while ( pch != NULL){
		ii++;
		pch = parse_string(NULL);

		if (pch != NULL){
			strcpy(strarray[ii],pch);
		}
	}
	return strarray;
}
static void set_input_par(struct param_input *ipar, struct Parameters* par, int num_par, int *status) {

   	char **pfiles = NULL;
    char **pnames = NULL;
    char **pnumval = NULL;
    char **plog = NULL;

    int ii;

    pfiles = parse_string2array(par->ParamFiles,num_par,status);
    CHECK_NULL_VOID(pfiles,*status,"Error when parsing ParamFiles");

    pnames = parse_string2array(par->ParamNames,num_par,status);
    CHECK_NULL_VOID(pnames,*status,"Error when parsing ParamNames");

    for (ii=0; ii < num_par ; ii++ ){
        	ipar[ii].param_files = pfiles[ii];
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

}

static void get_par_id(int id_array[], const SimputImg *ctsImg, SimputImg *parImgs[],
		const struct param_input ipar[], const int num_param, const int ii, const int jj,
		int *status){

	int kk;
	int outside = 0;
	double interp_value[num_param];
	for ( kk=0; kk < num_param; kk++) {
		interp_value[kk] = get_interp_from_img(ctsImg, parImgs[kk], ii, jj, &outside, status);

		// mark if one parameter value is outside the image
		// (pixel value has zero brightness by default)
		if (outside){
			id_array[kk] = -1;
		} else {
			// get the parameter ID
			id_array[kk] = get_single_par_id(ipar[kk].minPar,ipar[kk].maxPar,ipar[kk].num_values,
					interp_value[kk],ipar[kk].logScale);
		}
	}

} // end parameter loop


static void alloc_leaves(param_node *node, int num_leaves, int *status){

	// allocate array of indices for each possible parameter value
	node->next = (param_node**) malloc( num_leaves*(sizeof(param_node*)));
	CHECK_NULL_VOID(node, *status,"memory allocation failed");
	CHECK_MALLOC_VOID(node->next);

	int ii;
	for(ii=0;ii<num_leaves;ii++){
		node->next[ii] = NULL;
	}

}

param_node *insert_node(int param_num, int id, par_info *par, int *status){

	// alloc and init a node

	param_node *node = (param_node*) malloc( sizeof(param_node));
	CHECK_NULL_RET(node, *status,"memory allocation failed",NULL);

	node->param_num = param_num; // level of the tree (root=-1)

	node->pind = id; // param index for a certain depth

	node->par = par;

	node->img = NULL;

	node->next = NULL;

	return node;
}

static void add_value_to_param_img(param_node *node, SimputImg* img, int indX, int indY, int *status){

	// initialize if it does not exist
	if (node->img == NULL){
		int nx = img->naxis1;
		int ny = img->naxis2;

		// alloc image
		node->img = (double **) malloc ( nx*sizeof(double*));
		CHECK_NULL_VOID(node->img, *status,"memory allocation failed");

		// init image with zeros
		for (int ii=0; ii < nx; ii++){
			node->img[ii] = (double *) malloc ( ny*sizeof(double));
			CHECK_NULL_VOID(node->img[ii], *status,"memory allocation failed");
			for (int jj=0; jj < ny; jj++){
				node->img[ii][jj] = 0.0;
			}
		}
	}
	// set the value from the brightness map
	node->img[indX][indY] = img->dist[indX][indY];
	// add image also to the list of images
}

static param_node *init_root_node(int *status){
	// we define root as a "normal" node, without parameter values
	return insert_node(-1,-1,NULL,status);
}

static param_node *get_pointer_to_leave(int *id_array, par_info par[], param_node **root, int *status){

	if (*root == NULL)
		*root = init_root_node(status);

	// start out from the root node
	param_node *ptr=*root;

	// loop through all stages of the tree
	// (using a loop here as we always need to go to the leaves of the tree
	int ii=0;
	while( ii < par[0].num_param){

		// get the parameter id
		int id = id_array[ii];

		// check if already leaves exist and create them if necessary
		if (ptr->next == NULL){
			alloc_leaves(ptr, par[ii].num_pvals, status);
		}

		// does the node for the par[ii].ind[id] already exist?
		if ( ptr->next[id] == NULL){
			// insert a new node (this is the node for par[ii].pval[id]
			ptr->next[id] = insert_node(ii, id, &(par[ii]), status);
		}

		// as we can now be sure that the node exists, we move to the next node
		ptr = ptr->next[id];
		ii++;
	}

	if (*status!=EXIT_SUCCESS)
		return NULL;
	else
		return ptr;

}

static void add_img_to_list(img_list **li,  double **img, int *pind_ar, int num_param, int *status){

	img_list *newElem;
	img_list *start = *li;

	newElem = (img_list *) malloc(sizeof(img_list));
	CHECK_NULL_VOID(newElem, *status,"memory allocation failed");
	newElem->img = img;
	newElem->next=NULL;
	newElem->num_param=num_param;
	newElem->pval_ar = (int *) malloc(num_param * sizeof(int));
	for (int ii=0;ii<num_param;ii++)
		newElem->pval_ar[ii] = pind_ar[ii];
	CHECK_NULL_VOID(newElem->pval_ar, *status,"memory allocation failed");

	// does the list already exist?
	if (start == NULL){
		*li = newElem;
	} else {
		while (start->next != NULL)
			start = start->next;
		start->next = newElem;
	}

}

void loop_img_for_node(param_node *node, img_list **li, par_info *par, int *pind_ar, int num_param, int *status){

	// CHECK IF WE ARE ON A LEAF: this defines a leaf
	if (node->next == NULL){
		if (node->img != NULL){
			add_img_to_list(li, node->img,pind_ar,num_param,status);
			return;
		} else {
			printf("Error: Dead End in the tree\n");
			*status = EXIT_FAILURE;
			return;
		}
	}

	// NOT ON A LEAF -> GO TO THE NEXT LEVEL
	int param_num = node->param_num+1; // works also with root, as there param_num=-1 !
	if (param_num >= num_param ){
		printf("Error! Number of Parameters exceeded in Loop. Should never happen!!\n");
		*status = EXIT_FAILURE;
		return;
	}

	int num_pvals_next = par[param_num].num_pvals;

	int ii;
	for (ii=0; ii < num_pvals_next; ii++){
		if ( node->next[ii] != NULL){
			// next node exists, so set its value
			pind_ar[param_num] = node->next[ii]->pind;
			// and start a loop for the next node
			loop_img_for_node(node->next[ii],li,par,pind_ar,num_param,status);
		}
	}

}

static img_list* create_img_list(param_node *root, par_info *par, int num_param, int *status){
	img_list* start = NULL;

	int pind_ar[num_param];

	int ii;
	for (ii=0;ii<num_param;ii++)
		pind_ar[ii] = -1;

	loop_img_for_node(root, &start,par, pind_ar, num_param, status);

	return start;

}

static void init_string(char **str,char *f0, char *fmt, int *status){
	int len = maxStrLenCat;
	*str=(char*) malloc( len*sizeof(char) );
	CHECK_NULL_VOID(*str, *status,"memory allocation failed");
	sprintf(*str,fmt,f0);
}

static void get_filenames(char **fsrc, char **fspec, char **fimg, img_list *li,int *status){

	char fname[maxStrLenCat];

	strcpy(fname,"");
    // Get BASE name fomr the Param Combi (in the list)
	for (int ii=0;ii<li->num_param;ii++){
		sprintf(fname,"%s_%i",fname,li->pval_ar[ii]);
	}

	init_string(fsrc,fname,"src%s",status);
	init_string(fspec,fname,"spe%s",status);
	init_string(fimg,fname,"img%s",status);

}


static void ms_save_img(SimputImg* ctsImg, double ** img, struct Parameters par, char *fimg, int *status){


	// write the Image to the dummy simputImg file
	// (with correct naxis1/2 and wcs parameters)
	// todo: add a check?
	ctsImg->dist = img;

	// todo: should we use a "better" value here?
	int extver = 1;
	saveSimputImg_map(ctsImg, par.Simput, fimg, extver, status);

}

static void ms_save_src(SimputCtlg *ctlg, char *fname, char *fspec, char *fimg,
		struct Parameters par, int src_id, float sflux, int *status){

    char path_fspec[SIMPUT_MAXSTR];
    char path_fimg[SIMPUT_MAXSTR];

    sprintf(path_fspec,"[%s,%i][NAME=='%s']",extnameSpec,extver,fspec);
    sprintf(path_fimg,"[%s,%i]",fimg,extver);

    // Get a new source entry.
    SimputSrc *src = newSimputSrcV(src_id,fname, par.RA*M_PI/180., par.Dec*M_PI/180.,
		      0., 1., par.Emin, par.Emax, sflux,
		      path_fspec, path_fimg, "NULL", status);
    CHECK_STATUS_VOID(*status);


    // add source
    appendSimputSrc(ctlg,src, status);

    free(src);
}

static void get_setPar_string_isis(char **str, par_info *data_par,img_list *li, int *status){

	*str = (char*) malloc (SIMPUT_MAXSTR*sizeof(char));
	CHECK_NULL_VOID(*str,*status,"memory allocation failed");
	strcpy(*str,"");
	int ii;
	for (ii=0;ii<li->num_param;ii++){
		sprintf(*str,"%s set_par(\"%s\",%.8e);",
				*str,
				data_par[ii].par_names,
				data_par[ii].pvals[li->pval_ar[ii]]
			    );
	}
}

static void get_setPar_string_xspec(char **str, par_info *data_par,img_list *li, int *status){

	*str = (char*) malloc (SIMPUT_MAXSTR*sizeof(char));
	CHECK_NULL_VOID(*str,*status,"memory allocation failed");
	strcpy(*str,"");
	int ii;
	for (ii=0;ii<li->num_param;ii++){
		sprintf(*str,"%s newpar %s %.8e ;",
				*str,
				data_par[ii].par_names,
				data_par[ii].pvals[li->pval_ar[ii]]
			    );
	}
}


static void delete_tmp_spec_files(struct Parameters par){
	if (strlen(par.ISISFile)>0) {
		char filename[SIMPUT_MAXSTR];
		sprintf(filename, "%s.spec0", par.Simput);
		remove(filename);
	}
	if (strlen(par.XSPECFile)>0) {
		remove(par.TmpSpecFile);
	}
}

static void ms_save_spec_isis(SimputMIdpSpec *spec, struct Parameters par,
		par_info *data_par, img_list *li, int *status){


	char *ISISsetPar;
	get_setPar_string_isis(&ISISsetPar, data_par, li, status);
	printf("%s\n",ISISsetPar);
	write_isisSpec_fits_file(par.TmpSpecFile, par.ISISFile, par.ISISPrep,ISISsetPar,
			par.Elow, par.Eup, par.nbins, par.logegrid, 0.0, 0.0, 0.0, 0.0, 0.0, status);
	free(ISISsetPar);

	CHECK_STATUS_VOID(*status);

	read_isisSpec_fits_file(par.TmpSpecFile, spec, par.ISISFile, par.Emin, par.Emax,
			0.0,0.0,0.0,0.0,status);
	CHECK_STATUS_VOID(*status);

}

static void ms_save_spec_xspec(SimputMIdpSpec *spec, struct Parameters par,
		par_info *data_par, img_list *li, int *status){


	char *XSPECsetPar;
	get_setPar_string_xspec(&XSPECsetPar, data_par, li, status);
	printf("%s\n",XSPECsetPar);

	write_xspecSpec_file(par.TmpSpecFile, par.XSPECFile, par.XSPECPrep, XSPECsetPar,
			par.Elow, par.Eup, par.nbins, par.logegrid, status);


	free(XSPECsetPar);

	CHECK_STATUS_VOID(*status);

	read_xspecSpec_file(par.TmpSpecFile, spec ,status);
	CHECK_STATUS_VOID(*status);

}


static void ms_save_spec(struct Parameters par, char *fspec,
		par_info *data_par, img_list *li, int *status){



	// --- write and then read the spectral file with isis or XSPEC ---

	SimputMIdpSpec* spec = newSimputMIdpSpec(status);

	if (strlen(par.ISISFile)>0) {

		printf("Loading Parameter File in ISIS: %s \n",par.ISISFile);
		ms_save_spec_isis(spec,par, data_par,li,status);

	} else if (strlen(par.XSPECFile) > 0) {

		printf("Loading Parameter File in XSPEC: %s \n",par.XSPECFile);
		ms_save_spec_xspec(spec,par, data_par,li,status);

	} else {
		SIMPUT_ERROR("Error. No Parameter File exits. Should not happen.");
		*status = EXIT_FAILURE;
	}
	CHECK_STATUS_VOID(*status);

	// --------------------------------------------------------- //

	// set the name of the spectrum
	spec->name = (char*) malloc (maxStrLenCat*sizeof(char));
	CHECK_NULL_VOID(spec->name,*status,"memory allocation failed");
	strcpy(spec->name,fspec);

	// finally store it in the FITS table --------------------- //
	saveSimputMIdpSpec(spec,par.Simput,extnameSpec,extver,status);
	// and get the model flux
//	double srcFlux = getSimputMIdpSpecBandFlux(spec, par.Emin,par.Emax);

	// delete temp files and free memory
	delete_tmp_spec_files(par);
	free(spec);

	return;
}


static float get_img_cts(double ** img, long n1, long n2){
	double sumCts = 0.0;

	int ii,jj;

	for (ii=0;ii<n1;ii++){
		for (jj=0;jj<n2;jj++){
			sumCts += img[ii][jj];
		}
	}

	return (float) sumCts;
}


static void save_single_combi(img_list *li, SimputImg *ctsImg, struct Parameters par,
		par_info *data_par, SimputCtlg *ctlg, int src_counter, float totalFlux, int *status){

	char *fsrc,*fspec,*fimg;
	get_filenames(&fsrc,&fspec,&fimg, li, status);
	CHECK_STATUS_VOID(*status);

    // --- SPECTRUM --- //
	ms_save_spec(par, fspec, data_par, li, status);
	CHECK_STATUS_VOID(*status);

	// --- SOURCE --- ///
	float singleSrcFlux = get_img_cts(li->img,ctsImg->naxis1,ctsImg->naxis2)/totalFlux*par.srcFlux;
	ms_save_src(ctlg, fsrc, fspec, fimg, par, src_counter,singleSrcFlux,status);
	CHECK_STATUS_VOID(*status);


	// --- IMAGE --- //
	ms_save_img(ctsImg, li->img, par, fimg, status);
	CHECK_STATUS_VOID(*status);

}

static float get_total_img_flux(struct Parameters par, int *status){
	SimputImg* img = loadSimputImg_map(par.ImageFile, status);
	return get_img_cts(img->dist,img->naxis1,img->naxis2);
}

static void loop_images_fits(img_list *li, SimputImg *ctsImg, struct Parameters par,
		par_info *data_par, SimputCtlg *ctlg, int *status){

	float totalFlux = get_total_img_flux(par, status);

	int src_counter=1;
	save_single_combi(li,ctsImg,par,data_par,ctlg,src_counter,totalFlux,status);

	src_counter++;
	while (li->next != NULL){
		CHECK_STATUS_BREAK(*status);
		li = li->next;
		save_single_combi(li,ctsImg,par,data_par,ctlg,src_counter,totalFlux,status);
		src_counter++;
	}
}

static void check_if_output_exists(struct Parameters par, int * status){
	// Check if the SIMPUT file already exists and remove the old
	// one if necessary.
	int exists = 0;
  if (par.Simput)	fits_file_exists(par.Simput, &exists, status);
	CHECK_STATUS_VOID(*status);
	if (0!=exists) {
		if (0!=par.clobber) {
			// Delete the file.
			remove(par.Simput);
		} else {
			// Throw an error.
			char msg[SIMPUT_MAXSTR];
			sprintf(msg, "file '%s' already exists", par.Simput);
			SIMPUT_ERROR(msg);
			return;
		}
	}
}

static void create_simputmultispec_file(img_list *li, SimputImg *ctsImg,
			struct Parameters par, par_info *data_par, int *status){

	// check if file exists and/or if it should be deleted
	check_if_output_exists(par,status);
	CHECK_STATUS_VOID(*status);

    // Create a new SIMPUT catalog.
	SimputCtlg* cat=openSimputCtlg(par.Simput, READWRITE,
			maxStrLenCat,maxStrLenCat,maxStrLenCat,maxStrLenCat, status);
    CHECK_STATUS_VOID(*status);

    // loop over all images to create cataloge, spec and image entries in the FITS file
	loop_images_fits(li,ctsImg,par,data_par,cat,status);
}



int simputmultispec_main() {
  // Program parameters.
  struct Parameters input_par;

  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("simputmultispec");
  set_toolversion("0.00");


  do { // Beginning of ERROR HANDLING Loop.

   // Read the parameters using PIL.
    status=simputmultispec_getpar(&input_par);
    CHECK_STATUS_BREAK(status);
    // ---- END of Initialization ----

    // ************************************************************ //
    // TODO: needs to be given als input and done automatically!!
    //		 (and in a function)
    int num_param = parse_num_param(input_par.ParamFiles);
    struct param_input ipar[num_param];
    set_input_par(ipar, &input_par, num_param, &status);
    CHECK_STATUS_BREAK(status);

    // ************************************************************ //

    SimputImg* ctsImg = loadSimputImg_map(input_par.ImageFile, &status);
    SimputImg* parImgs[num_param];

    int kk;
    for ( kk=0; kk < num_param; kk++ ){
    	parImgs[kk] = loadSimputImg_map(ipar[kk].param_files, &status);
    }	// init tree root
    CHECK_STATUS_BREAK(status);

    // load all the parameter from the input (after loading the images to calculate
    // minPar and maxPar automatically from the images)
    par_info par[num_param];
    cal_par_arrays(par,ipar,parImgs,num_param,&status);
	CHECK_STATUS_BREAK(status);

	param_node *root = NULL;

	// loop over Pixels
	int ii, jj;
	for (ii = 0; ii < ctsImg->naxis1 ; ii++){
		for (jj = 0; jj <  ctsImg->naxis2 ; jj++){
			int outside=0;
			int kk;

			// get the IDs of the closest grid point
			// matching the given values in the parameter images
			int id_array[num_param];
			get_par_id(id_array, ctsImg, parImgs, ipar, num_param, ii, jj, &status);

			for(kk=0; kk<num_param; ++kk) {
				outside |= (id_array[kk]<0);
			}

			if(!outside) {
				param_node *tmp_ptr = get_pointer_to_leave(id_array, par, &root, &status);
				add_value_to_param_img(tmp_ptr, ctsImg, ii, jj, &status);
			}

			CHECK_STATUS_BREAK(status);
		}
		CHECK_STATUS_BREAK(status);
	} // END PIXEL LOOP
	CHECK_STATUS_BREAK(status);

    // ---- now gather the images in a linked list ---- //
    img_list *li = create_img_list(root,par,num_param,&status);
    CHECK_STATUS_BREAK(status);


    // ---- write the outfile with SIMPUT ---- //

    // check if file exists and/or if it should be deleted
    create_simputmultispec_file(li,ctsImg,input_par,par,&status);
    CHECK_STATUS_BREAK(status);


    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}



int simputmultispec_getpar(struct Parameters* const par) {

  // Error status.
  int status=EXIT_SUCCESS;

  query_simput_parameter_file_name("simput",&(par->Simput), &status);
  query_simput_parameter_string("ISISFile", &(par->ISISFile), &status );
  query_simput_parameter_string("XSPECFile", &(par->XSPECFile), &status );
  query_simput_parameter_string("XSPECPrep", &(par->XSPECPrep), &status );
  query_simput_parameter_string("ISISPrep", &(par->ISISPrep), &status );
  query_simput_parameter_file_name("ImageFile", &(par->ImageFile), &status );

  // Initialize name for temporary spectrum file.
  if (strlen(par->XSPECFile) > 0) {
    create_unique_tmp_name(par->TmpSpecFile, &status);
    strncat(par->TmpSpecFile, ".qdp", 4);
  }
  if (strlen(par->ISISFile) > 0) {
    create_unique_tmp_name(par->TmpSpecFile, &status);
    strncat(par->TmpSpecFile, ".sl", 3);
  }

  if ((strlen(par->ISISFile) == 0) && (strlen(par->XSPECFile) == 0) ) {
    SIMPUT_ERROR("Need to give either an ISIS or XSPEC parameter file");
    return EXIT_FAILURE;
  }
  

  query_simput_parameter_float("RA", &par->RA, &status );
  query_simput_parameter_float("Dec", &par->Dec, &status );

  query_simput_parameter_float("srcFlux", &par->srcFlux, &status );

  // *** PARAMETER information *** //
  query_simput_parameter_string("ParamFiles", &(par->ParamFiles), &status );
  query_simput_parameter_string("ParamNames", &(par->ParamNames), &status );
  query_simput_parameter_string("ParamsNumValues", &(par->ParamsNumValues), &status );
  query_simput_parameter_string("ParamsLogScale", &(par->ParamsLogScale), &status );

  query_simput_parameter_int("chatter", &par->chatter, &status );
  query_simput_parameter_bool("clobber", &par->clobber, &status );
  query_simput_parameter_bool("history", &par->history, &status );
  // ***************************** //

  // *** PARAMETER information *** //
  query_simput_parameter_float("Elow", &par->Elow, &status );
  query_simput_parameter_float("Eup", &par->Eup, &status );
  query_simput_parameter_float("Estep", &par->Estep, &status );

  query_simput_parameter_float("Emin", &par->Emin, &status );
  query_simput_parameter_float("Emax", &par->Emax, &status );

  query_simput_parameter_bool("logEgrid", &par->logegrid, &status );
  query_simput_parameter_int("Nbins", &par->nbins, &status );

    // ***************************** //

  // need to check how the energy grid for the spectrum should be calculated
  if (par->Estep > 1e-6){
  	SIMPUT_WARNING(" ** deprecated use of the Estep parameter ** \n    use Nbins instead to define the energy grid.");
  	par->nbins = (par->Eup - par->Elow) / par->Estep;
  	par->logegrid = 0;
  	printf(" -> given Estep=%.4e converted to nbins=%i on a linear grid \n",par->Estep,par->nbins);
  }


  return(status);
}
