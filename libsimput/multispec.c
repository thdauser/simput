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

#include "multispec.h"

void min_max(double* array, int n, double* minn, double* maxx)  {
	*minn = array[0];
	*maxx = array[0];
   for(int i = 1; i < n; i++ )  {
	   *minn = min(array[i],*minn);
	   *maxx = max(array[i],*maxx);
	}
}

SimputImg* loadSimputImg_map(const char* const parname, int* const status) {
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

char* parse_string(char *str){
	char *pch;
	pch = strtok (str,";");
	return pch;
}

int parse_num_param(char *paramstr){

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

char** parse_string2array(char *paramstr, int n, int *status){

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

void cal_single_par_array(double pmin, double pmax, int npar, double *pvals,
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

void print_par_array(par_info par, int ii, int logScale){
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

int get_single_par_id(double pmin, double pmax, int npar, double pval, int logScale){
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

param_node *init_root_node(int *status){
	// we define root as a "normal" node, without parameter values
	return insert_node(-1,-1,NULL,status);
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

void alloc_leaves(param_node *node, int num_leaves, int *status){

	// allocate array of indices for each possible parameter value
	node->next = (param_node**) malloc( num_leaves*(sizeof(param_node*)));
	CHECK_NULL_VOID(node, *status,"memory allocation failed");
	CHECK_MALLOC_VOID(node->next);

	int ii;
	for(ii=0;ii<num_leaves;ii++){
		node->next[ii] = NULL;
	}

}


param_node *get_pointer_to_leave(int *id_array, par_info par[], param_node **root, int *status){

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


void add_img_to_list(img_list **li,  double **img, int *pind_ar, int num_param, int *status){

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

/** Frees the memory allocated in a par_info array */
void free_par_info_array(par_info par[],int num_par){
	if(par!=NULL){
		for (int i=0;i<num_par;i++){
			free(par[i].par_names);
			free(par[i].pvals);
		}
	}
}

/** Frees memory of a param tree. Note that the images are not freed as there are going to be dealt with through the
 * image linked list. */
void free_param_tree(param_node **root,par_info par[]){
	if (*root!=NULL){
		if ((*root)->next!=NULL){
			for (int i=0;i<par[(*root)->param_num+1].num_pvals;i++){
				free_param_tree(&((*root)->next[i]),par);
			}
		}
		free((*root)->next);
		free(*root);
		*root = NULL;
	}
}

img_list* create_img_list(param_node *root, par_info *par, int num_param, int *status){
	img_list* start = NULL;

	int pind_ar[num_param];

	int ii;
	for (ii=0;ii<num_param;ii++)
		pind_ar[ii] = -1;

	loop_img_for_node(root, &start,par, pind_ar, num_param, status);

	return start;

}

/** Frees an img_list. free_img param determines whether the images in there should also be freed*/
void free_img_list(img_list** li,int free_img){
	if(*li!=NULL){
		if((*li)->next!=NULL){
			free_img_list(&((*li)->next),free_img);
		}
		if (free_img){
			free((*li)->img);
		}
		free((*li)->pval_ar);
		free(*li);
		*li=NULL;
	}
}

void check_if_output_exists(char * outfile, int clobber, int * status){
	// Check if the output file already exists and remove the old
	// one if necessary.
	int exists;
	fits_file_exists(outfile, &exists, status);
	CHECK_STATUS_VOID(*status);
	if (0!=exists) {
		if (0!=clobber) {
			// Delete the file.
			remove(outfile);
		} else {
			// Throw an error.
			char msg[SIMPUT_MAXSTR];
			sprintf(msg, "file '%s' already exists", outfile);
			SIMPUT_ERROR(msg);
			return;
		}
	}
}

void get_setPar_string_xspec(char **str, par_info *data_par,img_list *li, int *status){

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

void get_setPar_string_xspec_nohist(char **str, struct param_input *ipar, double* par_array,int num_par,int* status){
	*str = (char*) malloc (SIMPUT_MAXSTR*sizeof(char));
	CHECK_NULL_VOID(*str,*status,"memory allocation failed");
	strcpy(*str,"");
	int ii;
	for (ii=0;ii<num_par;ii++){
		sprintf(*str,"%s newpar %s %.8e ;",
				*str,
				ipar[ii].param_names,
				par_array[ii]);
	}
}
