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
*/

#ifndef MULTISPEC_H
#define MULTISPEC_H

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"
#include "parinput.h"

// FIXME: can we do this globally?
#ifndef min
#define min(a,b)        (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b)        (((a)>(b))?(a):(b))
#endif


#define CHECK_MALLOC_RET_NULL(a) \
  if (NULL==a) { \
    SIMPUT_ERROR("memory allocation failed"); \
    return NULL;\
  }

#define CHECK_MALLOC_VOID(a) \
  if (NULL==a) { \
    SIMPUT_ERROR("memory allocation failed"); \
    return;\
  }

const int DEFAULT_NUM_VALUES = 8;
const int maxStrLenCat = 64;

struct param_input{
	double minPar;
	double maxPar;
	int num_values;
	int logScale;
	char* param_files;
	char* param_names;

};

typedef struct{
	int num_param;

	// values for each parameter
	int num_pvals;
	char *par_names;
	double *pvals;
} par_info;


struct node{

	// number of the parameter
	int param_num;

	// value index of the parameter (pvals[ind])
	int pind;

	// pointer to the parameter structure
	par_info *par;

	// image for the certain parameter combination
	double **img;

	// pointer to the next array of elements (i.e., param_num+1!)
	struct node **next;

};
typedef struct node param_node;

struct li{

	// number of the parameter
	int num_param;

	// index array of the parameter values
	int *pval_ar;

	// image for the certain parameter combination
	double **img;

	// pointer to the next array of elements (i.e., param_num+1!)
	struct li *next;

};
typedef struct li img_list;

void min_max(double* array, int n, double* minn, double* maxx);
SimputImg* loadSimputImg_map(const char* const parname, int* const status);
char* parse_string(char *str);
int parse_num_param(char *paramstr);
char** parse_string2array(char *paramstr, int n, int *status);
void cal_single_par_array(double pmin, double pmax, int npar, double *pvals,
		int logScale);
void print_par_array(par_info par, int ii, int logScale);
int get_single_par_id(double pmin, double pmax, int npar, double pval, int logScale);
param_node *init_root_node(int *status);
param_node *insert_node(int param_num, int id, par_info *par, int *status);
void alloc_leaves(param_node *node, int num_leaves, int *status);
param_node *get_pointer_to_leave(int *id_array, par_info par[], param_node **root, int *status);
void add_img_to_list(img_list **li,  double **img, int *pind_ar, int num_param, int *status);
void loop_img_for_node(param_node *node, img_list **li, par_info *par, int *pind_ar, int num_param, int *status);
img_list* create_img_list(param_node *root, par_info *par, int num_param, int *status);
void check_if_output_exists(char * outfile, int clobber, int * status);
void get_setPar_string_xspec(char **str, par_info *data_par,img_list *li, int *status);
void get_setPar_string_xspec_nohist(char **str, struct param_input *ipar, double* par_array,int num_par,int* status);
void free_par_info_array(par_info* par,int num_par);
void free_param_tree(param_node **root,par_info par[]);
void free_img_list(img_list** li,int free_img);

#endif /* MULTISPEC_H */
