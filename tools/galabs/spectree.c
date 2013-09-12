#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spectree.h"

#define MAXSTRING 2048

static struct speclist *nil=NULL;

struct speclist *skew (struct speclist *root){

	if((root)->level!=0){
        	if((root)->link[0]->level == (root)->level){
			struct speclist *help = root;
			(root)=(root)->link[0];
			help->link[0]=(root)->link[1];
			(root)->link[1]=help;
		}
		(root)->link[1]=skew(((root)->link[1]));
	}
	return root;
}


struct speclist *split (struct speclist *root){

	if((root)->link[1]->link[1]->level == (root)->level &&(root)->level != 0){
		struct speclist *help=root;
		root=(root)->link[1];
		help->link[1]=(root)->link[0];
		(root)->link[0]=help;
		++((root)->level);
		(root)->link[1]=split(((root)->link[1]));
	}
	return root;
}


void nil_init(){

	puts("initialize nil..");
	nil=(struct speclist*)malloc(sizeof *nil);
	
	if(nil==NULL){
		puts("Memory allocation for nil failed. Exit.");
		exit(1);
	}

	nil->level=0;
	nil->link[0]=nil;
	nil->link[1]=nil;
}

struct speclist* return_nil(){

	return nil;
}

struct speclist* speclist_newnode(char* new_name, char* fullname, float e_min, float e_max, double integral, long nsteps, float *evalues, float *pfvalues, int level){
//Function allocates memory for a new node, sets the values and returns the pointer to the new node

	struct speclist *newnode=NULL;
	newnode=(struct speclist*)malloc(sizeof(struct speclist));
	if(newnode==NULL){
		puts("ERROR: Allocation of memory for new node failed. Exit.");
		exit(1);
	}

	newnode->name=(char*)malloc(MAXSTRING*sizeof(char));
	newnode->fullname=(char*)malloc(MAXSTRING*sizeof(char));

	if(newnode->name==NULL || newnode->fullname==NULL){
		puts("ERROR: Allocation of memory for new node names failed. Exit.");
		exit(1);
	}

	strcpy(newnode->name, new_name);
	strcpy(newnode->fullname, fullname);

	newnode->emin=e_min;
	newnode->emax=e_max;
	newnode->integral=integral;
	newnode->nentries=nsteps;
	
	newnode->energy=(float*)malloc(nsteps*sizeof(float));
	newnode->pflux=(float*)malloc(nsteps*sizeof(float));

	if(newnode->energy==NULL || newnode->pflux==NULL){
		puts("ERROR: Allocation of memory for new node spectrum arrays failed. Exit.");
		exit(1);
	}

	int ii;
	for(ii=0; ii<nsteps; ii++){
		newnode->energy[ii]=evalues[ii];
		newnode->pflux[ii]=pfvalues[ii];
	}

	newnode->link[0]=nil;
	newnode->link[1]=nil;
	newnode->level = level;

	return newnode;
}


struct speclist* insert_spec(struct speclist* root, char* new_name, char* fullname, float e_min, float e_max, double integral, long nsteps, float *evalues, float *pfvalues){
//This function inserts a new spectrum into the speclist pointed by root

	if(root==nil){
		root=speclist_newnode(new_name, fullname, e_min, e_max, integral, nsteps, evalues, pfvalues, 1);
	}else{
		int comp=strcmp(root->name, new_name);
		if(comp<0){
			comp=0;
		}else{
			comp=1;
		}
		root->link[comp]=insert_spec(root->link[comp], new_name, fullname, e_min, e_max, integral, nsteps, evalues, pfvalues);
		root=skew(root);
		root=split(root);
	}
	return root;
}

struct speclist* check_if_exists(struct speclist* list, char* name, int* status){
//This function returns a pointer to the speclist entry 
//with the searched name and sets status to 1 if its found,
//otherwise status is 0 and a NULL pointer is returned

	if(list==nil){
		*status=0;
		return NULL;
	}else{
		struct speclist* active = list;
		int comp;
		while(active!=nil){
			comp=strcmp(active->name, name);
			if(comp==0){
				*status=1;
				return active;
			}else if(comp<0){
				searchlev++;
				active=active->link[0];
			}else{
				searchlev++;
				active=active->link[1];
			}
		}
	}
	*status=0;
	return NULL;
}

struct speclist* check_fullname_if_exists(struct speclist* list, char* fullname, int* status){
//This function returns a pointer to the speclist entry 
//with the searched name and sets status to 1 if its found,
//otherwise status is 0 and a NULL pointer is returned

	if(list==nil){
		*status=0;
		return NULL;
	}else{
		struct speclist* active = list;
		int comp;
		while(active!=nil){
			comp=strcmp(active->fullname, fullname);
			if(comp==0){
				*status=1;
				return active;
			}else if(comp<0){
				searchlev++;
				active=active->link[0];
			}else{
				searchlev++;
				active=active->link[1];
			}
		}
	}
	*status=0;
	return NULL;
}
