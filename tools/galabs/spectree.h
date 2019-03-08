/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

extern int searchlev;

struct speclist{
	char* name;
	char* fullname;
	float emin;
	float emax;
	double integral;
	struct speclist *link[2];
	long nentries;
	float *energy;
	float *pflux;
	int level;
};

struct speclist *skew (struct speclist *root);

struct speclist *split (struct speclist *root);

void nil_init();

struct speclist* return_nil();

struct speclist* speclist_newnode(char* new_name, char* fullname, float e_min, float e_max, double integral, long nsteps, float *evalues, float *pfvalues, int level);

struct speclist* insert_spec(struct speclist* root, char* new_name, char* fullname, float e_min, float e_max, double integral, long nsteps, float *evalues, float *pfvalues);

struct speclist* check_if_exists(struct speclist* list, char* name, int* status);

struct speclist* check_fullname_if_exists(struct speclist* list, char* fullname, int* status);
