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

#ifndef SIMPUTMERGE_H
#define SIMPUTMERGE_H 1

#include "ape/ape_trad.h"

#include "simput.h"
#include "common.h"
#include "simput_tree.h"

#define TOOLSUB simputmerge_main
#include "headas_main.c"


struct Parameters {
  char* Infile1;
  char* Infile2;
  char* Infiles;
  char* Outfile;
  char FetchExtensions;

  char clobber;
};

typedef struct{

	char* spectrum;
	char* image;
	char* timing;

}simput_refs;

typedef struct{
	char* ref;
	char* orig_ref;
}tdata;

typedef struct{
	int type;
	char* ref;
	void* data;
}simput_data;

int simputmerge_getpar(struct Parameters* const par);
void merge_simput_files(int num_cat, SimputCtlg** incat, SimputCtlg* outcat,char* outfile, char FetchExtensions,
												long* src_id, int* status);



#endif /* SIMPUTMERGE_H */
