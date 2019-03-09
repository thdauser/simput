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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#ifndef S_TREE_H
#define S_TREE_H

#include <stdlib.h>
#include <stdio.h>

#define TREE_MINDEPDTH 1

/** Type definitions */

/** Node structure definition */
typedef struct node{

  /** pointer to parent */
  struct node *parent;

  /** pointer to left leaf */
  struct node *left;

  /** pointer to right leaf */
  struct node *right;

  /** data pointer */
  void *data;

}node;

/** Tree structure definition */
typedef struct{

  /** Comparison function. */
  /** Must return -1, 0, 1 */
  int (*cmp)(void*, void*);

  /** Free function for data type */
  void (*free_elmt)(void*);

  /** pointer to first node */
  node *treeptr;

  /** Number of elements in tree */
  long nelem;

  /** Maximum depdth of tree */
  long maxdepdth;

}tree;


/** Function declarations */

/** Function to find the element ressembling data */
node *find_elmt(tree *tree_ptr, void *data);

/** Function to add an element with data */
node *add_elmt(tree *tree_ptr, void *data);

/** Function to delete an element from the tree */
void remove_elmt(tree *tree_ptr, node *elmt);

/** Function to free the whole tree */
void free_tree(tree *tree_ptr);

/** Function to find the largest node below a node */
node *find_largest_node(node *elmt);

/** Function to find the smallest node below a node */
node *find_smallest_node(node *elmt);

/** Function to initialize a tree */
tree *get_tree(int (*cmp)(void*, void*),
	       void (*free_elmt)(void*));

#endif /* S_TREE_H */
