#ifndef S_TREE_H
#define S_TREE_H

#include <stdlib.h>
#include <stdio.h>

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