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

#include "simput_tree.h"

long iterate_tree_depdth(node *n);
node *find_node(node *nd, void *data, int (*cmp)(void*, void*), int *depdth);
long get_tree_depdth(tree *tree_ptr);
node *find_elmt_best(tree *tree_ptr, void *data, int *depdth);

long iterate_tree_depdth(node *n){

  long d_right=0;
  long d_left=0;

  if(n->right!=NULL){
    d_right=iterate_tree_depdth(n->right);
    d_right++;
  }
  if(n->left!=NULL){
    d_left=iterate_tree_depdth(n->left);
    d_left++;
  }
  if(d_left>d_right){
    return d_left;
  }else{
    return d_right;
  }
}

long get_tree_depdth(tree *tree_ptr){

  long d=0;
  node *ptr=tree_ptr->treeptr;

  d=iterate_tree_depdth(ptr);

  return d;
}


node *find_node(node *nd, void *data, int (*cmp)(void*, void*), int *depdth){

  node *ptr=nd;
  int cval=cmp(data, nd->data);

  if(cval<0 && nd->left!=NULL){
    (*depdth)++;
    ptr=find_node(nd->left, data, cmp, depdth);
  }else if(cval>0 && nd->right!=NULL){
    (*depdth)++;
    ptr=find_node(nd->right, data, cmp, depdth);
  }

  return ptr;
}

node *find_elmt_best(tree *tree_ptr, void *data, int *depdth){

  node *ptr=tree_ptr->treeptr;

  if(ptr!=NULL){
    (*depdth)++;
    ptr=find_node(ptr, data, tree_ptr->cmp, depdth);
  }

  return ptr;
}

node *find_elmt(tree *tree_ptr, void *data){

  if (tree_ptr->nelem==0){
	  return NULL;
  }

  int depdth=TREE_MINDEPDTH;
  node *ptr=find_elmt_best(tree_ptr, data, &depdth);

  if(tree_ptr->cmp(ptr->data, data)){
    ptr=NULL;
  }

  return ptr;
}

node* get_node(void *data, node *parent){

  node *ptr=(node*)malloc(sizeof(node));

  if(ptr!=NULL){
    ptr->parent=parent;
    ptr->left=NULL;
    ptr->right=NULL;
    ptr->data=data;
  }

  return ptr;
}

node *add_elmt(tree *tree_ptr, void *data){

  int depdth=TREE_MINDEPDTH;
  node *ptr=find_elmt_best(tree_ptr, data, &depdth);

  if(ptr==NULL){
    tree_ptr->treeptr=get_node(data, NULL);
    ptr=tree_ptr->treeptr;
  }else{

    int cval=tree_ptr->cmp(data, ptr->data);

    if(cval<0){
      if(ptr==tree_ptr->treeptr){
	ptr->left=get_node(data, tree_ptr->treeptr);
      }else{
	ptr->left=get_node(data, ptr);
      }
      if(ptr->left!=NULL){
	ptr=ptr->left;
      }
    }else if(cval>0){
      if(ptr==tree_ptr->treeptr){
	ptr->right=get_node(data, tree_ptr->treeptr);
      }else{
	ptr->right=get_node(data, ptr);
      }
      if(ptr->right!=NULL){
	ptr=ptr->right;
      }
    }
  }

  if(ptr!=NULL){
    tree_ptr->nelem++;
    if(depdth>tree_ptr->maxdepdth){
      tree_ptr->maxdepdth=depdth;
    }
  }

  return ptr;
}

node *find_largest_node(node *elmt){

  if(elmt->right!=NULL){
    node *ptr=elmt->right;
    ptr=find_largest_node(ptr);
    return ptr;
  }else{
    return elmt;
  }
}

node *find_smallest_node(node *elmt){

  if(elmt->left!=NULL){
    node *ptr=elmt->left;
    ptr=find_smallest_node(ptr);
    return ptr;
  }else{
    return elmt;
  }
}

void remove_elmt(tree *tree_ptr, node *elmt){

  node *par=elmt->parent;
  int cval;
  if(par!=NULL){
    cval=tree_ptr->cmp(elmt->data, par->data);
  }else{
    cval=0;
  }

  if(elmt->right!=NULL && elmt->left!=NULL){
    // element has two children.
    // find largest element beneath the element
    // to be removed, as it has less than 2 children.
    node *largest=find_largest_node(elmt);

    // if it is directly linked to the element to be
    // removed, substitute it
    if(largest==elmt->right){
      // new place for second child is the smallest end
      // of the tree beneath the replacement
      node *scnd=find_smallest_node(largest);
      scnd->left=elmt->left;
      // change parent's pointer
      if(cval<0){
	par->left=largest;
      }else if (cval>0){
	par->right=largest;
      }else{
	tree_ptr->treeptr=largest;
      }
      largest->parent=par;
    }else{
      // if it is deeper below the element to be removed,
      // first replace the pointer on it with its own left
      // pointer
      largest->parent->right=largest->left;
      if(largest->left!=NULL){
	largest->left->parent=largest->parent;
      }
      // then link the left child of the element to be removed
      // to the smallest branch below elmt->right
      node *scnd=find_smallest_node(elmt->right);
      scnd->left=elmt->left;
      elmt->left->parent=scnd;

      // then unlink the largest and substitute the element to
      // be removed with it
      largest->parent=elmt->parent;
      if(cval<0){
	par->left=largest;
      }else if(cval>0){
	par->right=largest;
      }else{
	tree_ptr->treeptr=largest;
      }
      largest->left=elmt->right;
      largest->left->parent=largest;
    }
  }else{
    // element has only one or zero children,
    // so substitute it by its one child or NULL
    if(cval<0){
      if(elmt->right!=NULL){
	par->left=elmt->right;
	elmt->parent=par;
      }else if(elmt->left!=NULL){
	par->left=elmt->left;
	elmt->parent=par;
      }else{
	par->left=NULL;
      }
    }else{
      if(elmt->right!=NULL){
	par->right=elmt->right;
	elmt->parent=par;
      }else if(elmt->left!=NULL){
	par->right=elmt->left;
	elmt->parent=par;
      }else{
	par->right=NULL;
      }
    }
  }

  tree_ptr->free_elmt(elmt->data);
  free(elmt);
  tree_ptr->nelem--;
  tree_ptr->maxdepdth=get_tree_depdth(tree_ptr);

  return;
}

void free_nodes_recursively(node *ptr, void (*free_elmt)(void*)){

  if(ptr->left!=NULL){
    free_nodes_recursively(ptr->left, free_elmt);
  }
  if(ptr->right!=NULL){
    free_nodes_recursively(ptr->right, free_elmt);
  }
  free_elmt(ptr->data);
  free(ptr);

  return;
}

void free_tree(tree *tree_ptr){

  node *ptr=tree_ptr->treeptr;

  if(ptr!=NULL){
    free_nodes_recursively(ptr, tree_ptr->free_elmt);
    tree_ptr->treeptr=NULL;
  }

  tree_ptr->cmp=NULL;
  tree_ptr->free_elmt=NULL;
  tree_ptr->nelem=0;
  tree_ptr->maxdepdth=0;

  return;
}

tree *get_tree(int (*cmp)(void*, void*),
	       void (*free_elmt)(void*)){

  tree *ptr=(tree*)malloc(sizeof(tree));

  if(ptr!=NULL){
    ptr->cmp=cmp;
    ptr->free_elmt=free_elmt;
    ptr->treeptr=NULL;
    ptr->nelem=0;
    ptr->maxdepdth=0;
  }

  return ptr;
}
