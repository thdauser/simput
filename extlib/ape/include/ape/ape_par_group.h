/** \file ape_par_group.h
    \brief Declaration of parameter group facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_par_group_h
#define ape_ape_par_group_h

#include "ape/ape_list.h"
#include "ape/ape_par.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ApeParGroup;
typedef struct ApeParGroup ApeParGroup;

ApeParGroup * ape_par_group_create(void);

void ape_par_group_destroy(ApeParGroup * ape_par_group);

ApeList * ape_par_group_get_list(ApeParGroup * group);

ApePar * ape_par_group_get_par(ApeListIterator position);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log$
*/
