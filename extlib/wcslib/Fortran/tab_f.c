/*============================================================================
  WCSLIB 7.7 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2021, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: tab_f.c,v 7.7 2021/07/12 06:36:49 mcalabre Exp $
*===========================================================================*/

#include <stdio.h>
#include <string.h>

#include <wcserr.h>
#include <wcsutil.h>
#include <tab.h>

// Fortran name mangling.
#include <wcsconfig_f77.h>
#define tabput_  F77_FUNC(tabput,  TABPUT)
#define tabptd_  F77_FUNC(tabptd,  TABPTD)
#define tabpti_  F77_FUNC(tabpti,  TABPTI)
#define tabget_  F77_FUNC(tabget,  TABGET)
#define tabgtd_  F77_FUNC(tabgtd,  TABGTD)
#define tabgti_  F77_FUNC(tabgti,  TABGTI)

#define tabini_  F77_FUNC(tabini,  TABINI)
#define tabmem_  F77_FUNC(tabmem,  TABMEM)
#define tabcpy_  F77_FUNC(tabcpy,  TABCPY)
#define tabcmp_  F77_FUNC(tabcmp,  TABCMP)
#define tabfree_ F77_FUNC(tabfree, TABFREE)
#define tabsize_ F77_FUNC(tabsize, TABSIZE)
#define tabprt_  F77_FUNC(tabprt,  TABPRT)
#define tabperr_ F77_FUNC(tabperr, TABPERR)
#define tabset_  F77_FUNC(tabset,  TABSET)
#define tabx2s_  F77_FUNC(tabx2s,  TABX2S)
#define tabs2x_  F77_FUNC(tabs2x,  TABS2X)

// Must match the value set in tab.inc.
#define TAB_FLAG     100
#define TAB_M        101
#define TAB_K        102
#define TAB_MAP      103
#define TAB_CRVAL    104
#define TAB_INDEX    105
#define TAB_COORD    106

#define TAB_NC       200
#define TAB_SENSE    201
#define TAB_P0       202
#define TAB_DELTA    203
#define TAB_EXTREMA  204
#define TAB_ERR      205

//----------------------------------------------------------------------------

int tabput_(
  int *tab,
  const int *what,
  const void *value,
  const int *m,
  const int *k)

{
  int k0, m0;
  const int  *ivalp;
  const double *dvalp;
  struct tabprm *tabp;

  // Cast pointers.
  tabp  = (struct tabprm *)tab;
  ivalp = (const int *)value;
  dvalp = (const double *)value;

  // Convert 1-relative FITS axis numbers to 0-relative C array indices.
  m0 = *m - 1;
  k0 = *k - 1;

  tabp->flag = 0;

  switch (*what) {
  case TAB_FLAG:
    tabp->flag = *ivalp;
    break;
  case TAB_M:
    tabp->M = *ivalp;
    break;
  case TAB_K:
    tabp->K[m0] = *ivalp;
    break;
  case TAB_MAP:
    tabp->map[m0] = *ivalp;
    break;
  case TAB_CRVAL:
    tabp->crval[m0] = *dvalp;
    break;
  case TAB_INDEX:
    tabp->index[m0][k0] = *dvalp;
    break;
  case TAB_COORD:
    tabp->coord[m0] = *dvalp;
    break;
  default:
    return 1;
  }

  return 0;
}

int tabptd_(int *tab, const int *what, const double *value, const int *m,
  const int *k)
{
  return tabput_(tab, what, value, m, k);
}

int tabpti_(int *tab, const int *what, const int *value, const int *m,
  const int *k)
{
  return tabput_(tab, what, value, m, k);
}

//----------------------------------------------------------------------------

int tabget_(const int *tab, const int *what, void *value)

{
  unsigned int l;
  int i, k, m, n;
  int    *ivalp;
  double *dvalp;
  const int *itabp;
  const struct tabprm *tabp;

  // Cast pointers.
  tabp  = (const struct tabprm *)tab;
  ivalp = (int *)value;
  dvalp = (double *)value;

  switch (*what) {
  case TAB_FLAG:
    *ivalp = tabp->flag;
    break;
  case TAB_M:
    *ivalp = tabp->M;
    break;
  case TAB_K:
    for (m = 0; m < tabp->M; m++) {
      *(ivalp++) = tabp->K[m];
    }
    break;
  case TAB_MAP:
    for (m = 0; m < tabp->M; m++) {
      *(ivalp++) = tabp->map[m];
    }
    break;
  case TAB_CRVAL:
    for (m = 0; m < tabp->M; m++) {
      *(dvalp++) = tabp->crval[m];
    }
    break;
  case TAB_INDEX:
    for (m = 0; m < tabp->M; m++) {
      for (k = 0; k < tabp->K[m]; k++) {
        *(dvalp++) = tabp->index[m][k];
      }
    }
    break;
  case TAB_COORD:
    // Don't rely on tabprm.nc being set.
    n = tabp->M;
    for (m = 0; m < tabp->M; m++) {
      n *= tabp->K[m];
    }

    for (i = 0; i < n; i++) {
      *(dvalp++) = tabp->coord[i];
    }
    break;
  case TAB_NC:
    *ivalp = tabp->nc;
    break;
  case TAB_SENSE:
    for (m = 0; m < tabp->M; m++) {
      *(ivalp++) = tabp->sense[m];
    }
    break;
  case TAB_P0:
    for (m = 0; m < tabp->M; m++) {
      *(ivalp++) = tabp->p0[m];
    }
    break;
  case TAB_DELTA:
    for (m = 0; m < tabp->M; m++) {
      *(dvalp++) = tabp->delta[m];
    }
    break;
  case TAB_EXTREMA:
    n = 2 * tabp->M;
    for (m = 1; m < tabp->M; m++) {
      n *= tabp->K[m];
    }

    for (i = 0; i < n; i++) {
      *(dvalp++) = tabp->extrema[i];
    }
    break;
  case TAB_ERR:
    // Copy the contents of the wcserr struct.
    if (tabp->err) {
      itabp = (int *)(tabp->err);
      for (l = 0; l < ERRLEN; l++) {
        *(ivalp++) = *(itabp++);
      }
    } else {
      for (l = 0; l < ERRLEN; l++) {
        *(ivalp++) = 0;
      }
    }
    break;
  default:
    return 1;
  }

  return 0;
}

int tabgtd_(const int *tab, const int *what, double *value)
{
  return tabget_(tab, what, value);
}

int tabgti_(const int *tab, const int *what, int *value)
{
  return tabget_(tab, what, value);
}

//----------------------------------------------------------------------------

int tabini_(const int *M, const int *K, int *tab)

{
  return tabini(1, *M, K, (struct tabprm *)tab);
}

//----------------------------------------------------------------------------

int tabmem_(int *tab)

{
  return tabmem((struct tabprm *)tab);
}

//----------------------------------------------------------------------------

int tabcpy_(const int *tabsrc, int *tabdst)

{
  return tabcpy(1, (const struct tabprm *)tabsrc, (struct tabprm *)tabdst);
}

//----------------------------------------------------------------------------

int tabcmp_(
  const int *cmp,
  const double *tol,
  const int *tab1,
  const int *tab2,
  int *equal)

{
  return tabcmp(*cmp, *tol, (const struct tabprm *)tab1,
                (const struct tabprm *)tab2, equal);
}

//----------------------------------------------------------------------------

int tabfree_(int *tab)

{
  return tabfree((struct tabprm *)tab);
}

//----------------------------------------------------------------------------

int tabsize_(const int *tab, int sizes[2])

{
  return tabsize((const struct tabprm *)tab, sizes);
}

//----------------------------------------------------------------------------

int tabprt_(const int *tab)

{
  // This may or may not force the Fortran I/O buffers to be flushed.  If
  // not, try CALL FLUSH(6) before calling TABPRT in the Fortran code.
  fflush(NULL);

  return tabprt((const struct tabprm *)tab);
}

//----------------------------------------------------------------------------

// If null-terminated (using the Fortran CHAR(0) intrinsic), prefix may be of
// length less than but not exceeding 72 and trailing blanks are preserved.
// Otherwise, it must be of length 72 and trailing blanks are stripped off.

int tabperr_(int *tab, const char prefix[72])

{
  char prefix_[73];
  wcsutil_strcvt(72, '\0', 1, prefix, prefix_);

  // This may or may not force the Fortran I/O buffers to be flushed.
  // If not, try CALL FLUSH(6) before calling TABPERR in the Fortran code.
  fflush(NULL);

  return wcserr_prt(((struct tabprm *)tab)->err, prefix_);
}

//----------------------------------------------------------------------------

int tabset_(int *tab)

{
  return tabset((struct tabprm *)tab);
}

//----------------------------------------------------------------------------

int tabx2s_(
  int *tab,
  const int *ncoord,
  const int *nelem,
  const double x[],
  double world[],
  int stat[])

{
  return tabx2s((struct tabprm *)tab, *ncoord, *nelem, x, world, stat);
}

//----------------------------------------------------------------------------

int tabs2x_(
  struct tabprm* tab,
  const int *ncoord,
  const int *nelem,
  const double world[],
  double x[],
  int stat[])

{
  return tabs2x((struct tabprm *)tab, *ncoord, *nelem, world, x, stat);
}
