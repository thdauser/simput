*=======================================================================
*
* WCSLIB 4.13 - an implementation of the FITS WCS standard.
* Copyright (C) 1995-2012, Mark Calabretta
*
* This file is part of WCSLIB.
*
* WCSLIB is free software: you can redistribute it and/or modify it
* under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* WCSLIB is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
* License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with WCSLIB.  If not, see http://www.gnu.org/licenses.
*
* Correspondence concerning WCSLIB may be directed to:
*   Internet email: mcalabre@atnf.csiro.au
*   Postal address: Dr. Mark Calabretta
*                   Australia Telescope National Facility, CSIRO
*                   PO Box 76
*                   Epping NSW 1710
*                   AUSTRALIA
*
* Author: Mark Calabretta, Australia Telescope National Facility
* http://www.atnf.csiro.au/~mcalabre/index.html
* $Id: tab.inc,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*=======================================================================

*     Functions.
      EXTERNAL  TABCPY, TABFREE, TABGET, TABGTD, TABGTI, TABINI, TABMEM,
     :          TABPRT, TABPTD, TABPTI, TABPUT, TABS2X, TABSET, TABX2S
      INTEGER   TABCPY, TABFREE, TABGET, TABGTD, TABGTI, TABINI, TABMEM,
     :          TABPRT, TABPTD, TABPTI, TABPUT, TABS2X, TABSET, TABX2S

*     Length of the TABPRM data structure (INTEGER array) on 64-bit
*     machines.  Only needs to be 24 on 32-bit machines.
      INTEGER   TABLEN
      PARAMETER (TABLEN = 40)

*     Codes for TAB data structure elements used by TABPUT and TABGET.
      INTEGER   TAB_COORD, TAB_CRVAL, TAB_FLAG, TAB_INDEX, TAB_K, TAB_M,
     :          TAB_MAP

      PARAMETER (TAB_FLAG    = 100)
      PARAMETER (TAB_M       = 101)
      PARAMETER (TAB_K       = 102)
      PARAMETER (TAB_MAP     = 103)
      PARAMETER (TAB_CRVAL   = 104)
      PARAMETER (TAB_INDEX   = 105)
      PARAMETER (TAB_COORD   = 106)

*     Codes for TAB data structure elements used by TABGET (only).
      INTEGER   TAB_DELTA, TAB_ERR, TAB_EXTREMA, TAB_NC, TAB_P0,
     :          TAB_SENSE

      PARAMETER (TAB_NC      = 200)
      PARAMETER (TAB_SENSE   = 201)
      PARAMETER (TAB_P0      = 202)
      PARAMETER (TAB_DELTA   = 203)
      PARAMETER (TAB_EXTREMA = 204)
      PARAMETER (TAB_ERR     = 205)

*     Error codes and messages.
      INTEGER   TABERR_BAD_PARAMS, TABERR_BAD_WORLD, TABERR_BAD_X,
     :          TABERR_MEMORY, TABERR_NULL_POINTER, TABERR_SUCCESS

      PARAMETER (TABERR_SUCCESS      = 0)
      PARAMETER (TABERR_NULL_POINTER = 1)
      PARAMETER (TABERR_MEMORY       = 2)
      PARAMETER (TABERR_BAD_PARAMS   = 3)
      PARAMETER (TABERR_BAD_X        = 4)
      PARAMETER (TABERR_BAD_WORLD    = 5)

      CHARACTER TAB_ERRMSG(0:5)*80
      COMMON /TAB_DATA/ TAB_ERRMSG
