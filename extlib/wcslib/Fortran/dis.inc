*=======================================================================
*
* WCSLIB 5.19 - an implementation of the FITS WCS standard.
* Copyright (C) 1995-2018, Mark Calabretta
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
* Direct correspondence concerning WCSLIB to mark@calabretta.id.au
*
* Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
* http://www.atnf.csiro.au/people/Mark.Calabretta
* $Id: dis.inc,v 5.19.1.1 2018/07/26 15:41:42 mcalabre Exp mcalabre $
*=======================================================================

*     Functions.
      EXTERNAL  DISALLOC,
     :          DISPUT, DISPTC,  DISPTD,  DISPTI,
     :          DISGET, DISGTC,  DISGTD,  DISGTI,
     :          DISNDP, DPFILL,
     :          DISINI, DISINIT, DISCPY,  DISFREE,
     :          DISPRT, DISPERR, DISHDO,
     :          DISSET, DISP2X,  DISX2P,  DISWARP

      INTEGER   DISALLOC,
     :          DISPUT, DISPTC,  DISPTD, DISPTI,
     :          DISGET, DISGTC,  DISGTD, DISGTI,
     :          DISNDP, DPFILL,
     :          DISINI, DISINIT, DISCPY, DISFREE,
     :          DISPRT, DISPERR, DISHDO,
     :          DISSET, DISP2X,  DISX2P, DISWARP

*     Length of the DISPRM and DPKEY data structures (INTEGER arrays)
*     on 64-bit machines.  DISLEN only needs to be 28 on 32-bit
*     machines.  DPLEN does not differ.
      INTEGER   DISLEN, DPLEN
      PARAMETER (DISLEN = 48)
      PARAMETER (DPLEN  = 22)

*     Codes for DIS data structure elements used by DISPUT and DISGET.
      INTEGER   DIS_FLAG, DIS_NAXIS, DIS_DTYPE, DIS_NDP, DIS_NDPMAX,
     :          DIS_DP, DIS_MAXDIS, DIS_TOTDIS

      PARAMETER (DIS_FLAG   = 100)
      PARAMETER (DIS_NAXIS  = 101)
      PARAMETER (DIS_DTYPE  = 102)
      PARAMETER (DIS_NDP    = 103)
      PARAMETER (DIS_NDPMAX = 104)
      PARAMETER (DIS_DP     = 105)
      PARAMETER (DIS_MAXDIS = 106)
      PARAMETER (DIS_TOTDIS = 107)

*     Codes for DIS data structure elements used by DISGET (only).
      INTEGER   DIS_AXMAP, DIS_NHAT, DIS_OFFSET, DIS_SCALE, DIS_IPARM,
     :          DIS_DPARM, DIS_INAXIS, DIS_NDIS
      PARAMETER (DIS_AXMAP  = 200)
      PARAMETER (DIS_NHAT   = 201)
      PARAMETER (DIS_OFFSET = 202)
      PARAMETER (DIS_SCALE  = 203)
      PARAMETER (DIS_IPARM  = 204)
      PARAMETER (DIS_DPARM  = 205)
      PARAMETER (DIS_INAXIS = 206)
      PARAMETER (DIS_NDIS   = 207)

*     Codes for DIS data structure elements used by DISGET (only).
      INTEGER   DIS_ERR

      PARAMETER (DIS_ERR    = 208)

*     Error codes and messages.
      INTEGER   DISERR_SUCCESS, DISERR_NULL_POINTER, DISERR_MEMORY,
     :          DISERR_BAD_PARAM, DISERR_DISTORT, DISERR_DEDISTORT

      PARAMETER (DISERR_SUCCESS      = 0)
      PARAMETER (DISERR_NULL_POINTER = 1)
      PARAMETER (DISERR_MEMORY       = 2)
      PARAMETER (DISERR_BAD_PARAM    = 3)
      PARAMETER (DISERR_DISTORT      = 4)
      PARAMETER (DISERR_DEDISTORT    = 5)

      CHARACTER DIS_ERRMSG(0:5)*80
      COMMON /DIS_DATA/ DIS_ERRMSG
