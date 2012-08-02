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
* $Id: lin.inc,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*=======================================================================

*     Functions.
      EXTERNAL  LINCPY, LINFREE, LINGET, LINGTD, LINGTI, LININI, LINP2X,
     :          LINPRT, LINPTD, LINPTI, LINPUT, LINSET, LINX2P
      INTEGER   LINCPY, LINFREE, LINGET, LINGTD, LINGTI, LININI, LINP2X,
     :          LINPRT, LINPTD, LINPTI, LINPUT, LINSET, LINX2P

*     Length of the LINPRM data structure (INTEGER array) on 64-bit
*     machines.  Only needs to be 18 on 32-bit machines.
      INTEGER   LINLEN
      PARAMETER (LINLEN = 28)

*     Codes for LIN data structure elements used by LINPUT and LINGET.
      INTEGER   LIN_CDELT, LIN_CRPIX, LIN_FLAG, LIN_NAXIS, LIN_PC

      PARAMETER (LIN_FLAG   = 100)
      PARAMETER (LIN_NAXIS  = 101)
      PARAMETER (LIN_CRPIX  = 102)
      PARAMETER (LIN_PC     = 103)
      PARAMETER (LIN_CDELT  = 104)

*     Codes for LIN data structure elements used by LINGET (only).
      INTEGER   LIN_ERR, LIN_IMGPIX, LIN_PIXIMG, LIN_UNITY

      PARAMETER (LIN_PIXIMG = 200)
      PARAMETER (LIN_IMGPIX = 201)
      PARAMETER (LIN_UNITY  = 202)
      PARAMETER (LIN_ERR    = 203)

*     Error codes and messages.
      INTEGER   LINERR_MEMORY, LINERR_NULL_POINTER, LINERR_SINGULAR_MTX,
     :          LINERR_SUCCESS

      PARAMETER (LINERR_SUCCESS      = 0)
      PARAMETER (LINERR_NULL_POINTER = 1)
      PARAMETER (LINERR_MEMORY       = 2)
      PARAMETER (LINERR_SINGULAR_MTX = 3)

      CHARACTER LIN_ERRMSG(0:3)*80
      COMMON /LIN_DATA/ LIN_ERRMSG
