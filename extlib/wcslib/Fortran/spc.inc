*=======================================================================
*
* WCSLIB 4.25 - an implementation of the FITS WCS standard.
* Copyright (C) 1995-2015, Mark Calabretta
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
* $Id: spc.inc,v 4.25.1.2 2015/01/06 01:02:17 mcalabre Exp mcalabre $
*=======================================================================

*     Functions.
      EXTERNAL  SPCAIPS, SPCFREE, SPCGET, SPCGTC, SPCGTD, SPCGTI,
     :          SPCINI, SPCPRT, SPCPTC, SPCPTD, SPCPTI, SPCPUT, SPCS2X,
     :          SPCSET, SPCSPXE, SPCTRNE, SPCTYPE, SPCX2S, SPCXPSE
      INTEGER   SPCAIPS, SPCFREE, SPCGET, SPCGTC, SPCGTD, SPCGTI,
     :          SPCINI, SPCPRT, SPCPTC, SPCPTD, SPCPTI, SPCPUT, SPCS2X,
     :          SPCSET, SPCSPXE, SPCTRNE, SPCTYPE, SPCX2S, SPCXPSE

*     Deprecated functions.
      EXTERNAL  SPCSPX, SPCTRN, SPCTYP, SPCXPS
      INTEGER   SPCSPX, SPCTRN, SPCTYP, SPCXPS

*     Length of the SPCPRM data structure (INTEGER array) on 64-bit
*     machines.  Only needs to be 44 on 32-bit machines.
      INTEGER   SPCLEN
      PARAMETER (SPCLEN = 50)

*     Codes for SPC data structure elements used by SPCPUT and SPCGET.
      INTEGER   SPC_CODE, SPC_CRVAL, SPC_FLAG, SPC_PV, SPC_RESTFRQ,
     :          SPC_RESTWAV, SPC_TYPE

      PARAMETER (SPC_FLAG    = 100)
      PARAMETER (SPC_TYPE    = 101)
      PARAMETER (SPC_CODE    = 102)
      PARAMETER (SPC_CRVAL   = 103)
      PARAMETER (SPC_RESTFRQ = 104)
      PARAMETER (SPC_RESTWAV = 105)
      PARAMETER (SPC_PV      = 106)

*     Codes for SPC data structure elements used by SPCGET (only).
      INTEGER   SPC_ERR, SPC_ISGRISM, SPC_W

      PARAMETER (SPC_W       = 200)
      PARAMETER (SPC_ISGRISM = 201)
      PARAMETER (SPC_ERR     = 202)

*     Error codes and messages.
      INTEGER   SPCERR_BAD_SPEC, SPCERR_BAD_SPEC_PARAMS, SPCERR_BAD_X,
     :          SPCERR_NULL_POINTER, SPCERR_SUCCESS

      PARAMETER (SPCERR_SUCCESS         = 0)
      PARAMETER (SPCERR_NULL_POINTER    = 1)
      PARAMETER (SPCERR_BAD_SPEC_PARAMS = 2)
      PARAMETER (SPCERR_BAD_X           = 3)
      PARAMETER (SPCERR_BAD_SPEC        = 4)

      CHARACTER SPC_ERRMSG(0:4)*80
      COMMON /SPC_DATA/ SPC_ERRMSG
