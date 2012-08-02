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
* $Id: fitshdr.inc,v 4.13.1.1 2012/03/14 07:40:38 cal103 Exp cal103 $
*=======================================================================

*     Functions.
      EXTERNAL  FITSHDR, FREEKEYS, KEYGET, KEYIDGET, KEYIDGTC, KEYIDGTI,
     :          KEYIDPTC, KEYIDPUT
      INTEGER   FITSHDR, FREEKEYS, KEYGET, KEYIDGET, KEYIDGTC, KEYIDGTI,
     :          KEYIDPTC, KEYIDPUT

*     Length of FITSKEY and FITSKEYID data structures (INTEGER arrays)
*     on 64-bit machines.  These are the same for 32-bit machines.
      INTEGER   KEYLEN, KEYIDLEN
      PARAMETER (KEYLEN = 48)
      PARAMETER (KEYIDLEN = 6)

*     Codes for KEYID data structure elements used by KEYIDPUT/GET.
      INTEGER   KEYID_COUNT, KEYID_IDX, KEYID_NAME

      PARAMETER (KEYID_NAME    = 100)
      PARAMETER (KEYID_COUNT   = 101)
      PARAMETER (KEYID_IDX     = 102)

*     Codes for KEY data structure elements used by KEYGET.
      INTEGER   KEY_COMMENT, KEY_KEYID, KEY_KEYNO, KEY_KEYVALUE,
     :          KEY_KEYWORD, KEY_STATUS, KEY_TYPE, KEY_ULEN

      PARAMETER (KEY_KEYNO     = 200)
      PARAMETER (KEY_KEYID     = 201)
      PARAMETER (KEY_STATUS    = 202)
      PARAMETER (KEY_KEYWORD   = 203)
      PARAMETER (KEY_TYPE      = 204)
      PARAMETER (KEY_KEYVALUE  = 205)
      PARAMETER (KEY_ULEN      = 206)
      PARAMETER (KEY_COMMENT   = 207)

      CHARACTER FITSHDR_ERRMSG(0:2)*80
      COMMON /FITSHDR_DATA/ FITSHDR_ERRMSG
