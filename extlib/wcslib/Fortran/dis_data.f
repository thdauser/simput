*=======================================================================
*
* WCSLIB 7.7 - an implementation of the FITS WCS standard.
* Copyright (C) 1995-2021, Mark Calabretta
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
* Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
* http://www.atnf.csiro.au/people/Mark.Calabretta
* $Id: dis_data.f,v 7.7 2021/07/12 06:36:49 mcalabre Exp $
*=======================================================================

      BLOCK DATA DIS_BLOCK_DATA

      CHARACTER DIS_ERRMSG(0:5)*80

      COMMON /DIS_DATA/ DIS_ERRMSG

      DATA DIS_ERRMSG /
     :  'Success',
     :  'Null disprm pointer passed',
     :  'Memory allocation failed',
     :  'Invalid parameter value',
     :  'Distort error',
     :  'De-distort error'/

      END
