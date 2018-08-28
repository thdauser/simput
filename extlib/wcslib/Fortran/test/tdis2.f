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
* $Id: tdis2.f,v 5.19.1.1 2018/07/26 15:41:42 mcalabre Exp mcalabre $
*=======================================================================

      PROGRAM TDIS2
*-----------------------------------------------------------------------
*
* TDIS2 tests the WCSLIB Fortran routines for manipulating the contents
* of a disprm struct.  Two disprm structs are used.  One, residing in a
* linprm struct, is referenced by "address" (DISP), and the other (DISL)
* resides locally in an INTEGER array.
*
* DISP and WCSP, which are meant to hold addresses, are declared as
* INTEGER arrays of length 2 to accomodate 64-bit machines for which
* sizeof(void *) = 2*sizeof(int).
*
* Test input comes from TPV7.fits.
*-----------------------------------------------------------------------
      INTEGER   DPMAX
      PARAMETER (DPMAX = 256)

      LOGICAL   GOTEND
      INTEGER   DISP(2), IBLOCK, IKEYRC, J, JK, JN, K, NAXIS, NDP,
     :          NKEYRC, NREJECT, NWCS, STATUS, WCSP(2)
      DOUBLE PRECISION MAXDIS(2), TOTDIS
      CHARACTER DTYPE(2)*72, KEYREC*80, HEADER*288000, INFILE*9

*     On some systems, such as Sun Sparc, the structs MUST be aligned
*     on a double precision boundary, done here using an equivalence.
*     Failure to do this may result in mysterious "bus errors".
      INCLUDE 'wcserr.inc'
      INCLUDE 'wcshdr.inc'
      INCLUDE 'wcs.inc'
      INCLUDE 'lin.inc'
      INCLUDE 'dis.inc'
      INTEGER   DP(DPMAX*DPLEN), WCS(WCSLEN), LIN(LINLEN), DISL(DISLEN)
      DOUBLE PRECISION DUMMY1, DUMMY2, DUMMY3, DUMMY4
      EQUIVALENCE (DUMMY1, DP)
      EQUIVALENCE (DUMMY2, WCS)
      EQUIVALENCE (DUMMY3, LIN)
      EQUIVALENCE (DUMMY4, DISL)

      DATA INFILE /'TPV7.fits'/
*-----------------------------------------------------------------------
      STATUS = WCSERR_ENABLE(1)

      WRITE (*, 10)
 10   FORMAT (
     :  'Testing WCSLIB Fortran disprm access routines (tdis2.f)',/,
     :  '-------------------------------------------------------')

*     Open the FITS WCS test header for formatted, direct I/O.
      OPEN (UNIT=1, FILE=INFILE, FORM='FORMATTED', ACCESS='DIRECT',
     :      RECL=80, IOSTAT=STATUS)
      IF (STATUS.NE.0) THEN
        WRITE (*, 20) STATUS, INFILE
 20     FORMAT ('ERROR',I3,' opening ',A)
        GO TO 999
      END IF

*     Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
      K = 1
      NKEYRC = 0
      GOTEND = .FALSE.
      DO 50 IBLOCK = 0, 100
        DO 40 IKEYRC = 1, 36
          READ (1, '(A80)', REC=36*IBLOCK+IKEYRC, IOSTAT=STATUS) KEYREC
          IF (STATUS.NE.0) THEN
            WRITE (*, 30) STATUS
 30         FORMAT ('ERROR',I5,' reading header.')
            GO TO 999
          END IF

          IF (KEYREC(:8).EQ.'        ') GO TO 40
          IF (KEYREC(:8).EQ.'COMMENT ') GO TO 40
          IF (KEYREC(:8).EQ.'HISTORY ') GO TO 40

          HEADER(K:) = KEYREC
          K = K + 80
          NKEYRC = NKEYRC + 1

          IF (KEYREC(:10).EQ.'END       ') THEN
*           An END keyrecord was read, but read the rest of the block.
            GOTEND = .TRUE.
          END IF
 40     CONTINUE

        IF (GOTEND) GO TO 60
 50   CONTINUE

 60   CLOSE (UNIT=1)


*     Parse the header.
      CALL FLUSH(6)
      STATUS = WCSPIH (HEADER, NKEYRC, WCSHDR_none, 2, NREJECT, NWCS,
     :                 WCSP)
      IF (STATUS.NE.0) GO TO 999

*     Copy into our WCSPRM struct.
      STATUS = WCSVCOPY (WCSP, 0, WCS)
      IF (STATUS.NE.0) GO TO 999

*     Translate the TPV "projection" into a sequent distortion.
      STATUS = WCSSET(WCS)
      IF (STATUS.NE.0) THEN
        STATUS = WCSPERR(WCS, CHAR(0))
        GO TO 999
      END IF

*     Extract the linprm struct.
      STATUS = WCSGET (WCS, WCS_LIN, LIN)
      IF (STATUS.NE.0) THEN
        STATUS = WCSPERR(WCS, CHAR(0))
        GO TO 999
      END IF

*     Extract a "pointer" to the disprm struct.
      STATUS = LINGET(LIN, LIN_DISSEQ, DISP)
      IF (STATUS.NE.0) THEN
        STATUS = LINPERR(LIN, CHAR(0))
        GO TO 999
      END IF

*     Add a little texture.
      STATUS = DISPTD (1, DISP, DIS_MAXDIS, 2.0010D0, 1, 0)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISPTD (1, DISP, DIS_MAXDIS, 2.0020D0, 2, 0)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISPTD (1, DISP, DIS_TOTDIS, 3.0000D0, 0, 0)
      IF (STATUS.NE.0) GO TO 997

*     Print its contents.
      CALL FLUSH(6)
      STATUS = DISPRT(1, DISP)
      IF (STATUS.NE.0) GO TO 997

      WRITE (*, 70)
 70   FORMAT (/,'------------------------------------',
     :          '------------------------------------')


*   Copy it the long way first to test the various PUT and GET routines.
*     Start by reading everything from the first.
      STATUS = DISGET (1, DISP, DIS_NAXIS,  NAXIS)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISGTC (1, DISP, DIS_DTYPE,  DTYPE)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISGTI (1, DISP, DIS_NDP,    NDP)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISGTD (1, DISP, DIS_DP,     DP)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISGTD (1, DISP, DIS_MAXDIS, MAXDIS)
      IF (STATUS.NE.0) GO TO 997
      STATUS = DISGTD (1, DISP, DIS_TOTDIS, TOTDIS)
      IF (STATUS.NE.0) GO TO 997

*     Initialize the destination.
      STATUS = DISPUT (0, DISL, DIS_FLAG, -1, 0, 0)
      IF (STATUS.NE.0) GO TO 998
      STATUS = DISINI (0, NAXIS, DISL)
      IF (STATUS.NE.0) GO TO 998

*     Copy to the destination.
      JK = 0
      JN = 0
      DO 110 J = 1, NAXIS
        STATUS = DISPTD (0, DISL, DIS_DTYPE, DTYPE, J, 0)
        IF (STATUS.NE.0) GO TO 998
 110  CONTINUE

      DO 120 K = 1, NDP*DPLEN, DPLEN
        STATUS = DISPTI (0, DISL, DIS_DP, DP(K), 0, 0)
 120  CONTINUE

      DO 130 J = 1, NAXIS
        STATUS = DISPTI (0, DISL, DIS_MAXDIS, MAXDIS(J), J, 0)
        IF (STATUS.NE.0) GO TO 998
 130  CONTINUE

      STATUS = DISPTI (0, DISL, DIS_TOTDIS, TOTDIS, J, 0)
      IF (STATUS.NE.0) GO TO 998

      STATUS = DISSET (0, DISL)
      IF (STATUS.NE.0) GO TO 998

*     Print the second.
      CALL FLUSH(6)
      STATUS = DISPRT(0, DISL)
      IF (STATUS.NE.0) GO TO 998

      WRITE (*, 70)


*   Now copy it the fast way.
*     Start with a clean slate.
      STATUS = DISINI (0, NAXIS, DISL)
      IF (STATUS.NE.0) GO TO 998

*     Here the source is a "pointer", and the destination a local array.
      STATUS = DISCPY (1, DISP, DISL)
      IF (STATUS.NE.0) GO TO 998

*     A reset is required before printing.
      STATUS = DISSET (0, DISL)
      IF (STATUS.NE.0) GO TO 998

*     Print the copy.
      CALL FLUSH(6)
      STATUS = DISPRT(0, DISL)
      IF (STATUS.NE.0) GO TO 998


*   Cleanup.
 997  IF (STATUS.NE.0) THEN
        STATUS = DISPERR(1, DISP, CHAR(0))
      END IF

      GO TO 999

 998  IF (STATUS.NE.0) THEN
        STATUS = DISPERR(0, DISL, CHAR(0))
      END IF

      STATUS = DISFREE (0, DISL)

*     Free memory allocated by WCSPIH, WCSINI, etc.
*     Should free DISP as well.
 999  STATUS = WCSVFREE (NWCS, WCSP)

      END
