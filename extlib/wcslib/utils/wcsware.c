/*============================================================================

  WCSLIB 5.19 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2018, Mark Calabretta

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

  Direct correspondence concerning WCSLIB to mark@calabretta.id.au

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcsware.c,v 5.19.1.1 2018/07/26 15:41:42 mcalabre Exp mcalabre $
*=============================================================================
* wcsware extracts the WCS keywords for an image from the specified FITS file,
* constructs wcsprm structs for each coordinate representation found and
* performs a variety of operations using them.
*---------------------------------------------------------------------------*/

char usage[] =
"Usage: wcsware [OPTION]... [FITSFILE]\n"
"\n"
"wcsware extracts the WCS keywords for an image from the specified FITS\n"
"file, constructs wcsprm structs for each coordinate representation found\n"
"and performs a variety of operations using them.  It features a \"lint\"\n"
"capability as one of its more useful operations.\n\n"
"By default, all known extensions to the FITS WCS standard are allowed,\n"
"including deprecated usage.  However, in \"lint\" mode, strict conformance\n"
"to the standard is enforced.\n\n"
"The FITS file may be specified according to the syntax understood by\n"
"cfitsio, for example \"file.fits.gz+1\" refers to the first extension of\n"
"a gzip'd FITS file.  Use \"-\" or omit the file name for input from stdin.\n"
"\n"
"Options:\n"
"  -a<alt>      Specify an alternate coordinate representation (ignored if\n"
"               there is only one).  Can also be specified as a 0-relative\n"
"               index in the range 0 to 26, where alternates are sequenced\n"
"               alphabetically following the primary representation.\n"
"  -b           Use wcsbth() for primary image headers, normally wcspih()\n"
"               is used.  (Implies -i.)\n"
"  -f           Apply wcsfix() to the header.\n"
"  -h<hdu>      Move to HDU number (1-relative) which is expected to\n"
"               contain an image array.  Overrides cfitsio extended\n"
"               filename syntax.  Also useful for input from stdin.\n"
"  -i           Allow image header WCS keywords in binary table headers.\n"
"  -l           Validate (lint) the WCS keyrecords in the specified FITS\n"
"               header for conformance to the WCS standard.  (Implies -s.)\n"
"  -o           Use wcshdo() to translate the wcsprm struct into a FITS\n"
"               header and print it.\n"
"  -p           Print the struct(s) using wcsprt() (default operation).\n"
"  -P           Same as -p but don't print a default struct.\n"
"  -r           Require strict adherence to the FITS WCS standard, though\n"
"               allowing the deprecated AIPS-convention keywords, CROTAn,\n"
"               EPOCH, and VELREF, and also some other deprecated usage.\n"
"               (Must follow -l if relaxed linting is required.)\n"
"  -s           Require strict adherence to the FITS WCS standard,\n"
"               disallowing all deprecated features.\n"
"  -t           Terse (with -l), report rejected WCS keyrecords only.\n"
"  -v           Verbose (with -l), report recognised WCS keyrecords as well.\n"
"  -w           Convert world coordinates, obtained from stdin, to pixel\n"
"               coordinates using wcss2p().\n"
"  -x           Convert pixel coordinates, obtained from stdin, to world\n"
"               coordinates using wcsp2s().\n";

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <fitsio.h>

#include <wcs.h>
#include <wcshdr.h>
#include <wcsfix.h>
#include <wcsprintf.h>
#include <getwcstab.h>

int main(int argc, char **argv)

{
  char *alt = 0x0, *header, *hptr, *infile, wcsname[72];
  int  allimg, alts[27], c, ctrl, dobth, dofix, dohdr, dolint, doprt, dopix,
       dowrld, hdunum, hdutype, i, ialt, j, k, keysel, nelem, nkeyrec,
       nreject, nwcs, ok, relax, *stat = 0x0, status, strict;
  double *imgcrd = 0x0, phi, *pixcrd = 0x0, theta, *world = 0x0;
  struct wcsprm *wcs, *wcsi = 0x0;
  fitsfile *fptr;

  /* Initialise control variables. */
  allimg =  0;
  ctrl   = -3;
  dobth  =  0;
  dofix  =  0;
  dohdr  =  0;
  dolint =  0;
  dopix  =  0;
  doprt  =  0;
  dowrld =  0;
  hdunum =  0;
  relax  = WCSHDR_all;
  strict =  0;

  /* Parse options. */
  for (i = 1; i < argc && argv[i][0] == '-'; i++) {
    if (!argv[i][1]) break;

    switch (argv[i][1]) {
    case 'a':
      /* Select an alternate WCS. */
      alt = argv[i]+2;
      break;

    case 'b':
      /* Use wcsbth() for image headers. */
      dobth  = 1;
      allimg = 1;
      break;

    case 'f':
      /* Apply wcsfix(). */
      dofix  = 1;
      break;

    case 'h':
      /* Move to HDU number. */
      hdunum = atoi(argv[i]+2);
      if (hdunum < 0) hdunum = 0;
      break;

    case 'i':
      /* Allow image header keywords in bintables. */
      allimg = 1;
      break;

    case 'l':
      /* Lint. */
      strict = 1;
      relax  = WCSHDR_strict;
      ctrl   = 3;
      dolint = 1;
      break;

    case 'o':
      /* Print header using wcshdo(). */
      dohdr = 1;
      break;

    case 'p':
      /* Print structs. */
      doprt  = 1;
      break;

    case 'P':
      /* Print non-default structs. */
      doprt  = -1;
      break;

    case 'r':
      /* Strict, but allow some deprecated usage. */
      strict = 0;
      relax  = WCSHDR_reject;
      break;

    case 's':
      /* Strict, really strict. */
      strict = 1;
      relax  = WCSHDR_strict;
      break;

    case 't':
      /* Terse. */
      ctrl   = 2;
      break;

    case 'v':
      /* Verbose. */
      ctrl   = 4;
      break;

    case 'w':
      dowrld = 1;
      break;

    case 'x':
      dopix  = 1;
      break;

    default:
      wcsfprintf(stderr, "%s", usage);
      return 1;
    }
  }

  if (allimg) {
    relax |= WCSHDR_ALLIMG;
  }

  if (i < argc) {
    infile = argv[i++];

    if (i < argc) {
      wcsfprintf(stderr, "%s", usage);
      return 1;
    }
  } else {
    infile = "-";
  }

  if (!(dolint || dohdr || doprt || dopix || dowrld)) doprt = 1;


  /* Open the FITS file and move to the required HDU. */
  status = 0;
  if (fits_open_file(&fptr, infile, READONLY, &status)) goto fitserr;
  if (hdunum) {
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) goto fitserr;
  } else {
    if (fits_get_hdu_type(fptr, &hdutype, &status)) goto fitserr;
  }

  /* Read in the FITS header, excluding COMMENT and HISTORY keyrecords. */
  if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
    goto fitserr;
  }


  /* Interpret the WCS keywords. */
  if ((hdutype == BINARY_TBL) ||
     ((hdutype == IMAGE_HDU) && dobth)) {
    keysel = (hdutype == BINARY_TBL) ? 0 : WCSHDR_IMGHEAD;
    if ((status = wcsbth(header, nkeyrec, relax, ctrl, keysel, 0x0, &nreject,
                    &nwcs, &wcs))) {
      wcsfprintf(stderr, "wcsbth ERROR %d: %s.\n", status,
                 wcshdr_errmsg[status]);
      return 1;
    }
  } else if (hdutype == IMAGE_HDU) {
    if ((status = wcspih(header, nkeyrec, relax, ctrl, &nreject,
                    &nwcs, &wcs))) {
      wcsfprintf(stderr, "wcspih ERROR %d: %s.\n", status,
                 wcshdr_errmsg[status]);
      return 1;
    }
  } else {
    wcsfprintf(stderr, "wcsware: Invalid FITS extension type.\n");
    return 1;
  }

  free(header);

  if (nreject) {
    if (ctrl <= 3) {
      wcsprintf("\n%d WCS keyrecords were rejected.\n", nreject);
    }
    wcsprintf("\nThe rejected keyrecords do not conform%s to the FITS WCS "
              "standard.\n", strict?" strictly":"");
  } else if (nwcs == 0) {
    if (2 < ctrl) wcsprintf("\n");
    wcsprintf("No world coordinate systems found.\n");
    fits_close_file(fptr, &status);
    return 0;

  } else if (2 < ctrl) {
    wcsprintf("\nNo invalid WCS keyrecords were found.\n");
  }


  /* Sort out alternates. */
  if (alt) {
    i = 0;

    if ('0' <= *alt && *alt <= '9') {
      if ((i = atoi(alt)) > nwcs-1) {
        wcsfprintf(stderr, "WARNING, no alternate coordinate "
          "representation \"%s\".\n", alt);
        return 1;
      }

    } else {
      wcsidx(nwcs, &wcs, alts);

      ialt = toupper(*alt);
      if (strlen(alt) > 1) {
        wcsfprintf(stderr, "WARNING, alternate specifier \"%s\" is "
          "invalid.\n", alt);
        return 1;

      } else if (*alt == ' ') {
        if (alts[0] == -1) {
          wcsfprintf(stderr, "WARNING, no primary coordinate "
            "representation.\n");
          return 1;
        }

      } else if (ialt < 'A' || ialt > 'Z') {
        wcsfprintf(stderr, "WARNING, alternate specifier \"%s\" is "
          "invalid.\n", alt);
        return 1;

      } else {
        if ((i = alts[ialt - 'A' + 1]) == -1) {
          wcsfprintf(stderr, "WARNING, no alternate coordinate "
            "representation \"%s\".\n", alt);
          return 1;
        }
      }
    }

    wcsi = wcs + i;
  }


  /* Initialize and possibly print the structs. */
  wcserr_enable(1);
  wcsprintf_set(stdout);

  for (i = 0; i < nwcs; i++) {
    if (wcsi) {
      if ((wcs+i) != wcsi) continue;
    } else if (i) {
      wcsprintf("\nType <CR> for next: ");
      fgetc(stdin);
      wcsprintf("\n");
    }

    /* Read -TAB arrays from the binary table extension (if necessary). */
    if (fits_read_wcstab(fptr, wcs[i].nwtb, (wtbarr *)wcs[i].wtb,
                         &status)) {
      goto fitserr;
    }

    /* Translate non-standard WCS keyvalues? */
    if (dofix) {
      stat = malloc(NWCSFIX * sizeof(int));
      if ((status = wcsfix(7, 0, wcs+i, stat))) {
        for (j = 0; j < NWCSFIX; j++) {
          if (stat[j] > 0) {
             wcsfprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
                        wcsfix_errmsg[stat[j]]);
          }
        }

        return 1;
      }
    }

    if ((status = wcsset(wcs+i))) {
      wcsperr(wcs+i, "");
      continue;
    }

    /* Get WCSNAME out of the wcsprm struct. */
    strcpy(wcsname, wcs[i].wcsname);

    /* Print the struct as a FITS header. */
    if (dohdr) {
      if (wcshdo(WCSHDO_all, wcs+i, &nkeyrec, &header)) {
        wcsperr(wcs+i, "");
        continue;
      }

      hptr = header;
      wcsprintf("\n");
      for (k = 0; k < nkeyrec; k++, hptr += 80) {
        wcsprintf("%.80s\n", hptr);
      }

      free(header);
    }

    /* Print the struct. */
    if (doprt) {
      if (doprt == 1 || strcmp(wcsname, "DEFAULTS")) {
        if (strlen(wcsname)) {
          wcsprintf("\n%s\n", wcsname);
        }

        wcsprintf("\n");
        wcsprt(wcs+i);
      }
    }

    /* Transform coordinates? */
    if (dopix || dowrld) {
      if (strlen(wcsname)) {
        wcsprintf("\n%s\n", wcsname);
      }

      nelem = wcs[i].naxis;
      world  = realloc(world,  nelem * sizeof(double));
      imgcrd = realloc(imgcrd, nelem * sizeof(double));
      pixcrd = realloc(pixcrd, nelem * sizeof(double));
      stat   = realloc(stat,   nelem * sizeof(int));

      if (dopix) {
        /* Transform pixel coordinates. */
        while (1) {
          wcsprintf("\nEnter %d pixel coordinate element%s: ", nelem,
            (nelem==1)?"":"s");
          c = fgetc(stdin);
          if (c == EOF || c == '\n') {
            if (c == EOF) wcsprintf("\n");
            break;
          }
          ungetc(c, stdin);

          if ((ok = (scanf("%lf", pixcrd) == 1))) {
            for (j = 1; j < nelem; j++) {
              if (scanf("%*[ ,]%lf", pixcrd+j) != 1) {
                ok = 0;
                break;
              }
            }
          }

          while (fgetc(stdin) != '\n');

          if (!ok) {
            wcsprintf("Input error, please try again.\n");
            continue;
          }

          wcsprintf("Pixel: ");
          for (j = 0; j < nelem; j++) {
            wcsprintf("%s%14.9g", j?", ":"", pixcrd[j]);
          }

          if ((status = wcsp2s(wcs+i, 1, nelem, pixcrd, imgcrd, &phi, &theta,
                               world, stat))) {
            wcsperr(wcs+i, "");

          } else {
            wcsprintf("\nImage: ");
            for (j = 0; j < nelem; j++) {
              if (j == wcs[i].lng || j == wcs[i].lat) {
                /* Print angles in fixed format. */
                wcsprintf("%s%14.6f", j?", ":"", imgcrd[j]);
              } else {
                wcsprintf("%s%14.9g", j?", ":"", imgcrd[j]);
              }
            }

            wcsprintf("\nWorld: ");
            for (j = 0; j < nelem; j++) {
              if (j == wcs[i].lng || j == wcs[i].lat) {
                /* Print angles in fixed format. */
                wcsprintf("%s%14.6f", j?", ":"", world[j]);
              } else {
                wcsprintf("%s%14.9g", j?", ":"", world[j]);
              }
            }
            wcsprintf("\n");
          }
        }
      }


      if (dowrld) {
        /* Transform world coordinates. */
        while (1) {
          wcsprintf("\nEnter %d world coordinate element%s: ", nelem,
            (nelem==1)?"":"s");
          c = fgetc(stdin);
          if (c == EOF || c == '\n') {
            if (c == EOF) wcsprintf("\n");
            break;
          }
          ungetc(c, stdin);

          if ((ok = (scanf("%lf", world) == 1))) {
            for (j = 1; j < nelem; j++) {
              if (scanf("%*[ ,]%lf", world+j) != 1) {
                ok = 0;
                break;
              }
            }
          }

          while (fgetc(stdin) != '\n');

          if (!ok) {
            wcsprintf("Input error, please try again.\n");
            continue;
          }

          wcsprintf("World: ");
          for (j = 0; j < nelem; j++) {
            if (j == wcs[i].lng || j == wcs[i].lat) {
              /* Print angles in fixed format. */
              wcsprintf("%s%14.6f", j?", ":"", world[j]);
            } else {
              wcsprintf("%s%14.9g", j?", ":"", world[j]);
            }
          }

          if ((status = wcss2p(wcs+i, 1, nelem, world, &phi, &theta, imgcrd,
                               pixcrd, stat))) {
            wcsperr(wcs+i, "");

          } else {
            wcsprintf("\nImage: ");
            for (j = 0; j < nelem; j++) {
              if (j == wcs[i].lng || j == wcs[i].lat) {
                /* Print angles in fixed format. */
                wcsprintf("%s%14.6f", j?", ":"", imgcrd[j]);
              } else {
                wcsprintf("%s%14.9g", j?", ":"", imgcrd[j]);
              }
            }

            wcsprintf("\nPixel: ");
            for (j = 0; j < nelem; j++) {
              wcsprintf("%s%14.9g", j?", ":"", pixcrd[j]);
            }
            wcsprintf("\n");
          }
        }
      }
    }
  }

  fits_close_file(fptr, &status);

  /* Defeat spurious reporting of memory leaks. */
  wcsvfree(&nwcs, &wcs);
  free(world);
  free(imgcrd);
  free(pixcrd);
  free(stat);

  return 0;

fitserr:
  fits_report_error(stderr, status);
  fits_close_file(fptr, &status);
  return 1;
}
