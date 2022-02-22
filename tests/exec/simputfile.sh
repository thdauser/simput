#!/bin/bash

. setup/setup.sh

simput1="athenacrab_test.fits"
simput2="webexample_test.fits"

simputfile \
  simput=$outdir/$simput1 \
  RA=0.0 Dec=0.0 \
  Emin=2. Emax=10. \
  plPhoIndex=2.1 NH=0.4 plFlux=2.052123e-10 \
  clobber=yes

simputfile \
  Simput=$outdir/$simput2 \
  RA=40.2 Dec=12.8 \
  XSPECFile=$indir/example_spectrum.xcm \
  LCFile=$indir/example_lightcurve.dat \
  MJDREF=50800.0 \
  Emin=0.5 Emax=10.0 \
  srcFlux=2.3e-12 \
  clobber=yes

