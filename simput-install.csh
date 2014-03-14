#!/bin/csh

# install for SIMPUT - csh version
#
# This script assumes that the SIMPUT environment variable is set
# and uses nifty tricks from HEADAS' setup scripts
#
# Author: Joern Wilms, joern.wilms@sternwarte.uni-erlangen.de
#

if (${?SIMPUT} == 0) then
  echo "simput-install.csh: ERROR -- set SIMPUT before sourcing simput-install.csh"
  exit 1
endif

if (! -d ${SIMPUT}) then
    echo "Directory ${SIMPUT} does not exist"
    exit 2
endif

set SIMPUT_BIN = ${SIMPUT}/bin
setenv PATH ${SIMPUT_BIN}:${PATH}

#
# setup parameter files
#
if (${?PFILES} == 0) then
    mkdir -p ${HOME}/pfiles
    setenv PFILES "${HOME}/pfiles;${SIMPUT}/share/simput/pfiles"
else
    setenv PFILES "${PFILES}:${SIMPUT}/share/simput/pfiles"
endif
