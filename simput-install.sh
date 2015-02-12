#!/bin/sh

# install for SIMPUT - sh version
#
# This script assumes that the SIMPUT environment variable is set
# and uses nifty tricks from HEADAS' setup scripts
#
# Author: Joern Wilms, joern.wilms@sternwarte.uni-erlangen.de
#

if [ "X${SIMPUT}" == X ];  then
  echo "simput-install.csh: ERROR -- set SIMPUT before sourcing simput-install.csh"
  exit 1
fi

if ! [ -d ${SIMPUT} ]; then
    echo "Directory ${SIMPUT} does not exist"
    exit 2
fi

simput_bin=${SIMPUT}/bin
PATH=${simput_bin}:${PATH}

#
# setup parameter files
#
if [ "X${PFILES}" == X ]; then
    mkdir -p ${HOME}/pfiles
    PFILES="${HOME}/pfiles;${SIMPUT}/share/simput/pfiles"
else
    PFILES="${PFILES}:${SIMPUT}/share/simput/pfiles"
fi

export DYLD_LIBRARY_PATH LD_LIBRARY_PATH PATH PFILES

SIMPUT_LIB=${SIMPUT}/lib
if [ "x$LD_LIBRARY_PATH" = x ]; then
  LD_LIBRARY_PATH="$SIMPUT_LIB"
else
  LD_LIBRARY_PATH=`echo ":$LD_LIBRARY_PATH:" | sed "s%:$SIMPUT_LIB:%:%g" | sed "s%::*$%%"`
  LD_LIBRARY_PATH="$SIMPUT_LIB$LD_LIBRARY_PATH"
fi

build_os=`uname`
if ["$build_os" == "Darwin"]; then
    if [ "x$DYLD_LIBRARY_PATH" = x ]; then
	DYLD_LIBRARY_PATH="$SIMPUT_LIB"
    else
	DYLD_LIBRARY_PATH=`echo ":$DYLD_LIBRARY_PATH:" | sed "s%:$SIMPUT_LIB:%:%g" | sed "s%::*$%%"`
	DYLD_LIBRARY_PATH="$SIMPUT_LIB$DYLD_LIBRARY_PATH"
    fi
fi
