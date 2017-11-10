#!/bin/bash

#################################################################
# Correcting the dependencies in SIMPUT libraries
# /!\ MAKE SURE $SIMPUT IS DEFINED (see www.sternwarte.uni-erlangen.de/research/sixte/simulation.php)
# Use after SIMPUT install and before SIXTE install
# by calling the following command: . simput_lib_correct.sh
# Author E. Cucchetti - IRAP
# If problems remains, contact the sixte support list
#################################################################

SIMPUT_PATH=$1

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/bin/simputfile
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/bin/simputfile

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputimg
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputimg

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputlc
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputlc

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputmerge
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputmerge

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputmulticell
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputmulticell

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputmultispec
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputmultispec

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputpsd
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputpsd

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputrotate
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputrotate

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputspec
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputspec

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/bin/simputsrc
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/bin/simputsrc

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib  $SIMPUT_PATH/bin/simputverify
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib  $SIMPUT_PATH/bin/simputverify

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/lib/libsimput.2.dylib
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libsimput.2.dylib

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/lib/libhdsp.1.dylib
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libhdsp.1.dylib

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/lib/libhdutils.1.dylib
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libhdutils.1.dylib

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/lib/libape.1.dylib
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libape.1.dylib

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/lib/libhdinit.1.dylib
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libhdinit.1.dylib

install_name_tool -change libcfitsio.2.dylib $SIMPUT_PATH/lib/libcfitsio.2.dylib $SIMPUT_PATH/lib/libhdio.1.dylib
install_name_tool -change libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libwcs.4.25.1.dylib $SIMPUT_PATH/lib/libhdio.1.dylib
